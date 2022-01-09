using JuMP, MosekTools, Combinatorics, FileIO


function read_input(filename)
  cmd = pipeline(`python3 -m kidney_solver.utils.convert`, stdin=filename)
  lines = readlines(cmd)
  sep = findall(x->x!="-1\t-1\t-1", lines)
  lines = lines[sep]
  lines = map(x->split(x, "\t"), lines)
  lines = [map(x->parse(Int,x), line) for line in lines]
  num_P, num_A = popfirst!(lines)

  cycle_arcs = [tuple(lines[i][1:2]...) for i=1:num_A]
  num_N, num_E = splice!(lines, num_A+1)
  ndd_edges = [tuple(lines[i+num_A][1:2]...) for i=1:num_E]
  ndd_edges = map(((x,y),)->(x+num_P,y), ndd_edges)
  G = (Set(i for i=0:num_P+num_N-1), union(Set(cycle_arcs), Set(ndd_edges))) 
  return G, num_P
end


function read_data(filename)
  lines = readlines(filename)
  lines = map(x->split(x, ","), lines)
  col_idx = findfirst(x->x=="%Pra", popfirst!(lines))
  pra_keys = map(x->parse(Int, x)-1, [line[1] for line in lines])
  pra_vals = map(x->parse(Float64, x), [line[col_idx] for line in lines])
  pra_dict = Dict(zip(pra_keys,pra_vals))
end


function copy_shortest_paths(G)
  (V, A) = G
  n = length(V)
  d = Dict((i,j,l)=>Inf for l=0:n-2 for i=l:n-1,j=l:n-1)
  v_type = Int64; e_type = Tuple{v_type, v_type}
  Γ = Dict{Int64, Tuple{Set{v_type},Set{e_type}}}()

  for l=n-1:-1:0
    V_l = Set(i for i=l:n-1)
    A_l = Set((i,j) for (i,j) in A if i in V_l && j in V_l)

    for (i,j) in A_l
      d[i,j,l] = 1
    end

    for i in V_l
      d[i,i,l] = 0
    end

    # find the shortest path between all pairs (i,j)
    # Known as Floyd-Warshall algorithm
    for k=l:n-1
      for i=l:n-1
        for j=l:n-1
          if d[i,j,l] > d[i,k,l] + d[k,j,l]
            d[i,j,l] = d[i,k,l] + d[k,j,l]
          end
        end
      end
    end

    Γ[l] = (V_l, A_l)
  end
  return Γ, d
end


function get_positions(G, d, L)
  (V,_) = G
  n = length(V)
  K = Dict{NTuple{3,Int64}, Set{Int64}}()

  for l=0:n-2
    for i=l:n-1
      for j=l:n-1
        K[i,j,l] = Set(k for k=1:L if d[l,i,l] < k && d[j,l,l] <= L-k)
      end
    end
  end

  return K
end


function columngen_params(model, unique, num_P)
  @variable(model, z)
  @variable(model, abs_var[i=0:num_P-1])
  @constraint(model, num_P * z == objective_function(model))
  @constraint(model, abs_cons[i=0:num_P-1], abs_var[i] == constraint_object(unique[i]).func - z)
  optimize!(model)

  if termination_status(model) != MOI.OPTIMAL && (termination_status(model) != MOI.TIME_LIMIT || !has_values(model)) 
    error("The model was not solved correctly")
  end

  feasible = Set{Int64}()
  for i=0:num_P-1
    temp = @constraint(model, constraint_object(unique[i]).func == 1)
    optimize!(model)
    delete(model, temp)
    if termination_status(model) == MOI.INFEASIBLE
      continue
    else
      push!(feasible, i)
    end
  end

  delete(model, z)
  unregister(model, :z)
  for i=0:num_P-1
    delete(model, abs_var[i])
  end
  vcount = Dict([i=>round(Int, value(unique[i])) for i in feasible])
  return vcount 
end


function sub_problem(model, obj_expr, x)
  @objective(model, Min, obj_expr)
  optimize!(model)

  if termination_status(model) == MOI.OPTIMAL
    ζ = objective_value(model)
    if ζ > -1e-12
      return false 
    else
      return true 
    end
  elseif termination_status(model) == MOI.TIME_LIMIT
    error("The time limit was reached for the subproblem")
  else
    error("The subproblem could not be solved correctly")
  end 
end


function master_problem(num_P, Γ, K, L, pra_dict)
  # start time
  stime = time() 
  # Set of sensitized patients 
  if pra_dict != nothing
    sensitized = Set(i for (i,pra) in pra_dict if pra >= 0.8)
  else
    sensitized = Set()
  end
  n = length(Γ)
  
  # Model i.e Subproblem
  model = Model(Mosek.Optimizer)
  set_optimizer_attribute(model, "MSK_IPAR_LOG", 0)

  @variable(model, x[l=0:n-2,(i,j)=Γ[l][2],K[i,j,l]], Bin)

  @constraint(model, unique[i=0:n-1], sum(x[l,(j,i),k] for l=0:i, (j,v) in Γ[l][2], k in K[j,v,l] if v == i) <= 1)
  @constraint(model, flow[l=0:n-2,i=l+1:n-1,k=1:L-1], sum(x[l,(j,i),k] for (j,v) in Γ[l][2] if v == i && k in K[j,i,l]) == sum(x[l,(i,j),k+1] for (v,j) in Γ[l][2] if v == i && (k+1) in K[i,j,l]))
  f_1 = sum(constraint_object(unique[i]).func for i=0:num_P-1)
  @objective(model, Max, f_1)
  
  # obtain the first solution and the nadir point
  vcount = columngen_params(model, unique, num_P)
  transplants = objective_value(model)
  P = Set{Int64}(keys(vcount))

  # delete the constraints that will not be used from the model
  for i=0:num_P-1
    if haskey(vcount, i)
      continue
    else
      delete(model, unique[i])
    end
  end

  for l=0:min(num_P-1,n-2)
    if haskey(vcount, l)
      continue
    else
      for (i,j) in Γ[l][2]
        if !isempty(K[i,j,l]) 
          delete(model, x[l,(i,j),K[i,j,l]])
        end
      end
    end
    
    for i=l+1:n-1
      for k=1:L-1
        delete(model, flow[l,i,k])
      end
    end
  end
  num_P = length(keys(vcount))

  # store solutions
  z = Dict(i=>vcount[i] for i=P)
  solutions = 1

  _unique = Dict(i=>constraint_object(unique[i]).func for i=P)
  # add solution cuts and reoptimize
  while true
    @constraint(model, sum(value(unique[i]) > 0.5 ? 1.0 - _unique[i] : _unique[i] for i=P) >= 1)
    optimize!(model)
    
    status = termination_status(model)
    if status == MOI.OPTIMAL
      for i=P
        z[i] += abs(value(_unique[i]))
      end
      solutions += 1
    elseif status == MOI.INFEASIBLE
      break
    elseif status == MOI.TIME_LIMIT
      error("The solver could not terminate in the allocated time")
    else
      error("The model could not be solved")
    end
  end
     
  z0 = sum(values(z)) / solutions / num_P
  L1 = sum(abs(z[i] / solutions - z0) for i=P)
  Minprob = minimum(z[i] for i=P) / solutions
  SumLog = sum(log(z[i]) for i=P)
  println("Number of transplants: ", transplants)
  println("L1: ", L1)
  println("Minprob: ", Minprob)
  println("SumLog: ", SumLog)
  etime = time() - stime
  stats = Dict("Utility:L1"=>L1, "Utility:Minprob"=>Minprob, "Utility:SumLog"=>SumLog, "Utility:transplants"=>transplants, "Utility:solutions"=>solutions, "Utility:time"=>etime)
  return stats
end


function main()
  path = ARGS[1]
  L = parse(Int64, ARGS[2])
  
  G, num_P = read_input(path*".wmd")
  pra_dict = isfile(path*".dat") ? read_data(path*".dat") : nothing

  (Γ, d) = copy_shortest_paths(G)
  K = get_positions(G, d, L)

  filename = split(path, '/')[end]

  stats = master_problem(num_P, Γ, K, L, pra_dict)

  if !isfile("$filename.jld2") 
    save("$filename.jld2", stats)
  else
    _stats = load("$filename.jld2")
    for key in keys(stats)
      if !haskey(_stats, key)
        _stats[key] = stats[key]
      end
    end
    save("$filename.jld2", _stats)
  end

end


main()
