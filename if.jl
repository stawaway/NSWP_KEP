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
  d = Dict((i,j,l)=>Inf for l=0:n-2 for i=l:n-1,j=l:n-1)  #Dict((i,j,l)=>Inf for i=0:n-1,j=0:n-1,l=0:n-2)
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


function get_tuples(s, n, L)
  if L==1
    return [Set(i) for i=s:n-1]
  end
  arr = Array{Set{Int64}, 1}()
  for i=s:n-L
    temp = get_tuples(i+1, n, L-1)
    for cycle in temp
      push!(cycle, i)
    end
    arr = cat(arr, temp, dims=1)
  end
  return arr
end


function is_cycle(tup, A)
  target = first(tup)
  n = length(tup)
  perms = permutations(collect(tup))
  for p in perms
    flag = true
    for (i,x) in enumerate(p)
      next = (i % n) + 1
      if !((x, p[next]) in A)
        flag = false
        break
      end
    end
    if flag
      return flag
    end
  end
  return false
end


function cycle_counts(G, L)
  (V, A) = G
  n = length(V)
  comb = []
  for k=L:-1:2
    comb = cat(comb, get_tuples(0, n, k), dims=1)
  end
  cycles = [c for c in comb if is_cycle(c, A)]
  counts = Dict(i=>count(x->i in x, cycles) for i=0:n-1)
end


#function f_1(x)
#  sum(x[idx] for idx in eachindex(x))
#end


#function f_2(x, sensitized)
#  expr = (x[l,(i,j),k] for (l,(i,j),k) in eachindex(x) if j in sensitized)
#  isempty(expr) ? GenericAffExpr{Float64, VariableRef}(0) : sum(expr)
#end


#function choose_obj(fname, x, sensitized)
#  if fname == "f_1"
#    f_1(x)
#  elseif fname == "f_2"
#    f_2(x,sensitized)
#  else
#    throw("Wrong function name for the objective. Try another one.")
#  end
#end


function columngen_params(P, Γ, K, L, pra_dict)
  if pra_dict != nothing
    sensitized = Set(i for (i,pra) in pra_dict if pra >= 0.8)
  else
    sensitized = Set()
  end
  n = length(Γ) 
  model = Model(Mosek.Optimizer)
  set_optimizer_attribute(model, "MSK_IPAR_LOG", 0)

  @variable(model, x[l=0:n-2,(i,j)=Γ[l][2],K[i,j,l]], Bin)

  @constraint(model, unique[i=0:n-1], sum(x[l,(j,i),k] for l=0:i, (j,v) in Γ[l][2], k in K[j,i,l] if v == i) <= 1)
  @constraint(model, flow[l=0:n-2,i=l+1:n-1,k=1:L-1], sum(x[l,(j,i),k] for (j,v) in Γ[l][2] if v == i && k in K[j,i,l]) == sum(x[l,(i,j),k+1] for (v,j) in Γ[l][2] if v == i && (k+1) in K[i,j,l]))

  @variable(model, z)
  @variable(model, abs_var[i=0:P-1])
  f_1 = sum(constraint_object(unique[i]).func for i=0:P-1)
  @constraint(model, P * z == f_1)
  @constraint(model, abs_cons[i=0:P-1], abs_var[i] == constraint_object(unique[i]).func - z)
  @objective(model, Max, f_1)
  optimize!(model)
  opt_f1 = objective_value(model)

  if termination_status(model) == MOI.OPTIMAL
    # obtain the nadir d_2
    d_2 = -sum(abs(value(abs_var[i])) for i=0:P-1)
    d_1 = 0.0
    vcount = Dict([i=>round(Int, value(unique[i])) for i=0:P-1])
    return vcount, [d_1, d_2]
  elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
    suboptimal_solution = value.(x)
    suboptimal_objective = objective_value(model)
  else
    error("The model was not solved correctly")
  end
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


function master_problem(vcount, nadir, P, Γ, K, L, pra_dict)
  # start time
  stime = time() 
  # Set of sensitized patients 
  if pra_dict != nothing
    sensitized = Set(i for (i,pra) in pra_dict if pra >= 0.8)
  else
    sensitized = Set()
  end
  n = length(Γ)
  d_1, d_2 = nadir

  # Main model i.e. Master Problem
  model = Model(Mosek.Optimizer) 
  set_optimizer_attribute(model, "MSK_IPAR_LOG", 0)
  # define variables
  δ = [@variable(model, lower_bound=0)]
  @variable(model, y[i=1:2])
  @variable(model, z_0)
  @variable(model, z[i=0:P-1])
  @variable(model, T)
  @variable(model, t[i=0:P-1])
  @variable(model, r)

  # define complicating constraints
  A = [Dict(i=>abs(vcount[i]) for i=0:P-1)]
  @constraint(model, c1, sum(δ) == 1)
  @constraint(model, c2, y[1] == sum(A[1][i] * δ[1] for i=0:P-1) - d_1) 
  @constraint(model, c3, y[2] == -T - d_2) 
  @constraint(model, c4, P * z_0 == sum(A[1][i] * δ[1] for i=0:P-1))
  @constraint(model, c5[i=0:P-1], z[i] == sum(A[1][i] * δ[1]) - z_0)
  @constraint(model, c6, [y[1], y[2], r] in RotatedSecondOrderCone())
  @constraint(model, c7[i=0:P-1], [t[i], z[i]] in SecondOrderCone())
  @constraint(model, c8, sum(t) == T)

  @objective(model, Min, -r)
  optimize!(model)

  status = termination_status(model)
  if status == MOI.OPTIMAL || status == MOI.SLOW_PROGRESS 
    opt_dist = value.(δ)
    L1 = objective_value(model)
  elseif termination_status(model) == MOI.TIME_LIMIT
    error("Time limit was reached")
  else
    println(termination_status(model))
    error("The model could not be solved correctly")
  end


  # Submodel i.e Subproblem
  submodel = Model(Mosek.Optimizer)
  set_optimizer_attribute(submodel, "MSK_IPAR_LOG", 0)

  @variable(submodel, x[l=0:n-2,(i,j)=Γ[l][2],K[i,j,l]], Bin)

  @constraint(submodel, unique[i=0:n-1], sum(x[l,(j,i),k] for l=0:i, (j,v) in Γ[l][2], k in K[j,v,l] if v == i) <= 1)
  con = Dict(i=>constraint_object(unique[i]).func for i=0:n-1)
  @constraint(submodel, flow[l=0:n-2,i=l+1:n-1,k=1:L-1], sum(x[l,(j,i),k] for (j,v) in Γ[l][2] if v == i && k in K[j,i,l]) == sum(x[l,(i,j),k+1] for (v,j) in Γ[l][2] if v == i && (k+1) in K[i,j,l]))

  while true
    π_0 = dual(c1)
    π_1 = dual(c2)
    π_2 = dual(c3)
    π_3 = dual(c4)
    β = dual.(c5)

    _vcount_s = Dict([i=>constraint_object(unique[i]).func for i=0:P-1])
    obj_expr = -π_0 + sum(_vcount_s[i] * (π_1 + π_3 + β[i]) for i=0:P-1)
    if !sub_problem(submodel, obj_expr, x)
      println("NSWP objective: ", objective_value(model))
      println("L1: ", value(T))
      println("Number of transplants: ", value(y[1]) + d_1)
      etime = time() - stime
      stats = Dict("L1"=>value(T), "NSWP"=>objective_value(model), "Number of transplants"=>value(y[1] + d_1), "time"=>etime)
      return stats
      break
    end
    
    # count vertices in solutions
    vcount_s = Dict([i=>round(Int, value(constraint_object(unique[i]).func)) for i=0:P-1])

    # add new coefficient and column
    push!(δ, @variable(model, lower_bound=0)) 
    push!(A, vcount_s)

    # modify the constraints by adding the new column
    for i=0:P-1
      set_normalized_coefficient(c5[i], δ[end], -vcount_s[i])
    end
    set_normalized_coefficient(c1, δ[end], 1)
    set_normalized_coefficient(c2, δ[end], -sum(vcount_s[i] for i=0:P-1))
    set_normalized_coefficient(c4, δ[end], -sum(vcount_s[i] for i=0:P-1))

    # reoptimize
    optimize!(model)
  end
end


function main()
  path = ARGS[1]
  L = parse(Int64, ARGS[2])
  
  G, num_P = read_input(path*".wmd")
  pra_dict = isfile(path*".dat") ? read_data(path*".dat") : nothing

  (Γ, d) = copy_shortest_paths(G)
  K = get_positions(G, d, L)

  filename = split(path, '/')[end]

  vcount, nadir = columngen_params(num_P, Γ, K, L, pra_dict)
  stats = master_problem(vcount, nadir, num_P, Γ, K, L, pra_dict)

  save("$filename.jld2", stats)

end


main()
