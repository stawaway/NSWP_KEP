using JuMP, MosekTools, Combinatorics, FileIO


abstract type Vertex end


struct PDP <: Vertex 
    index::Int
    pra::Float64
end


struct NDD <: Vertex 
    index::Int
end


mutable struct Node
    vertex::Vertex
    next::Union{Node, Nothing}
end
Node(v, w) = new(v, w)
Node(v) = Node(v, nothing)


# function to retrieve the next node
function nextnode(node::Node)
    return node.next
end


function setnextnode!(node::Node, next)
    node.next = next
end


# function to retrieve the vertex of the node
function vertex(node::Node)
    return node.vertex
end


# function to retrieve the pra of a PDD
function pra(vertex::PDP)
    vertex.pra
end


# function to retrieve the index of a vertex
function vindex(vertex::Vertex)::Int
    return vertex.index
end


mutable struct Graph
    numV::Int
    adjlist::Vector{Node}
end
Graph(numV) = Graph(numV, Vector{Node}(undef, numV))


function graphsize(G::Graph)
    G.numV
end


function adjlist(G::Graph)
    G.adjlist
end


function addarc!(G, v, w)
    nodew = Node(w)
     
    setnextnode!(nodew, nextnode(adjlist(G)[vindex(v)]))
    setnextnode!(adjlist(G)[vindex(v)], nodew)
end


# function to get all the arcs that satisfy predicate fn
function getarcs(G; fn = (x, y) -> true)
    arcs = Set{Tuple{Int, Int}}() 
    for i = 1:graphsize(G)
        node = adjlist(G)[i]

        while !isnothing(nextnode(node))
            next = nextnode(node)
            j = vindex(vertex(next))

            if fn(i, j)
                push!(arcs, (i, j))
            end

            node = nextnode(node)
        end
    end

    return arcs
end


function read_input(filename)
    cmd = pipeline(`python3 -m kidney_solver.utils.convert`, stdin=filename * ".wmd")
    lines = readlines(cmd)
    sep = findall(x->x!="-1\t-1\t-1", lines)
    lines = lines[sep]
    lines = map(x->split(x, "\t"), lines)
    lines = [map(x->parse(Int,x), line) for line in lines]
    numP, numA = popfirst!(lines)
    numN, numE = splice!(lines, numA + 1)

    pra_dict = read_data(filename * ".dat")

    G = Graph(numP + numN)
    
    for i = 1:numP
        adjlist(G)[i] = Node(PDP(i, pra_dict[i]))
    end

    for i = 1:numN
        adjlist(G)[i] = Node(NDD(i))
    end

    for i = 1:numA
        v, w = map(x -> PDP(x + 1, pra_dict[x + 1]), lines[i][1:2])
        addarc!(G, v, w)
    end

    for i = 1:numE
        _v, _w = map(x -> x + 1, lines[i][1:2])
        v, w = NDD(_v), PDP(_w, pra_dict[_w])
        addarc!(G, v, w)
    end

    return G
end


function read_data(filename)
    lines = readlines(filename)
    lines = map(x->split(x, ","), lines)
    col_idx = findfirst(x->x=="%Pra", popfirst!(lines))
    pra_keys = map(x->parse(Int, x), [line[1] for line in lines])
    if isnothing(col_idx)
        pra_vals = Float64[]
    else
        pra_vals = map(x->parse(Float64, x), [line[col_idx] for line in lines])
    end
    pra_dict = Dict(zip(pra_keys,pra_vals))
end


function copy_shortest_paths(G)
    d = Dict{NTuple{3, Int}, Float64}((i, j, l) => Inf for l = 1:G.numV - 1 for i = l:G.numV for j = l:G.numV)
    Γ = Dict{Int, Graph}

    for l = graphsize(G):-1:1
        # compute the distances from each vertex to the other
        for i = 1:graphsize(G)
            node = adjlist(G)[i]
            d[i, i, l] = 0

            while !isnothing(nextnode(node))
                next = nextnode(node)
                j = vindex(vertex(next))
                d[i, j, l] = 1
                node = next
            end
        end

        for k = l:graphsize(G)
            for i = l:graphsize(G)
                for j = l:graphsize(G)
                    if d[i, j, l] > d[i, k, l] + d[k, j ,l]
                        d[i, j, l] = d[i, k, l] + d[k, j, l]
                    end
                end
            end
        end
    end

    return d
end


function get_positions(G, d, L)
    K = Dict{NTuple{3, Int}, Set{Int64}}()
    
    for l = 1:graphsize(G) - 1
        for i = l:graphsize(G)
            for j = l:graphsize(G)
                K[i, j, l] = Set(k for k = 1:L if d[l, i, l] < k && d[j, l, l] <= L - k)
            end
        end
    end

    return K
end


function build_PIEF(G, d, K, L; threshold = 0.8)
    sensitized = Set{Int}()
    for node in adjlist(G)
        if pra(vertex(node)) >= threshold
            push!(sensitized, vindex(vertex(node)))
        end
    end

    # define the constraints of the PIEF submodel
    submodel = Model(Mosek.Optimizer)
    set_optimizer_attribute(submodel, "MSK_IPAR_LOG", 0)

    @variable(submodel, x[l = 1:graphsize(G), (i, j) = getarcs(G; fn = (x, y) -> min(x, y) >= l), K[i, j, l]], Bin)

    unique_expr = [AffExpr(0.0) for i = 1:graphsize(G)]  
    for i = 1:graphsize(G)
        for l = 1:i
            for (v, w) in getarcs(G; fn = (x, y) -> min(x, y) >= l)
                if w == i
                    for k in K[v, i, l]
                        add_to_expression!(unique_expr[i], x[l, (v, i), k])
                    end
                end
            end
        end
    end
    @constraint(submodel, unique[i = 1:graphsize(G)], unique_expr[i] <= 1)

    flow_expr1 = AffExpr(0.0)
    flow_expr2 = AffExpr(0.0)
    for l = 1:graphsize(G)
        for (v, w) in getarcs(G; fn = (x, y) -> min(x, y) >= l)
            if w >= l + 1 
                for k in K[v, w, l]
                    add_to_expression!(flow_expr1, x[l, (v, w), k])
                end
            end

            if v >= l + 1
                for k in K[v, w, l]
                    add_to_expression!(flow_expr2, x[l, (v, w), k])
                end
            end
        end
    end
    @constraint(submodel, flow[l = 0:graphsize(G) - 1, i = l + 1:graphsize(G), k = 1:L - 1], flow_expr1 == flow_expr2)

    f_1 = sum(constraint_object(unique[i]).func for i = 1:graphsize(G))
    @objective(submodel, Max, f_1)
    return submodel
    
    """
    # obtain the first solution and the nadir point
    vcount = columngen_params(submodel, unique, num_P)
    P = Set{Int64}(keys(vcount))

    # delete the constraints that will not be used from the submodel
    for i=0:num_P-1
      if haskey(vcount, i)
        continue
      else
        temp = unique[i]
        delete(submodel, temp)
        unregister(submodel, :temp)
      end
    end

    for l=0:min(num_P-1,n-2)
      if haskey(vcount, l)
        continue
      else
        for (i,j) in Γ[l][2]
          if !isempty(K[i,j,l])
            delete(submodel, x[l,(i,j),K[i,j,l]])
          end
        end
      end
      
      for i=l+1:n-1
        for k=1:L-1
          delete(submodel, flow[l,i,k])
        end
      end
    end
    num_P = length(keys(vcount))
    """
end


function feasible_subgraph!(submodel)
    x = submodel[:x]
    unique = submodel[:unique]

    feasible = Set{Int}()
    size_feasible = 0

    optimize!(submodel)
    status = termination_status(submodel)
    if status == MOI.INFEASIBLE
        error("The model is infeasible")
    end
    size_feasible = length(feasible)

    while true
        # add constraints to find other vertices that can be matched
        feasible_expr = AffExpr(0.0)
        for i = 1:graphsize(G)
            if !(i in feasible)
                add_to_expression!(feasible_expr, constraint_object(unique[i]).func)
            end
        end
        temp = @constraint(submodel, feasible_expr >= 1.0)
        
        _size_feasible = size_feasible
        optimize!(submodel)
        
        status = termination_status(submodel)
        if status == MOI.INFEASIBLE
            break
        end
        for i = 1:graphsize(G)
            if value(constraint_object(unique[i]).func) == 1.0
                push!(feasible, i)
            end
        end

        if length(feasible) > _size_feasible
            size_feasible = length(feasible)
        else
            break
        end
        
        delete(submodel, temp)
    end

    for l = 1:graphsize(G)
        for (v, w) in getarcs(G; fn = (x, y) -> min(x, y) >= l)
            if !(v in feasible) || !(w in feasible)
                for k in K[v, w, l]
                    delete(submodel, x[l, (v, w), k])
                end
            end
        end
    end

    for i = 1:graphsize(G)
        if !(i in feasible)
            delete(submodel, unique[i])
        end
    end

    return feasible
end


#G = read_input("kidney/MD-00001-00000080")
#d = copy_shortest_paths(G)
#K = get_positions(G, d, 3)
#submodel = build_PIEF(G, d, K, 3)
#feasible = feasible_subgraph!(submodel)
