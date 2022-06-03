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


# function to retrieve the node at index i
function getnode(G, i::Int)
    if i > length(adjlistN(G))
        return adjlistP(G)[i - length(adjlistN(G))]
    end
    return adjlistN(G)[i]
end


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
function vindex(vertex::T)::Int where T <: Vertex
    return vertex.index
end


mutable struct Graph
    adjlistN::Vector{Node}
    adjlistP::Vector{Node}
end
Graph(numN, numP) = Graph(Vector{Node}(undef, numN), Vector{Node}(undef, numP))


function graphsize(G::Graph; N = true, P = true)
    total = 0
    @assert P || N
    if N
        total += length(adjlistN(G))
    end
    if P
        total += length(adjlistP(G))
    end
    return total
end


function adjlistN(G::Graph)
    G.adjlistN
end


function adjlistP(G::Graph)
    G.adjlistP
end


function addarc!(G, v, w)
    nodew = Node(w)
    
    setnextnode!(nodew, nextnode(getnode(G, vindex(v))))
    setnextnode!(getnode(G, vindex(v)), nodew)
end


# function to get all the arcs that satisfy predicate fn
function getarcs(G; fn = (x, y) -> true)
    arcs = Set{Tuple{Int, Int}}() 
    for i = 1:graphsize(G)
        node = getnode(G, i)

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


function read_input(dir, filename)
    path = joinpath(dir, filename) 
    cmd = pipeline(`python3 -m kidney_solver.utils.convert`, stdin=path * ".wmd")
    lines = readlines(cmd)
    sep = findall(x->x!="-1\t-1\t-1", lines)
    lines = lines[sep]
    lines = map(x->split(x, "\t"), lines)
    lines = [map(x->parse(Int,x), line) for line in lines]
    numP, numA = popfirst!(lines)
    numN, numE = splice!(lines, numA + 1)

    pra_dict = read_data(path * ".dat")

    G = Graph(numN, numP)
    
    for i = 1:numP
        j = i + numN
        adjlistP(G)[i] = Node(PDP(j, pra_dict[j]))
    end

    for i = 1:numN
        adjlistN(G)[i] = Node(NDD(i))
    end

    for i = 1:numA
        v, w = map(x -> PDP(x + 1 + numN, pra_dict[x + 1]), lines[i][1:2])
        addarc!(G, v, w)
    end

    for i = 1:numE
        _v, _w = map(x -> x + 1, lines[i][1:2])
        v, w = NDD(_v), PDP(_w + numN, pra_dict[_w])
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
    numV = graphsize(G)
    d = Dict{NTuple{3, Int}, Float64}((i, j, l) => Inf for l = 1:numV - 1 for i = l:numV for j = l:numV)
    Γ = Dict{Int, Graph}

    for l = graphsize(G):-1:1
        # compute the distances from each vertex to the other
        for i = 1:graphsize(G)
            node = getnode(G, i)
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
    K_ = Dict{NTuple{2, Int}, Set{Int64}}()
    
    for l = 1:graphsize(G) - 1
        for i = l:graphsize(G)
            for j = l:graphsize(G)
                K[i, j, l] = Set(k for k = 1:L if d[l, i, l] < k && d[j, l, l] <= L - k)
            end
        end
    end

    for (i, j) in getarcs(G)
        if i > graphsize(G; P = false)
            K_[i, j] = Set(2:L)
        else
            K_[i, j] = Set(1)
        end
    end

    return K, K_
end


# TODO normalize the squared distance between i1/f1 and i2/f2
# otherwise we get a skewed idea of the distance
function distance_to_ideal(sol, ideal, nadir)
    i1, i2 = ideal
    d1, d2 = nadir
    f1, f2 = sol

    norm = i1 - d1 > 0 ? (i1 - d1)^2 : 0.0
    norm += i2 - d2 > 0 ? (i2 - d2)^2 : 0.0
    norm = sqrt(norm)

    dist = i1 - f1 > 0 ? (i1 - f1)^2 : 0.0
    dist += i2 - f2 > 0 ? (i2 - f2)^2 : 0.0
    dist = sqrt(dist)

    return norm > 0 ? dist / norm : 0.0
end


function distance_to_nadir(sol, ideal, nadir)
    i1, i2 = ideal
    d1, d2 = nadir
    f1, f2 = sol

    norm = i1 - d1 > 0 ? (i1 - d1)^2 : 0.0
    norm += i2 - d2 > 0 ? (i2 - d2)^2 : 0.0
    norm = sqrt(norm)

    dist = f1 - d1 > 0 ? (f1 - d1)^2 : 0.0
    dist += f2 - d2 > 0 ? (f2 - d2)^2 : 0.0
    dist = sqrt(dist)

    return norm > 0 ? dist / norm : 0.0
end


function price_of_fairness(sol, ideal, nadir)
    i1, i2 = ideal
    d1, d2 = nadir
    f1, f2 = sol
    return i1 > 0 ? (i1 - f1) / i1 : 0.0
end


function price_of_utility(sol, ideal, nadir)
    i1, i2 = ideal
    d1, d2 = nadir
    f1, f2 = sol

    if i2 <= 0
        return d2 < 0 ? (i2 - f2) / abs(d2) : 0.0 
    else
        return i2 > 0 ? (i2 - f2) / i2 : 0.0
    end
end


function support_size(model, A)
    size = 0
    for i = 1:length(A)
        var = variable_by_name(model, "δ[$i]")
        if value(var) > 0
            size += 1
        end
    end

    return size
end


mutable struct Stats 
    time::Float64
    support::Int
end


# TODO
# Function to retrieve the constraints that are part of the same container.
# This is useful since after as model is saved to a file, we lose the container
# references inside the model.
function constraint(model, name, index_set)

end
