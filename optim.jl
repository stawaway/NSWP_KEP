# This file defines functions that are used to build the PIEF model
# and do column generation on an arbitrary model


function build_HPIEF(G, d, K, K_, L; threshold = 0.8)
    sensitized = Set{Int}()
    for node in adjlistP(G)
        @assert typeof(vertex(node)) == PDP
        if pra(vertex(node)) >= threshold
            push!(sensitized, vindex(vertex(node)))
        end
    end

    # define the constraints of the PIEF submodel
    submodel = Model(Mosek.Optimizer)
    submodel[:sensitized] = sensitized
    set_optimizer_attribute(submodel, "MSK_IPAR_LOG", 0)

    # create variables x and y for the submodel
    @variable(submodel, 
        x[l = 1 + graphsize(G; P = false):graphsize(G), 
        (i, j) = getarcs(G; fn = (x, y) -> min(x, y) >= l),
        K[i, j, l]], 
        Bin)
    @variable(submodel, 
        y[(i, j) = getarcs(G), k = K_[i, j]], 
        Bin)

    unique_expr = [AffExpr(0.0) for i = 1:graphsize(G; P = false)]  
    for i = 1:graphsize(G; P = false)
        for (_, j) in getarcs(G; fn = (x, y) -> x == i)
            add_to_expression!(unique_expr[i], y[(i, j), 1]) 
        end
    end
    @constraint(submodel, unique[i = 1:graphsize(G; P = false)], unique_expr[i] <= 1)

    # build capacity constraint ensuring each pair is in one chain or cycle
    capacity_expr = Dict{Int, AffExpr}(i => AffExpr(0.0) for i = 1 + graphsize(G; P = false):graphsize(G))
    for i = 1 + graphsize(G; P = false):graphsize(G)
        for l = 1 + graphsize(G; P = false):graphsize(G)
            for (j, _) in getarcs(G; fn = (x, y) -> min(x, y) >= l && y == i)
                for k in K[j, i, l]
                    add_to_expression!(capacity_expr[i], x[l, (j, i), k])
                end
            end
        end

        for (j, _) in getarcs(G; fn = (x, y) -> y == i)
            for k in K_[j, i]
                add_to_expression!(capacity_expr[i], y[(j, i), k])
            end
        end
    end
    @constraint(submodel, capacity[i = 1 + graphsize(G; P = false):graphsize(G)], capacity_expr[i] <= 1) 

    # build chain constraint to ensure chains involve positions in sequence
    chain_expr1 = Dict{NTuple{2, Int}, AffExpr}()
    chain_expr2 = Dict{NTuple{2, Int}, AffExpr}()
    for i = 1 + graphsize(G; P = false):graphsize(G)
        for k = 1:L - 1
            chain_expr1[i, k] = AffExpr(0.0)
            chain_expr2[i, k] = AffExpr(0.0)
        end
    end
    for i = 1 + graphsize(G; P = false):graphsize(G)
        for k = 1:L - 1
            for (v, w) in getarcs(G)
                if w == i && k in K_[v, i]
                    add_to_expression!(chain_expr1[i, k], y[(v, i), k])
                end

                if v == i
                    add_to_expression!(chain_expr2[i, k], y[(i, w), k + 1])
                end
            end
        end
    end
    @constraint(submodel,
        chain[i = 1 + graphsize(G; P = false):graphsize(G), k = 1:L - 1],
        chain_expr1[i, k] >= chain_expr2[i, k])

    flow_expr1 = Dict((l, i, k) => AffExpr(0.0) for l = 1:graphsize(G) for i = l + 1:graphsize(G) for k = 1:L-1)
    flow_expr2 = Dict((l, i, k) => AffExpr(0.0) for l = 1:graphsize(G) for i = l + 1:graphsize(G) for k = 1:L-1)
    for l = 1:graphsize(G)
        for (v, w) in getarcs(G; fn = (x, y) -> min(x, y) >= l)
            if w >= l + 1 
                for k in K[v, w, l]
                    add_to_expression!(flow_expr1[l, w, k], x[l, (v, w), k])
                end
            end

            if v >= l + 1
                for k in K[v, w, l]
                    add_to_expression!(flow_expr2[l, v, k - 1], x[l, (v, w), k])
                end
            end
        end
    end
    @constraint(submodel, 
        flow[l = 1 + graphsize(G; P = false):graphsize(G) - 1, i = l + 1:graphsize(G), k = 1:L - 1], 
        flow_expr1[l, i, k] == flow_expr2[l, i, k])

    f_1 = sum(constraint_object(capacity[i]).func for i = 1 + graphsize(G; P = false):graphsize(G))
    @objective(submodel, Max, f_1)
    return submodel
end


function optimizer_status(model)
    status = termination_status(model)
    if status == MOI.OPTIMAL
        return status
    elseif status == MOI.INFEASIBLE
        error("The model is infeasible")
    elseif status == MOI.TIME_LIMIT
        error("Time limit was reached")
    elseif status == MOI.SLOW_PROGRESS
        return status
    else
        println(status)
        error("The model could not be solved correctly")
    end
end


function feasible_subgraph!(submodel, G, d, K, K_, L)
    x, y = submodel[:x], submodel[:y]
    unique = submodel[:unique]
    capacity = submodel[:capacity]
    chain = submodel[:chain]
    flow = submodel[:flow]

    feasible = Set{Int}()
    size_feasible = 0

    optimize!(submodel)
    optimizer_status(submodel)
    size_feasible = length(feasible)

    while true
        # add constraints to find other vertices that can be matched
        feasible_expr = AffExpr(0.0)
        for i = 1 + graphsize(G; P = false):graphsize(G)
            if !(i in feasible)
                add_to_expression!(feasible_expr, constraint_object(capacity[i]).func)
            end
        end
        temp = @constraint(submodel, feasible_expr >= 1.0)
        
        _size_feasible = size_feasible
        optimize!(submodel)
        if termination_status(submodel) == MOI.INFEASIBLE
            delete(submodel, temp)
            break
        end
        for i = 1 + graphsize(G; P = false):graphsize(G)
            if value(constraint_object(capacity[i]).func) == 1.0
                push!(feasible, i)
            end
        end

        if length(feasible) > _size_feasible
            size_feasible = length(feasible)
        else
            delete(submodel, temp)
            break
        end
        
        delete(submodel, temp)
    end

    # delte x constraints that are irrelevant
    for l = 1 + graphsize(G; P = false):graphsize(G)
        for (v, w) in getarcs(G; fn = (x, y) -> min(x, y) >= l)
            if !(v in feasible) || !(w in feasible)
                for k in K[v, w, l]
                    delete(submodel, x[l, (v, w), k])
                end
            end
        end
    end

    # do the same for the y variables
    for (i, j) in getarcs(G)
        for k in K_[i, j]
            if !(i in feasible) && i > graphsize(G; P = false) || !(j in feasible)
                delete(submodel, y[(i, j), k])
            end
        end
    end

    # delete constraints that will not be used
    for i = 1:graphsize(G)
        if !(i in feasible) && i <= graphsize(G; P = false)
            delete(submodel, unique[i])
        end

        if !(i in feasible) && i > graphsize(G; P = false)
            delete(submodel, capacity[i])

            for k = 1:L -1
                delete(submodel, chain[i, k])
            end
        end

        l = i
        if l <= graphsize(G; P = false)
            continue
        end
        for i = l + 1:graphsize(G)
            if !(l in feasible) || !(i in feasible)
                for k = 1:L - 1
                    delete(submodel, flow[l, i, k])
                end
            end
        end
    end

    submodel[:feasible] = filter(x -> x > graphsize(G; P = false), feasible)
    return feasible
end


function generate_column_master!(model, submodel_fn, update_constr_fn, A)
    submodel = model[:submodel] # TODO add function to retrieve solution instead
    optimize!(model)
    optimizer_status(model)

    while true
        # solve the submodel
        # add the solution to the Matrix A if it exists
        if ! submodel_fn(model)
            sol = Dict(i => round(value(submodel[:capacity][i])) for i = submodel[:feasible])
        else
            break
        end

        # update the constraints using the new set of solutions
        if sol != A[end]
            push!(A, sol)
            update_constr_fn(model, sol)

            optimize!(model)
            optimizer_status(model)
        end
    end

    return A
end


function solve!(model, submodel_fn, update_constr_fn)
    generate_column_master!(model, submodel_fn, update_constr_fn, model[:A]) 
end


#end
