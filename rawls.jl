module Rawls
using JuMP, MosekTools
import ..Main: optimizer_status, generate_column_master!
export solve_subproblem!, update_constr!, init, fairness_objective!, reference_point!, partial_build_master_problem
export build_master_problem!, build_linear_combination!


function solve_subproblem!(model)
    submodel = model[:submodel]
    feasible = submodel[:feasible]
    α0 = dual(constraint_by_name(model, "α0"))
    α1 = dual(constraint_by_name(model, "α1"))
    β = Dict(i => dual(constraint_by_name(model, "β[$i]")) for i = feasible)

    # set the new objective
    pair_count = Dict(i => (constraint_object(constraint_by_name(submodel, "capacity[$i]"))).func for i = feasible)
    @objective(submodel, Min, -α0 + sum((α1 + β[i]) *  pair_count[i] for i = feasible))
    optimize!(submodel)
    optimizer_status(submodel)
    ζ = objective_value(submodel)
    if ζ > -1e-8
        return true 
    else
        return false
    end
end


function update_constr!(model, new_sol)
    P_ = keys(new_sol)
    α0 = constraint_by_name(model, "α0")
    α1 = constraint_by_name(model, "α1")
    β = Dict(i => constraint_by_name(model, "β[$i]") for i = P_)
    A = model[:A]
    
    # modify the constraints by adding the new column
    δ = @variable(model, lower_bound = 0.0, base_name = "δ[$(length(A) + 1)]")
    for i = P_
        set_normalized_coefficient(β[i], δ, -new_sol[i])
    end
    set_normalized_coefficient(α0, δ, 1)
    set_normalized_coefficient(α1, δ, -sum(new_sol[i] for i = P_))
end


function init(submodel)
    sol = Dict(i => (round ∘ value ∘ constraint_by_name)(submodel, "capacity[$i]") for i = submodel[:feasible])
    return sol
end


function fairness_objective!(model)
    return @expression(model, 1.0 * variable_by_name(model, "T")) 
end


function reference_point!(model, submodel, A)
    y1 = variable_by_name(model, "y[1]")

    # compute i1
    f1 = 1.0 * y1
    optimize!(submodel)
    optimizer_status(submodel)
    i1 = objective_value(submodel)

    # add temporary constraint
    temp = @constraint(model, f1 == objective_value(submodel))

    # change the objective to f2 to compute d2
    f2 = fairness_objective!(model)
    @objective(model, Max, f2)
    generate_column_master!(model, solve_subproblem!, update_constr!, A)
    d2 = objective_value(model)

    # delete temp constraint and optimize to get f2
    delete(model, temp)
    generate_column_master!(model, solve_subproblem!, update_constr!, A)
    i2 = objective_value(model)

    # add temp constraint and optimize to get d1
    temp = @constraint(model, f2 == i2)
    @objective(model, Max, 1.0 * y1)
    generate_column_master!(model, solve_subproblem!, update_constr!, A)
    d1 = objective_value(model)
    delete(model, temp)

    ideal = (i1, i2)
    nadir = (d1, d2)

    return ideal, nadir
end


function partial_build_master_problem(submodel, init_sol)
    P_ = submodel[:feasible]
    
    # Main model i.e. Master Problem
    model = Model(Mosek.Optimizer) 
    set_optimizer_attribute(model, "MSK_IPAR_LOG", 0)
    model[:submodel] = submodel

    # define variables
    δ = [@variable(model, lower_bound=0)]
    model[:δ] = δ
    @variable(model, y[i = 1:2])
    @variable(model, z[i = P_])
    @variable(model, T)
    @variable(model, t[i = P_])

    # define basic constraints
    A = [Dict(i => init_sol[i] for i = P_)]
    model[:A] = A
    @constraint(model, α0, sum(δ) == 1)
    @constraint(model, β[i = P_], z[i] == sum(A[1][i] * δ[1]))
    @constraint(model, γ[i = P_], z[i] >= T)

    # partially build complicating constraints
    @constraint(model, α1, y[1] == sum(A[j][i] * δ[j] for j = 1:length(A) for i = P_)) 

    return model
end


function build_master_problem!(model, ideal, nadir)
    i1, i2 = ideal
    d1, d2 = nadir
    y1 = variable_by_name(model, "y[1]")
    y2 = variable_by_name(model, "y[2]")
    α1 = constraint_by_name(model, "α1")
    T = variable_by_name(model, "T")

    set_normalized_rhs(α1, -d1)
    @variable(model, r)
    @constraint(model, α2, y2 == T - d2) 
    @constraint(model, u, [y1, y2, r] in RotatedSecondOrderCone())

    @objective(model, Max, r)
end


function build_linear_combination!(model, ideal, nadir)
    T = variable_by_name(model, "T")
    y1 = variable_by_name(model, "y[1]")
    A = model[:A]

    # weight the objectives
    i1, i2 = ideal
    d1, d2 = nadir
    if i1 != d1
        if i2 == d2
            w1 = 1.0
            w2 = 0.0
        else
            w1 = 1.0 / abs(i1 - d1)
            w2 = 1.0 / abs(i2 - d2)
        end
    else
        if i2 != d2
            w1 = 0.0
            w2 = 1.0
        else
            w1 = 0.0
            w2 = 0.0
        end
    end

    @objective(model, Max, w1 * y1 + w2 * T)
end


end