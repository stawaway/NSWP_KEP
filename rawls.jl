module Rawls
using JuMP, MosekTools
import ..Main: optimizer_status, generate_column_master!
export solve_subproblem!, update_constr!, init, fairness_objective!, reference_point!, partial_build_master_problem
export build_master_problem!, build_linear_combination


function solve_subproblem!(model)
    submodel = model[:submodel]
    feasible = submodel[:feasible]
    α0 = dual(model[:α0])
    α1 = dual(model[:α1])
    β = dual.(model[:β])

    # set the new objective
    pair_count = Dict(i => constraint_object(submodel[:capacity][i]).func for i = feasible)
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
    α0, α1, β, δ = model[:α0], model[:α1], model[:β], model[:δ]
    # modify the constraints by adding the new column
    P_ = keys(new_sol)
    push!(δ, @variable(model, lower_bound = 0.0))
    for i = P_
        set_normalized_coefficient(β[i], δ[end], -new_sol[i])
    end
    set_normalized_coefficient(α0, δ[end], 1)
    set_normalized_coefficient(α1, δ[end], -sum(new_sol[i] for i = P_))
end


function init(submodel)
    sol = Dict(i => round(value(submodel[:capacity][i])) for i = submodel[:feasible])
    return sol
end


function fairness_objective!(model)
    return @expression(model, 1.0 * model[:T]) 
end


function reference_point!(model, submodel, A)
    y = model[:y]
    δ = model[:δ]

    # compute i1
    f1 = 1.0 * y[1]
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
    @objective(model, Max, 1.0 * y[1])
    generate_column_master!(model, solve_subproblem!, update_constr!, A)
    d1 = objective_value(model)
    delete(model, temp)

    ideal = (i1, i2)
    nadir = (d1, d2)

    return ideal, nadir
end


function partial_build_master_problem(submodel, init_sol)
    x = submodel[:x]
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
    y, α1, T = model[:y], model[:α1], model[:T]

    set_normalized_rhs(model[:α1], -d1)
    @variable(model, r)
    @constraint(model, α2, y[2] == T - d2) 
    @constraint(model, u, [y[1], y[2], r] in RotatedSecondOrderCone())

    @objective(model, Max, r)
end


function build_linear_combination(model, ideal, nadir)

end


end
