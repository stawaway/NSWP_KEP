using ArgParse, JLD2, MosekTools


s = ArgParseSettings()
@add_arg_table s begin
    "--exp_path", "-e"
        help = "path to the directory used to store experimental data"
        default = "experiments"
    "nswp_module"
        help = "NSWP module that is used"
        required = true
    "filename"
        help = "The kidney file name"
        required = true
end

# parse the command line arguments
parsed_args = parse_args(ARGS, s)
filename = parsed_args["filename"]
nswp_module = parsed_args["nswp_module"]
exp_path = parsed_args["exp_path"]

using JuMP, MosekTools
if nswp_module == "IF"
    modfile = "if.jl"
elseif nswp_module == "Nash"
    modfile = "nash.jl"
elseif nswp_module == "Rawls"
    modfile = "rawls.jl"
elseif nswp_module == "Aristotle"
    modfile = "aristotle.jl"
else
    error("Invalid NSWP module")
end
include("utils.jl")
include("optim.jl")
include("$modfile")
expr = Meta.parse("using .$(nswp_module)")
eval(expr)


function main()
    # Open file that has the experimental data
    fid = jldopen(joinpath(exp_path, filename * ".mof.json"), "a+")
    
    # retrieve feasible and sensitized sets
    feasible = fid["stats/feasible"]
    sensitized = fid["stats/sensitized"]

    # retrieve the initial set of solutions
    A = fid["models/model/solutions/$nswp_module"]
    
    # retrieve the model and submodel
    submodel = load_HPIEF(fid["models/submodel/path"], feasible, sensitized)
    set_optimizer(submodel, Mosek.Optimizer)
    set_optimizer_attribute(submodel, "MSK_IPAR_LOG", 0)
    model = load_model(fid["models/model/path/$nswp_module"], submodel, A)
    set_optimizer(model, Mosek.Optimizer)
    set_optimizer_attribute(model, "MSK_IPAR_LOG", 0)

    # solve the model
    starttime = time()
    solve!(model, solve_subproblem!, update_constr!)

    ideal = fid["stats/ideal/$nswp_module"]
    nadir = fid["stats/nadir/$nswp_module"]
    i1, i2 = ideal
    d1, d2 = nadir

    # save the experiment's stats
    # TODO also include the POF and POU
    sol = solution(model, reference = nadir)
    if ! haskey(fid, "stats/support_size/model")
        fid["stats/support_size/model"] = length(model[:A])
    end
    fid["stats/time/$nswp_module/model"] = time() - starttime
    fid["stats/ideal_distance/model/$nswp_module"] = distance_to_ideal(sol, ideal, nadir)
    fid["stats/nadir_distance/model/$nswp_module"] = distance_to_nadir(sol, ideal, nadir)
    fid["stats/pof/model/$nswp_module"] = price_of_fairness(sol, ideal, nadir)
    fid["stats/pou/model/$nswp_module"] = price_of_utility(sol, ideal, nadir)
#    println("Size of the support is ", length(model[:A]))
#    println("Distance to ideal is ", distance_to_ideal(sol, ideal, nadir))
#    println("Distance to nadir is ", distance_to_nadir(sol, ideal, nadir))
#    println("Time spent by the method ", time() - starttime)
#    println("POF is ", price_of_fairness(sol, ideal, nadir))
#    println("POU is ", price_of_utility(sol, ideal, nadir))
    
    close(fid)
end


main()

