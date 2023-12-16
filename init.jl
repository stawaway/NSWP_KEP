using ArgParse, JLD2, MosekTools


s = ArgParseSettings()
@add_arg_table s begin
    "--dir", "-d"
        help = "directory of the kidney data"
        default = "kidney"
    "--model_path", "-m"
        help = "path to the directory used to store models"
        default = "models"
    "--exp_path", "-e"
        help = "path to the directory used to store models"
        default = "experiments"
    "--skip", "-s"
        help = "skip rebuilding the graph and HPIEF submodel"
        action = :store_true
    "--cycle_size"
        help = "The cap on cycle size"
        arg_type = Int
        default = 3
    "--path_size"
        help = "The cap on path starting with an NDD"
        arg_type = Int
        default = 3
    "nswp_module"
        help = "NSWP module that is used"
        required = true
    "filename"
        help = "The kidney file name"
        required = true
end

# parse the command line arguments
parsed_args = parse_args(ARGS, s)
dir = parsed_args["dir"]
filename = parsed_args["filename"]
nswp_module = parsed_args["nswp_module"]
model_path = parsed_args["model_path"]
cycle_size = parsed_args["cycle_size"]
path_size = parsed_args["path_size"]
exp_path = parsed_args["exp_path"]
skip = parsed_args["skip"]

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


macro setup()
    ex1 = quote
        L = cycle_size
        G = read_input(dir, filename)
        d = copy_shortest_paths(G)
        K, K_ = get_positions(G, d, L)
        submodel = build_HPIEF(G, d, K, K_, L)
        feasible = feasible_subgraph!(submodel, G, d, K, K_, L)
        sensitized = submodel[:sensitized]

        # save stats
        fid["stats/G"] = G
        fid["stats/K"] = K
        fid["stats/K_"] = K_
        fid["stats/d"] = d
        fid["stats/L"] = L

        # save submodel to file
        save_HPIEF(submodel, fid, joinpath(model_path, "submodels", filename * ".mof.json"))
    end

    ex2 = quote
        # retrieve groups
        stats_group = fid["stats"]

        # obtain the stats and graph information
        feasible = stats_group["feasible"]
        sensitized = stats_group["sensitized"]
        G = stats_group["G"]
        K = stats_group["K"]
        K_ = stats_group["K_"]
        d = stats_group["d"]
        L = stats_group["L"]

        # obtain the submodel
        submodel = load_HPIEF(fid["models/submodel/path"], fid["stats/feasible"], fid["stats/sensitized"])
        set_optimizer(submodel, Mosek.Optimizer)
        set_optimizer_attribute(submodel, "MSK_IPAR_LOG", 0)
    end

    if ! skip
        return esc(ex1)
    else
        return esc(ex2)
    end
end


function module_specific!(fid, module_name, submodel)
    stats = Stats(0.0, 1, Dict{Int, Float64}())
    single_time = time()
    optimize!(submodel)
    init_sol = init(submodel)
    
    # partially build the main model and find ideal / reference point
    model = partial_build_master_problem(submodel, init_sol)
    ideal, nadir = reference_point!(model, submodel, model[:A], stats = stats)
 
    # save the time and support size
    fid["stats/time/single/$nswp_module"] = time() - single_time
    fid["stats/reference_time/$nswp_module"] = stats.time
    fid["stats/support_size/single/$nswp_module"] = stats.support
    fid["stats/probs/single/$nswp_module"] = stats.solution

    # build linear model and save it
    build_linear_combination!(model, ideal, nadir)
    save_linear(model, fid, module_name, joinpath(model_path, nswp_module, filename * "_linear.mof.json"))
  
    # build master problem and save it
    build_master_problem!(model, ideal, nadir)
    save_model(model, fid, module_name, joinpath(model_path, nswp_module, filename * ".mof.json"))

    # save stats
    fid["models/submodel/solutions/$module_name"] = init_sol
    fid["stats/ideal/$nswp_module"] = ideal
    fid["stats/nadir/$nswp_module"] = nadir
end


function main()
    # create paths to directories and subdirectories if they don't exist
    mkpath(model_path)
    mkpath(joinpath(model_path, "submodels"))
    mkpath(joinpath(model_path, nswp_module))
    mkpath(exp_path)

    # open JLD2 file and create it if it doesn't exist
    fid = jldopen(joinpath(exp_path, filename * ".mof.json"), "a+")

    # check if a file exists that stores the Graph and submodel information
    @setup()

    module_specific!(fid, nswp_module, submodel)
    close(fid)
end

main()

