using ArgParse, JLD2


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


L = 3
G = read_input(dir, filename)
d = copy_shortest_paths(G)
K, K_ = get_positions(G, d, 3)
submodel = build_HPIEF(G, d, K, K_, L)
feasible = feasible_subgraph!(submodel, G, d, K, K_, L)


optimize!(submodel)
init_sol = init(submodel)
model = partial_build_master_problem(submodel, init_sol)
ideal, nadir = reference_point!(model, submodel, model[:A])
build_master_problem!(model, ideal, nadir)

# create paths to directories and subdirectories
mkpath(model_path)
mkpath(joinpath(model_path, "submodels"))
mkpath(joinpath(model_path, "models"))
mkpath(exp_path)

# write the models to file
write_to_file(submodel, joinpath(model_path, "submodels", filename * ".cbf"))
write_to_file(model, joinpath(model_path, "models", filename * ".cbf"))

# open JLD2 file and create it if it doesn't exist
fid = jldopen(joinpath(exp_path, filename * ".hdf"), "w")

# create groups to store relevant information
model_group = JLD2.Group(fid, "models")
stats_group = JLD2.Group(fid, "stats")

# store models
model_group["submodel"] = joinpath(model_path, "submodels", filename * ".cbf")
model_group["model"] = joinpath(model_path, "models", filename * ".cbf")

# store some stats
stats_group["feasible"] = feasible
stats_group["G"] = G
stats_group["K"] = K
stats_group["K_"] = K_
stats_group["d"] = d
stats_group["L"] = 3

close(fid)

