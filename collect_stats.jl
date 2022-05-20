using JLD2, HDF5, ArgParse
include("utils.jl")


s = ArgParseSettings()
@add_arg_table s begin
    "--dir", "-d"
        help = "directory containing the statistics file for each graph"
        default = "experiments"
    "output"
        help = "output file name containing the experimental data"
        required = true
end

# parse the command line arguments
parsed_args = parse_args(ARGS, s)
dir = parsed_args["dir"]
output = parsed_args["output"]


function file_data(filename)
    fid = jldopen(filename, "r")
    linear_dict = Dict{String, Float64}()
    model_dict = Dict{String, Float64}()

    if haskey(fid, "stats/G")
        G = fid["stats/G"]
        linear_dict["|P|"] = length(adjlistP(G))
        model_dict["|P|"] = length(adjlistP(G))

        linear_dict["|N|"] = length(adjlistN(G))
        model_dict["|N|"] = length(adjlistN(G))
    else
        println(fid) # TODO Remove this else after you are done
    end

    schemes = ["IF", "Rawls", "Nash", "Aristotle"]
    for scheme in schemes
        # get the time
        if haskey(fid, "stats/time/linear/$scheme")
            time_linear = fid["stats/time/linear/$scheme"]
            linear_dict["time"] = time_linear
        end

        if haskey(fid, "stats/time/model/$scheme")
            time_model = fid["stats/time/model/$scheme"]
            model_dict["time"] = time_model
        end

        # get the distance to ideal
        if haskey(fid, "stats/ideal_distance/linear/$scheme")
            ideal_dist_linear = fid["stats/ideal_distance/linear/$scheme"]
            linear_dict["ideal_distance"] = ideal_dist_linear 
        end

        if haskey(fid, "stats/ideal_distance/model/$scheme")
            ideal_dist_model = fid["stats/ideal_distance/model/$scheme"]
            model_dict["ideal_distance"] = ideal_dist_model
        end

        # get the distance to nadir
        if haskey(fid, "stats/nadir_distance/linear/$scheme")
            nadir_dist_linear = fid["stats/nadir_distance/linear/$scheme"]
            linear_dict["nadir_distance"] = nadir_dist_linear
        end

        if haskey(fid, "stats/nadir_distance/model/$scheme")
            nadir_dist_model = fid["stats/nadir_distance/model/$scheme"]
            model_dict["nadir_distance"] = nadir_dist_model
        end

        # get the POF
        if haskey(fid, "stats/pof/linear/$scheme")
            pof_linear = fid["stats/pof/linear/$scheme"]
            linear_dict["POF"] = pof_linear
        end

        if haskey(fid, "stats/pof/model/$scheme")
            pof_model = fid["stats/pof/model/$scheme"]
            model_dict["POF"] = pof_model
        end

        # get the POU
        if haskey(fid, "stats/pou/linear/$scheme")
            pou_linear = fid["stats/pof/linear/$scheme"]
            linear_dict["POU"] = pou_linear
        end

        if haskey(fid, "stats/pou/model/$scheme")
            pou_model = fid["stats/pof/model/$scheme"]
            model_dict["POU"] = pou_model
        end
    end
    
    close(fid)
    return linear_dict, model_dict
end


function collect_stats(dir, output)
    # open file to store the data
    fid = h5open(output, "w") # TODO take care of this line
    linear = Dict{String, Vector{Float64}}("time" => [], "ideal_distance" => [], "nadir_distance" => [],
        "POF" => [], "POU" => [], "|P|" => [], "|N|" => []) 
    model = Dict{String, Vector{Float64}}("time" => [], "ideal_distance" => [], "nadir_distance" => [], 
        "POF" => [], "POU" => [], "|P|" => [], "|N|" => []) 

    # loop over files and collect data
    for filename in Base.Filesystem.readdir(dir)
        if split(filename, '.', limit = 2)[end] != "mof.json"
            continue
        end
        # store the data in the file
        linear_dict, model_dict = file_data(joinpath(dir, filename))

        for (key, val) in linear_dict
            push!(linear[key], val)
        end

        for (key, val) in model_dict
            push!(model[key], val)
        end
    end

    # save the experimental data for the model with linear combination of objectives
    for (key, val) in linear
        fid["linear/$key"] = val
    end

    # save the experimental data for the main model
    for (key, val) in model
        fid["model/$key"] = val
    end

    close(fid)
end


collect_stats(dir, output)
