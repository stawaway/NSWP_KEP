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


function restrict(val, lb, ub)
    @assert lb <= ub
    if val == -1.0
        return val
    end
    retval = lb <= val ? val : lb
    retval = retval <= ub ? retval : ub

    return retval
end


function price_if(probs, ideal, nadir)
    avg = sum(val for (_, val) in probs) / length(probs)
    f2 = sum(abs(val - avg) for (_, val) in probs)
    f1 = sum(val for (_, val) in probs)
    pof = price_of_fairness((f1, f2), ideal, nadir)

    return pof, f2
end


function price_rawls(probs, ideal, nadir)
    f2 = minimum(values(probs))
    f1 = sum(val for (_, val) in probs)
    pof = price_of_fairness((f1, f2), ideal, nadir)

    return pof, f2
end


function price_aristotle(probs, ideal, nadir, sensitized)
    f2 = sum([val for (key, val) in probs if key in sensitized], init = 0.0)
    f1 = sum(val for (_, val) in probs)
    pof = price_of_fairness((f1, f2), ideal, nadir)

    return pof, f2
end


function price_nash(probs, ideal, nadir)
    f2 = sum(val > 0.0 ? log(val) : -Inf for (_, val) in probs)
    f1 = sum(val for (_, val) in probs)
    pof = price_of_fairness((f1, f2), ideal, nadir)

    return pof, f2
end


function price_utility(scheme, ideal, nadir, numP)
    total = round(ideal[1])
    if scheme == "IF"
        avg = total / numP
        return sum(1.0 - avg for i = 1:total; init = 0.0) + sum(abs(0.0 - avg) for i = total + 1:numP; init = 0.0)
    end
    if scheme == "Rawls"
        return total < numP ? 0.0 : 1.0
    end
    if scheme == "Aristotle"
        return nadir[2]
    end
    if scheme == "Nash"
        return total < numP ? -Inf : 0.0
    end

    return 0.0
end


function _price_objectives!(fid, dicts, group, scheme, price_fn)
    single_dict, linear_dict, model_dict = dicts

    if haskey(fid, "stats/ideal/$group") && haskey(fid, "stats/probs/single/$scheme")
        ideal = fid["stats/ideal/$group"]
        nadir = fid["stats/nadir/$group"]
        probs = fid["stats/probs/single/$scheme"]
        pof_single, pou_single = price_fn(probs, ideal, nadir)
        pof_single = restrict(pof_single, 0.0, 1.0)

        single_dict["POF-$group/$scheme"] = pof_single
        single_dict["POU-$group/$scheme"] = pou_single
    else
        single_dict["POF-$group/$scheme"] = -1.0
        single_dict["POU-$group/$scheme"] = -1.0
    end

    if haskey(fid, "stats/ideal/$group") && haskey(fid, "stats/probs/linear/$scheme")
        ideal = fid["stats/ideal/$group"]
        nadir = fid["stats/nadir/$group"]
        probs = fid["stats/probs/linear/$scheme"]
        pof_linear, pou_linear = price_fn(probs, ideal, nadir)
        pof_linear = restrict(pof_linear, 0.0, 1.0)
        linear_dict["POF-$group/$scheme"] = pof_linear
        linear_dict["POU-$group/$scheme"] = pou_linear
    else
        linear_dict["POF-$group/$scheme"] = -1.0
        linear_dict["POU-$group/$scheme"] = -1.0
    end

    if haskey(fid, "stats/ideal/$group") && haskey(fid, "stats/probs/model/$scheme")
        ideal = fid["stats/ideal/$group"]
        nadir = fid["stats/nadir/$group"]
        probs = fid["stats/probs/model/$scheme"]
        pof_model, pou_model = price_fn(probs, ideal, nadir)
        pof_model = restrict(pof_model, 0.0, 1.0)
        model_dict["POF-$group/$scheme"] = pof_model
        model_dict["POU-$group/$scheme"] = pou_model
    else
        model_dict["POF-$group/$scheme"] = -1.0
        model_dict["POU-$group/$scheme"] = -1.0
    end
end


function price_objectives!(fid, dicts, scheme)
    if haskey(fid, "stats/sensitized")
        sensitized = fid["stats/sensitized"]
    else
        sensitized = Set{Int}()
    end

    # IF is the reference
    _price_objectives!(fid, dicts, "IF", scheme, price_if)

    # Rawls is the reference
    _price_objectives!(fid, dicts, "Rawls", scheme, price_rawls)

    # Aristotle is the reference
    _price_objectives!(fid, dicts, "Aristotle", scheme, (probs, ideal, nadir) -> price_aristotle(probs, ideal, nadir, sensitized))

    # Nash is the reference
    _price_objectives!(fid, dicts, "Nash", scheme, price_nash) 
end


function file_data(filename)
    fid = jldopen(filename, "r")
    single_dict = Dict{String, Float64}()
    linear_dict = Dict{String, Float64}()
    model_dict = Dict{String, Float64}()

    if haskey(fid, "stats/G")
        G = fid["stats/G"]
        single_dict["|P|"] = length(adjlistP(G))
        linear_dict["|P|"] = length(adjlistP(G))
        model_dict["|P|"] = length(adjlistP(G))

        single_dict["|N|"] = length(adjlistN(G))
        linear_dict["|N|"] = length(adjlistN(G))
        model_dict["|N|"] = length(adjlistN(G))
    else
        println(fid)
        single_dict["|P|"] = -1.0
        linear_dict["|P|"] = -1.0
        model_dict["|P|"] = -1.0

        single_dict["|N|"] = -1.0
        linear_dict["|N|"] = -1.0
        model_dict["|N|"] = -1.0
    end

    schemes = ["IF", "Rawls", "Nash", "Aristotle"]
    for scheme in schemes
        # get the time
        if haskey(fid, "stats/time/single/$scheme")
            time_single = fid["stats/time/single/$scheme"]
            single_dict["time/$scheme"] = time_single
        else
            single_dict["time/$scheme"] = -1.0
        end

        if haskey(fid, "stats/time/linear/$scheme")
            time_linear = fid["stats/time/linear/$scheme"]
            linear_dict["time/$scheme"] = time_linear + time_single
        else
            linear_dict["time/$scheme"] = -1.0
        end

        if haskey(fid, "stats/time/model/$scheme")
            time_model = fid["stats/time/model/$scheme"]
            model_dict["time/$scheme"] = time_model + time_single
        else
            model_dict["time/$scheme"] = -1.0
        end

        # get the distance to ideal
        if haskey(fid, "stats/ideal/$scheme")
            ideal = fid["stats/ideal/$scheme"]
            nadir = fid["stats/nadir/$scheme"]
            sol = (ideal[1], nadir[2])
            val = distance_to_ideal(sol, ideal, nadir) 
            single_dict["ideal_distance/$scheme"] = 0.0 <= val <= 1.0 ? val : -1.0
        else
            single_dict["ideal_distance/$scheme"] = -1.0
        end

        if haskey(fid, "stats/ideal_distance/linear/$scheme")
            ideal_dist_linear = fid["stats/ideal_distance/linear/$scheme"]
            linear_dict["ideal_distance/$scheme"] = 0.0 <= ideal_dist_linear <= 1.0 ? ideal_dist_linear : -1.0
        else
            linear_dict["ideal_distance/$scheme"] = -1.0
        end

        if haskey(fid, "stats/ideal_distance/model/$scheme")
            ideal_dist_model = fid["stats/ideal_distance/model/$scheme"]
            model_dict["ideal_distance/$scheme"] = 0.0 <= ideal_dist_model <= 1.0 ? ideal_dist_model : -1.0
        else
            model_dict["ideal_distance/$scheme"] = -1.0
        end

        # get the distance to nadir
        if haskey(fid, "stats/nadir/$scheme")
            ideal = fid["stats/ideal/$scheme"]
            nadir = fid["stats/nadir/$scheme"]
            sol = (ideal[1], nadir[2])
            val = distance_to_nadir(sol, ideal, nadir)
            single_dict["nadir_distance/$scheme"] = 0.0 <= val <= 1.0 ? val : -1.0
        else
            single_dict["nadir_distance/$scheme"] = -1.0
        end

        if haskey(fid, "stats/nadir_distance/linear/$scheme")
            nadir_dist_linear = fid["stats/nadir_distance/linear/$scheme"]
            linear_dict["nadir_distance/$scheme"] = 0.0 <= nadir_dist_linear <= 1.0 ? nadir_dist_linear : -1.0
        else
            linear_dict["nadir_distance/$scheme"] = -1.0
        end

        if haskey(fid, "stats/nadir_distance/model/$scheme")
            nadir_dist_model = fid["stats/nadir_distance/model/$scheme"]
            model_dict["nadir_distance/$scheme"] = 0.0 <= nadir_dist_model <= 1.0 ? nadir_dist_model : -1.0
        else
            model_dict["nadir_distance/$scheme"] = -1.0
        end

        # get the POF & POU
        price_objectives!(fid, (single_dict, linear_dict, model_dict), scheme)
        if haskey(fid, "stats/ideal/$scheme")
            ideal = fid["stats/ideal/$scheme"]
            nadir = fid["stats/nadir/$scheme"]
            numP = length(fid["stats/feasible"])
            pof_single = 0.0
            pou_single = price_utility(scheme, ideal, nadir, numP)
            
            single_dict["POF-$scheme/Utilitarian"] = pof_single
            single_dict["POU-$scheme/Utilitarian"] = pou_single
        else
            single_dict["POF-$scheme/Utilitarian"] = -1.0
            single_dict["POU-$scheme/Utilitarian"] = -1.0
        end


        # get the support size
        if haskey(fid, "stats/support_size/single/$scheme")
            support_size_model = fid["stats/support_size/single/$scheme"]
            single_dict["support_size/$scheme"] = support_size_model
        else
            single_dict["support_size/$scheme"] = -1.0
        end

        if haskey(fid, "stats/support_size/model/$scheme")
            support_size_model = fid["stats/support_size/model/$scheme"]
            model_dict["support_size/$scheme"] = support_size_model
        else
            model_dict["support_size/$scheme"] = -1.0
        end

        if haskey(fid, "stats/support_size/linear/$scheme")
            support_size_linear = fid["stats/support_size/linear/$scheme"]
            linear_dict["support_size/$scheme"] = support_size_linear
        else
            linear_dict["support_size/$scheme"] = -1.0
        end

    end

    close(fid)
    return single_dict, linear_dict, model_dict
end


function collect_stats(dir, output)
    # open file to store the data
    fid = h5open(output, "w") # TODO take care of this line
    
    single = Dict{String, Vector{Float64}}()
    linear = Dict{String, Vector{Float64}}()
    model = Dict{String, Vector{Float64}}()

    # loop over files and collect data
    for filename in Base.Filesystem.readdir(dir)
        if split(filename, '.', limit = 2)[end] != "mof.json"
            continue
        end
        # store the data in the file
        try
            single_dict, linear_dict, model_dict = file_data(joinpath(dir, filename))

            for (key, val) in single_dict 
                if ! haskey(single, key)
                    single[key] = Vector{Float64}()
                end
                push!(single[key], val)
            end

            for (key, val) in linear_dict
                if ! haskey(linear, key)
                    linear[key] = Vector{Float64}()
                end
                push!(linear[key], val)
            end

            for (key, val) in model_dict
                if ! haskey(model, key)
                    model[key] = Vector{Float64}()
                end
                push!(model[key], val)
            end
        catch e 
            if isa(e, EOFError)
                continue
            end
            throw(e)
        end
    end

    # save the experimental data for the model with linear combination of objectives
    for (key, val) in single 
        fid["single/$key"] = val
    end

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
