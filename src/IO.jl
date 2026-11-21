using JSON
using Logging
using Printf

#=
create directory name
=#
function create_dirname_ic(type_orbit::String, E::Float64, L::Float64)

    ## generate icfile path name
    if type_orbit == "2D"
        icname = "E_"*@sprintf("%.3f", E)*"_L_"*@sprintf("%.3f", L)
        return icname
    else
        icname = "E_"*@sprintf("%.3f", E)
        return icname
    end
end

"""
    ensure_mkpath(pname::String; msg_create::String="", msg_exists::String="")

Ensure that the given path exists.  
If it does not exist, create it using `mkpath(pname)`.

# Keyword arguments
- `msg_create`: optional message logged when a new path is created  
- `msg_exists`: optional message logged when the path already exists
"""
function ensure_mkpath(pname::String;
     msg_create::String = "\t> Created path: $pname", 
     msg_exist::String  = "\t> Path already exists: $pname")

    if !ispath(pname)
        mkpath(pname)
        @debug msg_create
    else
        @debug msg_exist
    end
end

"""
    check_isfile(pname::String; msg_exist::String="", msg_notexist::String="")

Check whether a file exists at the given path.

# Keyword arguments
- `msg_exist`: message logged when the file exists
- `msg_notexist`: message logged when the file does not exist

# Returns
- `true` if the file exists
- `false` otherwise
"""
function check_isfile(pname::String;
        msg_exist::String    = "\t> File already exists: $pname",
        msg_notexist::String = "\t> File does not exist: $pname")

    if isfile(pname)
        @debug msg_exist
        return true
    else
        @debug msg_notexist
        return false
    end
end

"""
    IO functions
"""

"""
    struct_to_dict(x)

Convert a Julia struct to a dictionary mapping field names to their values.
"""
function struct_to_dict(x)
    return Dict(Symbol(k) => getfield(x, Symbol(k)) for k in fieldnames(typeof(x)))
end


#=
save parameter file
=#
function save_params(params, path_data::String, fname::String)
    
    ## ensure directory exists
    ensure_mkpath(path_data)

    ## parameter file
    json_path = joinpath(path_data, fname*".json")	 
    
    ## save parameter file
    open(json_path, "w") do file
        write(file, JSON.json(struct_to_dict(params)))
    end
    @info "\t> Parameters saved to $json_path"
end

#=
load parameter file
=#
function load_params(::Type{T}, path_data::String, fname::String) where T
    
    ## paramter file
    json_path = joinpath(path_data, fname*".json")

    ## load parameter file
    data = JSON.parsefile(json_path)

    return T(;(Symbol(k) => v for (k, v) in data)...)
end
