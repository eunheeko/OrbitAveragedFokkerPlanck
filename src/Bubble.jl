using Random
using Distributions
using Printf
using HDF5
using Logging
# using StaticArrays

include("Param.jl")

#=
Force fluctuations from SBs (3D)
=#
function cal_force_bubble(x1::Float64, y1::Float64, z1::Float64, x2::Float64, y2::Float64, z2::Float64, rbub_in::Float64, rbub_out::Float64, G::Float64, ρbub_in::Float64, ρbub_out::Float64)

    ## distance between two points
    r = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
    C = -4 * pi * G / 3

    ## density as a function of distance from the center
    if r <= rbub_in                        # density in inner cavity

        f = C * ρbub_in * r 

    elseif (r > rbub_in) & (r <= rbub_out) # density in outer shell

        f =  C * ( ρbub_in * rbub_in^3/r^2 + ρbub_out * (r - rbub_in^3/r^2) )

    else                                   # density outside of SB               
        f = 0.0
    end

    return f
end


#=
PDF of SB positions
=#
function pdf_r_inverselinear(b::Float64, rmax::Float64, rmin::Float64)
    
    ## PDF of SB positions
    fr = r -> 1/(1 + r/b)

    ## max(PDF)
    frmax = fr(rmin)

    return fr, frmax
end

#= 
rejection sampling for SB positions
=#
function mcsample_bubble_position(FSIZE::Float64, FSIZE_Z::Float64, Nsample::Int, fr::Function, frmax::Float64)
    
    ## initialize tables for saving sampled SB positions
    sample_x = fill(NaN, Nsample)
    sample_y = fill(NaN, Nsample)
    sample_z = fill(NaN, Nsample)

    ## sample n = Nsample SB positions
    for isample in 1:Nsample

        ## initialize threshold and positions
        u = 0.
        x = 0.
        y = 0.
        z = 0.

        ## iterate for rejection sampling
        while true
            ## sample positions following uniform PDF between -1 and 1
            x = rand(Uniform(  -FSIZE,    FSIZE))
            y = rand(Uniform(  -FSIZE,    FSIZE))
            z = rand(Uniform(-FSIZE_Z,  FSIZE_Z))
            
            ## distance from the galactic center (mid-plane)
            r = sqrt(x^2 + y^2)

            ## acceptance test
            u = rand() # threshold
            if u < fr(r) / frmax # reject condition with uniform control PDF
                break
            end
        end
        
        ## assign sampled position
        sample_x[isample] = x
        sample_y[isample] = y
        sample_z[isample] = z
    end

    return sample_x, sample_y, sample_z
end


#= 
create force field data
=#
function create_bubble_size(params::BubbleParameters, seed::Int, path_data::String, fname::String, LRANGE::Vector, LRANGE_Z::Vector, TIMES::Vector)

    path_data_seed = joinpath(path_data, fname*"-SEED-"*@sprintf("%d", seed)*".h5")
    if check_isfile(path_data_seed; msg_exist = "\t> File already exists: seed = $seed. Skipping.", msg_notexist = "\t> Start creating seed = $seed.")
        return
    end

    ######### params
    ## constant
    G = params.G

    ## bubbles
    RBUB_IN     = params.RBUB_IN
    RBUB_OUT    = params.RBUB_OUT
    TBUB        = params.TBUB
    ρBUB_IN     = params.ρBUB_IN
    ρBUB_OUT    = params.ρBUB_OUT

    ## gaseous disk
    dnBUB       = params.dnBUB
    RDISK       = params.RDISK

    ## simulation setting
    NGRID       = params.NGRID
    NGRID_Z     = params.NGRID_Z
    DEL_TSIM    = params.DEL_TSIM
    TSIM        = params.TSIM

    ## create simulation gread
    XGRID = repeat(reshape(  LRANGE, NGRID,     1,       1),     1, NGRID, NGRID_Z) # [L]
    YGRID = repeat(reshape(  LRANGE,     1, NGRID,       1), NGRID,     1, NGRID_Z) # [L]
    ZGRID = repeat(reshape(LRANGE_Z,     1,     1, NGRID_Z), NGRID, NGRID,       1) # [L]
    RGRID = @. sqrt(XGRID^2 + YGRID^2) # distance from the center (mid-plane) [L]

    ## generate random seed
    rng = Xoshiro(seed)

    ## random position of bubbles ~ exposnential profile
    rmax = FSIZE # maximum distance from the center [L]
    rmin = 0.    # minimum distance from the center [L]

    ## load SB distribution
    fr, frmax = pdf_r_inverselinear(RDISK, rmax, rmin)

    ## total number of SBs
    Vol_box = (2 * FSIZE)^2 * (2 * FSIZE_Z)                                 # simulation box volume 
    Vol_eff = sum(@. 1 / (1 + RGRID/RDISK)) / (NGRID^2 * NGRID_Z) * Vol_box # effective volume weigthed by 1/ (1+r/b)
    NBUB = floor(Int, dnBUB * Vol_eff * TSIM)                               # total number of SBs

    ## sample random position of SBs
    xcenters, ycenters, zcenters = mcsample_bubble_position(FSIZE, FSIZE_Z, NBUB, fr, frmax) # [L]

    ## sample random birth time of SBs
    t_births = rand(rng, Uniform(0,TSIM), NBUB) # sample birth time       [T]
    t_deaths = t_births .+ TBUB;                # time when SB radius = 0 [T] 

    ## pre-allocation
    # mat_dΨdr = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in r
    # dΨdr     = zeros(NGRID, NGRID, NGRID_Z)                # force field in r at t

    mat_dΨdx = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in x
    dΨdx     = zeros(NGRID, NGRID, NGRID_Z)                # force field in x at t
    mat_dΨdy = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in y
    dΨdy     = zeros(NGRID, NGRID, NGRID_Z)                # force field in y at t
    mat_dΨdz = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in z
    dΨdz     = zeros(NGRID, NGRID, NGRID_Z)                # force field in z at t

    r_sizes = zeros(NBUB)                                  # initialize a table to save SB size [L]

    ## force field at each timestep
    @time @views for (it, tt) in enumerate(TIMES)

        ## growing SBs
        tmask_grow = (t_births .<= tt) .&& (tt .< t_births .+ TBUB/2)
        @. r_sizes[tmask_grow] += RBUB_OUT/TBUB*2*DEL_TSIM

        ## collapsing SBs
        tmask_collapse = (t_deaths .>= tt).&& (tt .>= t_births .+ TBUB/2)
        @. r_sizes[tmask_collapse] -= RBUB_OUT/TBUB*2*DEL_TSIM

        ## initialize a new force matrix at t.
        # @. dΨdr .= 0
        @. dΨdx .= 0
        @. dΨdy .= 0
        @. dΨdz .= 0

        ## mask to select (1) either growing or collapsing bubbles and (2) exclude negative or zero radius bubble
        tmask = (tmask_grow .| tmask_collapse) .& (r_sizes .> 0)
        
        ## calculate force for alive SBs only
        @inbounds for idx in findall(tmask)

            ## load position of alive SB
            xx = xcenters[idx] # [L]
            yy = ycenters[idx] # [L]
            zz = zcenters[idx] # [L]
            
            ## load size of alive SBs
            r_out = r_sizes[idx]              # outer shell size  [L]
            r_in = r_out * RBUB_IN / RBUB_OUT # inner cavity size [L]

            ## calculate the force field by SBs
            # dΨdr .+= cal_force_bubble.(Ref(xx), Ref(yy), Ref(zz), XGRID, YGRID, ZGRID, r_in, r_out, G, ρBUB_IN, ρBUB_OUT)
            
            _f = cal_force_bubble.(Ref(xx), Ref(yy), Ref(zz), XGRID, YGRID, ZGRID, r_in, r_out, G, ρBUB_IN, ρBUB_OUT)
            
            DX = XGRID .- xx
            DY = YGRID .- yy
            DZ = ZGRID .- zz
            DR = sqrt.(DX.^2 .+ DY.^2 .+ DZ.^2)

            dΨdx .+= _f .* DX ./ DR
            dΨdy .+= _f .* DY ./ DR
            dΨdz .+= _f .* DZ ./ DR
            # dΨdr .+= _f
        end

        ## update force field at t
        # @. mat_dΨdr[:, :, :, it] = dΨdr

        @. mat_dΨdx[:, :, :, it] = dΨdx
        @. mat_dΨdy[:, :, :, it] = dΨdy
        @. mat_dΨdz[:, :, :, it] = dΨdz

    end 

    # # save force field
    # path_save = joinpath(path_data, "dphidr-SEED-"*@sprintf("%d", seed)*".h5")
    # h5open(path_save, "w") do file
    #     write(file, "data", mat_dΨdr)
    # end

    path_save = joinpath(path_data, "dphidx-SEED-"*@sprintf("%d", seed)*".h5")
    h5open(path_save, "w") do file
        write(file, "data", mat_dΨdx)
    end

    path_save = joinpath(path_data, "dphidy-SEED-"*@sprintf("%d", seed)*".h5")
    h5open(path_save, "w") do file
        write(file, "data", mat_dΨdy)
    end

    path_save = joinpath(path_data, "dphidz-SEED-"*@sprintf("%d", seed)*".h5")
    h5open(path_save, "w") do file
        write(file, "data", mat_dΨdz)
    end


    @info "\t> 3D-disk dϕdr is and saved: seed = $seed"

end

#= 
wrapper
=#
function run_create_bubble_size(params::BubbleParameters, path_data::String, fname::String, seed_start::Int, seed_end::Int)

    ## load simulation patameters
    FSIZE    = params.FSIZE
    FSIZE_Z  = params.FSIZE_Z
    NGRID    = params.NGRID
    NGRID_Z  = params.NGRID_Z
    DEL_TSIM = params.DEL_TSIM
    TSIM     = params.TSIM

    ## create space and time grid
    LRANGE   = vcat(LinRange(-FSIZE,FSIZE,NGRID))       # spatial grid (x, y) [L]
    LRANGE_Z = vcat(LinRange(-FSIZE_Z,FSIZE_Z,NGRID_Z)) # spatial grid (z)    [L]
    TIMES    = collect(0:DEL_TSIM:TSIM)                 # time grid [T]

    ## path for saving realizations
    ensure_mkpath(path_data)

    ## run
    Threads.@threads for seed in seed_start:seed_end
        create_bubble_size(params, seed, path_data, fname, LRANGE, LRANGE_Z, TIMES)
    end

end


#= 
create force field data
=#
function create_bubble_size_inplane(params::BubbleParameters, seed::Int, path_data::String, fname::String, LRANGE::Vector, LRANGE_Z::Vector, TIMES::Vector)

    path_data_seed = joinpath(path_data, fname*"-SEED-"*@sprintf("%d", seed)*".h5")
    if check_isfile(path_data_seed; msg_exist = "\t> File already exists: seed = $seed. Skipping.", msg_notexist = "\t> Start creating seed = $seed.")
        return
    end

    ######### params
    ## constant
    G = params.G

    ## bubbles
    RBUB_IN     = params.RBUB_IN
    RBUB_OUT    = params.RBUB_OUT
    TBUB        = params.TBUB
    ρBUB_IN     = params.ρBUB_IN
    ρBUB_OUT    = params.ρBUB_OUT

    ## gaseous disk
    dnBUB       = params.dnBUB
    RDISK       = params.RDISK

    ## simulation setting
    NGRID       = params.NGRID
    NGRID_Z     = params.NGRID_Z
    DEL_TSIM    = params.DEL_TSIM
    TSIM        = params.TSIM

    ## create simulation gread
    XGRID = repeat(reshape(  LRANGE, NGRID,     1,       1),     1, NGRID, NGRID_Z) # [L]
    YGRID = repeat(reshape(  LRANGE,     1, NGRID,       1), NGRID,     1, NGRID_Z) # [L]
    ZGRID = repeat(reshape(LRANGE_Z,     1,     1, NGRID_Z), NGRID, NGRID,       1) # [L]
    RGRID = @. sqrt(XGRID^2 + YGRID^2) # distance from the center (mid-plane) [L]

    ## generate random seed
    rng = Xoshiro(seed)

    ## random position of bubbles ~ exposnential profile
    rmax = FSIZE # maximum distance from the center [L]
    rmin = 0.    # minimum distance from the center [L]

    ## load SB distribution
    fr, frmax = pdf_r_inverselinear(RDISK, rmax, rmin)

    ## total number of SBs
    Vol_box = (2 * FSIZE)^2 * (2 * FSIZE_Z)                                 # simulation box volume 
    Vol_eff = sum(@. 1 / (1 + RGRID/RDISK)) / (NGRID^2 * NGRID_Z) * Vol_box # effective volume weigthed by 1/ (1+r/b)
    NBUB = floor(Int, dnBUB * Vol_eff * TSIM)                               # total number of SBs

    ## sample random position of SBs
    xcenters, ycenters, zcenters = mcsample_bubble_position(FSIZE, FSIZE_Z, NBUB, fr, frmax) # [L]
    zcenters .= 0. # set z = 0 for in-plane bubbles

    ## sample random birth time of SBs
    t_births = rand(rng, Uniform(0,TSIM), NBUB) # sample birth time       [T]
    t_deaths = t_births .+ TBUB;                # time when SB radius = 0 [T] 

    ## pre-allocation
    # mat_dΨdr = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in r
    # dΨdr     = zeros(NGRID, NGRID, NGRID_Z)                # force field in r at t

    mat_dΨdx = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in x
    dΨdx     = zeros(NGRID, NGRID, NGRID_Z)                # force field in x at t
    mat_dΨdy = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in y
    dΨdy     = zeros(NGRID, NGRID, NGRID_Z)                # force field in y at t
    mat_dΨdz = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in z
    dΨdz     = zeros(NGRID, NGRID, NGRID_Z)                # force field in z at t

    r_sizes = zeros(NBUB)                                  # initialize a table to save SB size [L]

    ## force field at each timestep
    @time @views for (it, tt) in enumerate(TIMES)

        ## growing SBs
        tmask_grow = (t_births .<= tt) .&& (tt .< t_births .+ TBUB/2)
        @. r_sizes[tmask_grow] += RBUB_OUT/TBUB*2*DEL_TSIM

        ## collapsing SBs
        tmask_collapse = (t_deaths .>= tt).&& (tt .>= t_births .+ TBUB/2)
        @. r_sizes[tmask_collapse] -= RBUB_OUT/TBUB*2*DEL_TSIM

        ## initialize a new force matrix at t.
        # @. dΨdr .= 0
        @. dΨdx .= 0
        @. dΨdy .= 0
        @. dΨdz .= 0

        ## mask to select (1) either growing or collapsing bubbles and (2) exclude negative or zero radius bubble
        tmask = (tmask_grow .| tmask_collapse) .& (r_sizes .> 0)
        
        ## calculate force for alive SBs only
        @inbounds for idx in findall(tmask)

            ## load position of alive SB
            xx = xcenters[idx] # [L]
            yy = ycenters[idx] # [L]
            zz = zcenters[idx] # [L]
            
            ## load size of alive SBs
            r_out = r_sizes[idx]              # outer shell size  [L]
            r_in = r_out * RBUB_IN / RBUB_OUT # inner cavity size [L]

            ## calculate the force field by SBs
            # dΨdr .+= cal_force_bubble.(Ref(xx), Ref(yy), Ref(zz), XGRID, YGRID, ZGRID, r_in, r_out, G, ρBUB_IN, ρBUB_OUT)
            
            _f = cal_force_bubble.(Ref(xx), Ref(yy), Ref(zz), XGRID, YGRID, ZGRID, r_in, r_out, G, ρBUB_IN, ρBUB_OUT)
            
            DX = XGRID .- xx
            DY = YGRID .- yy
            DZ = ZGRID .- zz
            DR = sqrt.(DX.^2 .+ DY.^2 .+ DZ.^2)

            dΨdx .+= _f .* DX ./ DR
            dΨdy .+= _f .* DY ./ DR
            dΨdz .+= _f .* DZ ./ DR
            # dΨdr .+= _f
        end

        ## update force field at t
        # @. mat_dΨdr[:, :, :, it] = dΨdr

        @. mat_dΨdx[:, :, :, it] = dΨdx
        @. mat_dΨdy[:, :, :, it] = dΨdy
        @. mat_dΨdz[:, :, :, it] = dΨdz

    end 

    # # save force field
    # path_save = joinpath(path_data, "dphidr-SEED-"*@sprintf("%d", seed)*".h5")
    # h5open(path_save, "w") do file
    #     write(file, "data", mat_dΨdr)
    # end

    path_save = joinpath(path_data, "dphidx-SEED-"*@sprintf("%d", seed)*".h5")
    h5open(path_save, "w") do file
        write(file, "data", mat_dΨdx)
    end

    path_save = joinpath(path_data, "dphidy-SEED-"*@sprintf("%d", seed)*".h5")
    h5open(path_save, "w") do file
        write(file, "data", mat_dΨdy)
    end

    path_save = joinpath(path_data, "dphidz-SEED-"*@sprintf("%d", seed)*".h5")
    h5open(path_save, "w") do file
        write(file, "data", mat_dΨdz)
    end


    @info "\t> 3D-disk dϕdr is and saved: seed = $seed"

end

#= 
wrapper
=#
function run_create_bubble_size_inplane(params::BubbleParameters, path_data::String, fname::String, seed_start::Int, seed_end::Int)

    ## load simulation patameters
    FSIZE    = params.FSIZE
    FSIZE_Z  = params.FSIZE_Z
    NGRID    = params.NGRID
    NGRID_Z  = params.NGRID_Z
    DEL_TSIM = params.DEL_TSIM
    TSIM     = params.TSIM

    ## create space and time grid
    LRANGE   = vcat(LinRange(-FSIZE,FSIZE,NGRID))       # spatial grid (x, y) [L]
    LRANGE_Z = vcat(LinRange(-FSIZE_Z,FSIZE_Z,NGRID_Z)) # spatial grid (z)    [L]
    TIMES    = collect(0:DEL_TSIM:TSIM)                 # time grid [T]

    ## path for saving realizations
    ensure_mkpath(path_data)

    ## run
    Threads.@threads for seed in seed_start:seed_end
        create_bubble_size_inplane(params, seed, path_data, fname, LRANGE, LRANGE_Z, TIMES)
    end

end


# #= 
# create force field data
# =#
# function create_bubble_size_periodic(params::BubbleParameters, seed::Int, path_data::String, fname::String, LRANGE::Vector, LRANGE_Z::Vector, TIMES::Vector)

#     path_data_seed = joinpath(path_data, fname*"-SEED-"*@sprintf("%d", seed)*".h5")
#     if check_isfile(path_data_seed; msg_exist = "\t> File already exists: seed = $seed. Skipping.", msg_notexist = "\t> Start creating seed = $seed.")
#         return
#     end

#     ######### params
#     ## constant
#     G = params.G

#     ## bubbles
#     RBUB_IN     = params.RBUB_IN
#     RBUB_OUT    = params.RBUB_OUT
#     TBUB        = params.TBUB
#     ρBUB_IN     = params.ρBUB_IN
#     ρBUB_OUT    = params.ρBUB_OUT

#     ## gaseous disk
#     dnBUB       = params.dnBUB
#     RDISK       = params.RDISK

#     ## simulation setting
#     NGRID       = params.NGRID
#     NGRID_Z     = params.NGRID_Z
#     DEL_TSIM    = params.DEL_TSIM
#     TSIM        = params.TSIM

#     ## create simulation gread
#     XGRID = repeat(reshape(  LRANGE, NGRID,     1,       1),     1, NGRID, NGRID_Z) # [L]
#     YGRID = repeat(reshape(  LRANGE,     1, NGRID,       1), NGRID,     1, NGRID_Z) # [L]
#     ZGRID = repeat(reshape(LRANGE_Z,     1,     1, NGRID_Z), NGRID, NGRID,       1) # [L]
#     RGRID = @. sqrt(XGRID^2 + YGRID^2) # distance from the center (mid-plane) [L]

#     ## generate random seed
#     rng = Xoshiro(seed)

#     ## random position of bubbles ~ exposnential profile
#     rmax = FSIZE # maximum distance from the center [L]
#     rmin = 0.    # minimum distance from the center [L]

#     ## load SB distribution
#     # fr, frmax = pdf_r_inverselinear(RDISK, rmax, rmin)

#     ## total number of SBs
#     # Vol_box = (2 * FSIZE)^2 * (2 * FSIZE_Z)                                 # simulation box volume 
#     # Vol_eff = sum(@. 1 / (1 + RGRID/RDISK)) / (NGRID^2 * NGRID_Z) * Vol_box # effective volume weigthed by 1/ (1+r/b)
#     # NBUB = floor(Int, dnBUB * Vol_eff * TSIM)                               # total number of SBs
#     NBUB = 1

#     ## sample random position of SBs
#     # xcenters, ycenters, zcenters = mcsample_bubble_position(FSIZE, FSIZE_Z, NBUB, fr, frmax) # [L]
#     xcenters = zeros(NBUB)
#     ycenters = zeros(NBUB)
#     zcenters = zeros(NBUB) # rand(rng, Uniform(-FSIZE_Z, FSIZE_Z), NBUB)

#     ## sample random birth time of SBs
#     # t_births = rand(rng, Uniform(0,TSIM), NBUB) # zeros(NBUB) # rand(rng, Uniform(0,TSIM), NBUB)  # rand(rng, Uniform(0,TSIM/2), NBUB) # sample birth time       [T]
#     t_births = zeros(NBUB) # zeros(NBUB) # rand(rng, Uniform(0,TSIM), NBUB)  # rand(rng, Uniform(0,TSIM/2), NBUB) # sample birth time       [T]
#     t_deaths = t_births .+ TBUB;                # time when SB radius = 0 [T] 

#     ## pre-allocation
#     mat_dΨdr = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in r
#     dΨdr     = zeros(NGRID, NGRID, NGRID_Z)                # force field in r at t

#     mat_dΨdx = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in x
#     dΨdx     = zeros(NGRID, NGRID, NGRID_Z)                # force field in x at t
#     mat_dΨdy = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in y
#     dΨdy     = zeros(NGRID, NGRID, NGRID_Z)                # force field in y at t
#     mat_dΨdz = zeros(NGRID, NGRID, NGRID_Z, length(TIMES)) # force field in z
#     dΨdz     = zeros(NGRID, NGRID, NGRID_Z)                # force field in z at t

#     r_sizes = zeros(NBUB)                                  # initialize a table to save SB size [L]
    
#     ## force field at each timestep
#     @time @views for (it, tt) in enumerate(TIMES)

#         ## initialize a new force matrix at t.
#         # @. dΨdr .= 0
#         # @. dΨdx .= 0
#         # @. dΨdy .= 0
#         # @. dΨdz .= 0


#         xx = 0.
#         yy = 0.
#         zz = 0.

#         _f = cal_force_bubble.(Ref(xx), Ref(yy), Ref(zz), XGRID, YGRID, ZGRID, RBUB_IN, RBUB_OUT, G, ρBUB_IN, ρBUB_OUT)
        
#         DX = XGRID .- xx
#         DY = YGRID .- yy
#         DZ = ZGRID .- zz
        
#         DR = sqrt.(DX.^2 .+ DY.^2 .+ DZ.^2)

#         # dΨdx .+= _f .* DX ./ DR
#         # dΨdy .+= _f .* DY ./ DR
#         # dΨdz .+= _f .* DZ ./ DR
#         # dΨdr .+= _f

#         ## update force field at t
#         @. mat_dΨdr[:, :, :, it] = _f

#         @. mat_dΨdx[:, :, :, it] .= _f .* DX ./ DR
#         @. mat_dΨdy[:, :, :, it] .= _f .* DY ./ DR
#         @. mat_dΨdz[:, :, :, it] .= _f .* DZ ./ DR

#     end 

#     # ## force field at each timestep
#     # @time @views for (it, tt) in enumerate(TIMES)

#     #     ## growing SBs
#     #     tmask_grow = (t_births .<= tt) .&& (tt .< t_births .+ TBUB/2)
#     #     @. r_sizes[tmask_grow] += RBUB_OUT/TBUB*2*DEL_TSIM

#     #     ## collapsing SBs
#     #     tmask_collapse = (t_deaths .>= tt).&& (tt .>= t_births .+ TBUB/2)
#     #     @. r_sizes[tmask_collapse] -= RBUB_OUT/TBUB*2*DEL_TSIM

#     #     ## initialize a new force matrix at t.
#     #     @. dΨdr .= 0
#     #     @. dΨdx .= 0
#     #     @. dΨdy .= 0
#     #     @. dΨdz .= 0

#     #     ## mask to select (1) either growing or collapsing bubbles and (2) exclude negative or zero radius bubble
#     #     tmask = (tmask_grow .| tmask_collapse) .& (r_sizes .> 0)
        
#     #     ## calculate force for alive SBs only
#     #     @inbounds for idx in findall(tmask)

#     #         ## load position of alive SB
#     #         xx = xcenters[idx] # [L]
#     #         yy = ycenters[idx] # [L]
#     #         zz = zcenters[idx] # [L]
            
#     #         ## load size of alive SBs
#     #         r_out = r_sizes[idx]              # outer shell size  [L]
#     #         r_in = r_out * RBUB_IN / RBUB_OUT # inner cavity size [L]

#     #         ## calculate the force field by SBs
#     #         # dΨdr .+= cal_force_bubble.(Ref(xx), Ref(yy), Ref(zz), XGRID, YGRID, ZGRID, r_in, r_out, G, ρBUB_IN, ρBUB_OUT)
            
#     #         _f = cal_force_bubble.(Ref(xx), Ref(yy), Ref(zz), XGRID, YGRID, ZGRID, r_in, r_out, G, ρBUB_IN, ρBUB_OUT)
            
#     #         DX = XGRID .- xx
#     #         DY = YGRID .- yy
#     #         DZ = ZGRID .- zz
#     #         DR = sqrt.(DX.^2 .+ DY.^2 .+ DZ.^2)

#     #         dΨdx .+= _f .* DX ./ DR
#     #         dΨdy .+= _f .* DY ./ DR
#     #         dΨdz .+= _f .* DZ ./ DR
#     #         dΨdr .+= _f
#     #     end

#     #     ## update force field at t
#     #     @. mat_dΨdr[:, :, :, it] = dΨdr

#     #     @. mat_dΨdx[:, :, :, it] = dΨdx
#     #     @. mat_dΨdy[:, :, :, it] = dΨdy
#     #     @. mat_dΨdz[:, :, :, it] = dΨdz

#     # end 

#     # save force field
#     path_save = joinpath(path_data, "dphidr-SEED-"*@sprintf("%d", seed)*".h5")
#     h5open(path_save, "w") do file
#         write(file, "data", mat_dΨdr)
#     end

#     path_save = joinpath(path_data, "dphidx-SEED-"*@sprintf("%d", seed)*".h5")
#     h5open(path_save, "w") do file
#         write(file, "data", mat_dΨdx)
#     end

#     path_save = joinpath(path_data, "dphidy-SEED-"*@sprintf("%d", seed)*".h5")
#     h5open(path_save, "w") do file
#         write(file, "data", mat_dΨdy)
#     end

#     path_save = joinpath(path_data, "dphidz-SEED-"*@sprintf("%d", seed)*".h5")
#     h5open(path_save, "w") do file
#         write(file, "data", mat_dΨdz)
#     end


#     @info "\t> 3D-disk dϕdr is and saved: seed = $seed"

# end

# #= 
# wrapper
# =#
# function run_create_bubble_size_periodic(params::BubbleParameters, path_data::String, fname::String, seed_start::Int, seed_end::Int)

#     ## load simulation patameters
#     FSIZE    = params.FSIZE
#     FSIZE_Z  = params.FSIZE_Z
#     NGRID    = params.NGRID
#     NGRID_Z  = params.NGRID_Z
#     DEL_TSIM = params.DEL_TSIM
#     TSIM     = params.TSIM

#     ## create space and time grid
#     LRANGE   = vcat(LinRange(-FSIZE,FSIZE,NGRID))       # spatial grid (x, y) [L]
#     LRANGE_Z = vcat(LinRange(-FSIZE_Z,FSIZE_Z,NGRID_Z)) # spatial grid (z)    [L]
#     TIMES    = collect(0:DEL_TSIM:TSIM)                 # time grid [T]

#     ## path for saving realizations
#     ensure_mkpath(path_data)

#     ## run
#     Threads.@threads for seed in seed_start:seed_end
#         create_bubble_size_periodic(params, seed, path_data, fname, LRANGE, LRANGE_Z, TIMES)
#     end

# end
