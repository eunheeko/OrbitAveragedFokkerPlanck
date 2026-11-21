using NaNMath
using Interpolations
using Printf
using HDF5

include("Param.jl")

struct BubbleInterpolation
    itp_dphidx::Interpolations.GriddedInterpolation
    itp_dphidy::Interpolations.GriddedInterpolation
    itp_dphidz::Interpolations.GriddedInterpolation
end

#=
Load interpolation function of SB force fluctuation fields
=#
function load_bubble_interpolation(path_data::String, seed::Int, xgrid::Array, ygrid::Array, zgrid::Array, rgrid::Array, lrange::Vector, lrange_z::Vector, times::Vector) 


    # ## load SB force fluctuation fields (discrete space and time)
    # dphidr = h5open(path_data*"/dphidr/dphidr-SEED-"*@sprintf("%d", seed)*".h5", "r") do file
    #     read(file, "data")
    # end

    # ## SB force fluctuation field in (x, y, z)
    # dphidx = @. xgrid / rgrid * dphidr
    # dphidy = @. ygrid / rgrid * dphidr
    # dphidz = @. zgrid / rgrid * dphidr

    ## load SB force fluctuation fields (discrete space and time)
    dphidx = h5open(path_data*"/dphidr/dphidx-SEED-"*@sprintf("%d", seed)*".h5", "r") do file
        read(file, "data")
    end

    dphidy = h5open(path_data*"/dphidr/dphidy-SEED-"*@sprintf("%d", seed)*".h5", "r") do file
        read(file, "data")
    end

    dphidz = h5open(path_data*"/dphidr/dphidz-SEED-"*@sprintf("%d", seed)*".h5", "r") do file
        read(file, "data")
    end

    ## interpolated SB fluctuation field
    itp_dphidx =  interpolate((lrange, lrange, lrange_z, times), dphidx, Gridded(Linear()))
    itp_dphidy =  interpolate((lrange, lrange, lrange_z, times), dphidy, Gridded(Linear()))
    itp_dphidz =  interpolate((lrange, lrange, lrange_z, times), dphidz, Gridded(Linear()))

    return BubbleInterpolation(itp_dphidx, itp_dphidy, itp_dphidz)
end

#=
force by mean field halo potential
=#
function dΨ_halo!(pos::Vector, force_halo::Vector) 
    @views force_halo .= dΨ_cartesian(pos)
end

#=
force fluctuations by SBs
=#
function dΨ_bubble!(pos::Vector, force_bubble::Vector, t::Float64, fsize::Float64, fsize_z::Float64, ρfactor::Float64, itp_bubble::BubbleInterpolation, margin::Float64=1e-10)
    
    ## load position of a test particle
    x, y, z = pos

    ## load force fluctuations by SBs
    if (x <= -fsize + margin) || (x >= fsize - margin) || (y <= -fsize + margin) || (y >= fsize - margin) || (z <= -fsize_z + margin) || (z >= fsize_z - margin)
        
        ## no force outside of the simulation box
        force_bubble .= 0.0
    else

        force_bubble[1] = itp_bubble.itp_dphidx(x, y, z, t) * ρfactor
        force_bubble[2] = itp_bubble.itp_dphidy(x, y, z, t) * ρfactor
        force_bubble[3] = itp_bubble.itp_dphidz(x, y, z, t) * ρfactor
    end

end

#=
drift function for leapfrog
=#
function drift!(pos::Vector, vel::Vector, dt::Float64)
    @. pos += vel * dt
end

#=
kick function for leapfrog (non-Markovian)
=#
function kick!(pos::Vector, vel::Vector, force_bubble::Vector, force_halo::Vector, dt::Float64, t::Float64, fsize::Float64, fsize_z::Float64, ρfactor::Float64, itp_bubble::BubbleInterpolation)

    ## force by mean potential + SB force fluctuations
    dΨ_halo!(pos, force_halo)
    dΨ_bubble!(pos, force_bubble, t, fsize, fsize_z, ρfactor, itp_bubble)
    
    @. vel -= (force_halo + force_bubble) * dt

end

#=
kick function for leapfrog (Markovian)
=#
function kick_nobub!(pos::Vector, vel::Vector, force_halo::Vector, dt::Float64)

    ## force by mean potential
    dΨ_halo!(pos, force_halo)

    @. vel -= (force_halo) * dt

end

#=
kick function for leapfrog (periodic orbit)
=#
function kick_periodic!(pos::Vector, vel::Vector, force_bubble::Vector, force_halo::Vector, dt::Float64, t::Float64, fsize::Float64, fsize_z::Float64, ρfactor::Float64, itp_bubble::BubbleInterpolation, flag_markovian::Bool)

    if flag_markovian # Markovian process, no bubble force
        return nothing
    else              # non-Markovian process, with bubble force
        dΨ_bubble!(pos, force_bubble, t, fsize, fsize_z, ρfactor, itp_bubble)
        @. vel -= (force_bubble) * dt
    end

end

#=
relocate the position of test particle when leaving the simulation box (periodic orbit)
=#
function wrap_periodic!(pos, fsize::Float64, fsize_z::Float64)

    x, y, z = pos
    if x < -fsize
        pos[1] += 2 * fsize
    elseif x > fsize
        pos[1] -= 2 * fsize
    end

    if y < -fsize
        pos[2] += 2 * fsize
    elseif y > fsize
        pos[2] -= 2 * fsize
    end

    if z < -fsize_z
        pos[3] += 2 * fsize_z
    elseif z > fsize_z
        pos[3] -= 2 * fsize_z
    end

end

#=
wrapper function for orbit integration (non-Markovian)
=#
function leapfrog!(positions::Matrix, velocities::Matrix, pos::Vector, vel::Vector, leapfrog_param::LeapfrogParameter, itp_bubble::BubbleInterpolation, dΨ_cartesian::Function)    

    ## load leapfrog parameters
    dt      = leapfrog_param.dt
    tsim    = leapfrog_param.tsim
    times   = leapfrog_param.times
    nstep   = leapfrog_param.nstep
    nsave   = leapfrog_param.nsave
    ρfactor = leapfrog_param.ρfactor

    ## load simulation parameters
    fsize   = leapfrog_param.fsize
    fsize_z = leapfrog_param.fsize_z

    ## pre-allocation for intermediate calculations
    force_bubble = [0., 0., 0.] # SB force fluctuations at t
    force_halo   = [0., 0., 0.] # mean field force at t
    dt_half      = dt/2.        # half timestep

    ## initialize the initial position and velocities
    @views positions[1, :]  .= copy(pos)
    @views velocities[1, :] .= copy(vel)

    ## orbit integration
    for it in 1:nstep

        ## DKD scheme
        drift!(pos, vel, dt_half)
        kick!(pos, vel, force_bubble, force_halo, dt, times[it], fsize, fsize_z, ρfactor, itp_bubble)
        drift!(pos, vel, dt_half)

        ## save orbits for every nsave points 
        if mod(it, nsave) == 0
            idx_nsave = div(it, nsave) + 1

            @views positions[idx_nsave, :]  .= pos
            @views velocities[idx_nsave, :] .= vel
        end

    end
end


#=
wrapper function for orbit integration (Markovian)
=#
function leapfrog_nobub!(positions::Matrix, velocities::Matrix, pos::Vector, vel::Vector, leapfrog_param::LeapfrogParameter, dΨ_cartesian::Function)    

    ## load leapfrog parameters
    dt      = leapfrog_param.dt
    tsim    = leapfrog_param.tsim
    times   = leapfrog_param.times
    nstep   = leapfrog_param.nstep
    nsave   = leapfrog_param.nsave
    ρfactor = leapfrog_param.ρfactor


    ## load simulation parameters
    fsize   = leapfrog_param.fsize
    fsize_z = leapfrog_param.fsize_z

    ## pre-allocation for intermediate calculations
    force_halo   = [0., 0., 0.] # mean field force at t
    dt_half      = dt/2.        # half timestep

    ## initialize the initial position and velocities
    @views positions[1, :]  .= copy(pos)
    @views velocities[1, :] .= copy(vel)

    # add half-kick velocity to the initial velocity
    # dΨ_halo!(pos, force_halo) # upate the force by halo potential
    # @. vel -= (force_halo) * dt_half # update the half kick velocity to be symplectic, v(t + dt/2)

    ## orbit integration
    for it in 1:nstep

        ## DKD scheme
        drift!(pos, vel, dt_half)             # after half drift: x(t + dt/2), v(t)
        kick_nobub!(pos, vel, force_halo, dt) # after full kick: x(t + dt/2), v(t + dt)
        drift!(pos, vel, dt_half)             # after half drift: x(t + dt), v(t + dt)

        ## save orbits for every nsave points 
        if mod(it, nsave) == 0
            idx_nsave = div(it, nsave) + 1 

            @views positions[idx_nsave, :] .= pos
            @views velocities[idx_nsave, :] .= vel
        end

    end
end


#=
wrapper function for orbit integration (non-Markovian)
=#
function leapfrog_periodic!(positions::Matrix, velocities::Matrix, pos::Vector, vel::Vector, leapfrog_param::LeapfrogParameter, itp_bubble::BubbleInterpolation, Ψ::Function, dΨ::Function, dΨ_cartesian::Function, flag_markovian::Bool)    

    ## load leapfrog parameters
    dt      = leapfrog_param.dt
    tsim    = leapfrog_param.tsim
    times   = leapfrog_param.times
    nstep   = leapfrog_param.nstep
    nsave   = leapfrog_param.nsave
    ρfactor = leapfrog_param.ρfactor

    ## load simulation parameters
    fsize   = leapfrog_param.fsize
    fsize_z = leapfrog_param.fsize_z

    ## pre-allocation for intermediate calculations
    force_bubble = [0., 0., 0.] # SB force fluctuations at t
    force_halo   = [0., 0., 0.] # mean field force at t
    dt_half      = dt/2.        # half timestep

    ## initialize the initial position and velocities
    @views positions[1, :]  .= copy(pos)
    @views velocities[1, :] .= copy(vel)

    ## orbit integration
    for it in 1:nstep
        drift!(pos, vel, dt_half)           # after half drift: x(t + dt/2), v(t)
        # wrap_periodic!(pos, fsize, fsize_z) # wrap the position to be periodic before calculating the force
        kick_periodic!(pos, vel, force_bubble, force_halo, dt, times[it], fsize, fsize_z, ρfactor, itp_bubble, flag_markovian) # after full kick: x(t + dt/2), v(t + dt)
        drift!(pos, vel, dt_half)           # after half drift: x(t + dt), v(t + dt)
        # wrap_periodic!(pos, fsize, fsize_z) # wrap the position to be periodic before calculating the force

        ## save orbits for every nsave points 
        if mod(it, nsave) == 0
            idx_nsave = div(it, nsave) + 1 

            @views positions[idx_nsave, :] .= pos
            @views velocities[idx_nsave, :] .= vel
        end
    end
end

#=
wrapper function for orbit integration (Markovian)
=#
function leapfrog_periodic_nobub!(positions::Matrix, velocities::Matrix, pos::Vector, vel::Vector, leapfrog_param::LeapfrogParameter, Ψ::Function, dΨ::Function, dΨ_cartesian::Function, flag_markovian::Bool)    

    ## load leapfrog parameters
    dt      = leapfrog_param.dt
    tsim    = leapfrog_param.tsim
    times   = leapfrog_param.times
    nstep   = leapfrog_param.nstep
    nsave   = leapfrog_param.nsave
    ρfactor = leapfrog_param.ρfactor


    ## load simulation parameters
    fsize   = leapfrog_param.fsize
    fsize_z = leapfrog_param.fsize_z

    ## pre-allocation for intermediate calculations
    force_halo   = [0., 0., 0.] # mean field force at t
    dt_half      = dt/2.        # half timestep

    ## initialize the initial position and velocities
    @views positions[1, :]  .= copy(pos)
    @views velocities[1, :] .= copy(vel)
    
    ## orbit integration
    for it in 1:nstep
        drift!(pos, vel, dt)                # after one drift: x(t + dt), v(t)
        # wrap_periodic!(pos, fsize, fsize_z) # wrap the position to be periodic before calculating the force

        ## save orbits for every nsave points 
        if mod(it, nsave) == 0
            idx_nsave = div(it, nsave) + 1 

            @views positions[idx_nsave, :] .= pos
            @views velocities[idx_nsave, :] .= vel
        end
    end

end
########################