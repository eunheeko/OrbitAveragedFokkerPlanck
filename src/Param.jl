using JSON
using Logging

# export BubbleParameters
"""
    Parameters for bubble generation
"""
struct BubbleParameters

    ## constant 
    G::Float64 # Gravitational constant

    ## DM halo
    ρHALO::Float64 # Halo mass density
    RHALO::Float64 # Halo scale radius

    ## Bubbles
    RBUB_IN::Float64  # Bubble inner radius [L]
    RBUB_OUT::Float64 # Bubble outer radius [L]
    TBUB::Float64     # Bubble lifespan [T]
    ρBUB_IN::Float64  # Bubble inner mass density at r < RBUB_IN [M L-3]
    ρBUB_OUT::Float64 # Bubble outer mass density at RBUB_IN < r < RBUB_OUT [M L-3]

    ## Gaseous disk
    dnBUB::Float64 # Bubble number density per volume per time [L-3 T-1]
    RDISK::Float64 # Disk scale radius [L]

    ## Simulation setting
    FSIZE::Float64    # Half of the frame size (x, y) [L]
    FSIZE_Z::Float64  # Halof of the  frame thickness (z) [L]
    NGRID::Int        # number of voxels in (x, y)
    NGRID_Z::Int      # number of voxels in z
    DEL_TSIM::Float64 # timestep [T]
    TSIM::Float64     # total simulation time [T]
end

function BubbleParameters(
    ;G::Float64       = DEFAULT_G,
    ρHALO::Float64    = DEFAULT_ρHALO,
    RHALO::Float64    = DEFAULT_RHALO,
    RBUB_IN::Float64  = DEFAULT_RBUB_IN,
    RBUB_OUT::Float64 = DEFAULT_RBUB_OUT,
    TBUB::Float64     = DEFAULT_TBUB,
    ρBUB_IN::Float64  = DEFAULT_ρBUB_IN,
    ρBUB_OUT::Float64 = DEFAULT_ρBUB_OUT,
    dnBUB::Float64    = DEFAULT_dnBUB,
    RDISK::Float64    = DEFAULT_RDISK,
    FSIZE::Float64    = DEFAULT_FSIZE,
    FSIZE_Z::Float64  = DEFAULT_FSIZE_Z,
    NGRID::Int        = DEFAULT_NGRID,
    NGRID_Z::Int      = DEFAULT_NGRID_Z,
    DEL_TSIM::Float64 = DEFAULT_DEL_TSIM,
    TSIM::Float64     = DEFAULT_TSIM,
    )
    return BubbleParameters(G, ρHALO, RHALO, RBUB_IN, RBUB_OUT, TBUB, ρBUB_IN, ρBUB_OUT, dnBUB, RDISK, FSIZE, FSIZE_Z, NGRID, NGRID_Z, DEL_TSIM, TSIM)
end

"""
    Parameters for Leapfrog integrator
"""
### Leapfrog
struct LeapfrogParameter
    dt::Float64
    tsim::Float64
    times::Vector
    nstep::Int
    nsave::Int
    ϵ::Float64
    ρfactor::Float64
    fsize::Float64
    fsize_z::Float64
end

function LeapfrogParameter(
   ;dt::Float64,
   tsim::Float64,
   times::Vector,
   nstep::Int,
   nsave::Int,
   ϵ::Float64,
   ρfactor::Float64,
   fsize::Float64,
   fsize_z::Float64
    )
    return LeapfrogParameter(dt, tsim, times, nstep, nsave, ϵ, ρfactor, fsize, fsize_z)
end