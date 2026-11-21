# module BubModule
    
################################################  
## Default parameters
################################################ 
## Constant 
const DEFAULT_G = 1. # Gravitational constant

## DM halo
const DEFAULT_ρHALO = 1/(4*pi) # Halo mass density
const DEFAULT_RHALO = 1.       # Halo scale radius

## Bubbles
const DEFAULT_RBUB_IN  = 0.50 * DEFAULT_RHALO                                                            # Bubble inner radius [L]
const DEFAULT_RBUB_OUT = 2 * DEFAULT_RBUB_IN                                                             # Bubble outer radius [L]
const DEFAULT_TBUB     = 0.5                                                                             # Bubble lifespan [T]
const DEFAULT_ρBUB_IN  = -0.1 * abs(DEFAULT_ρHALO)                                                       # Bubble inner mass density at r < RBUB_IN [M L-3]
const DEFAULT_ρBUB_OUT = -DEFAULT_ρBUB_IN * DEFAULT_RBUB_IN^3 / (DEFAULT_RBUB_OUT^3 - DEFAULT_RBUB_IN^3) # Bubble outer mass density at RBUB_IN < r < RBUB_OUT [M L-3]

## Gaseous disk
const DEFAULT_dnBUB = 0.2           # Bubble number density per volume per time [L-3 T-1]
const DEFAULT_RDISK = DEFAULT_RHALO # Disk scale radius [L]

## Simulation setting
const DEFAULT_FSIZE    = DEFAULT_RHALO * 4    # Half of the frame size (x, y) [L]
const DEFAULT_FSIZE_Z  = DEFAULT_RHALO * 0.25 # Half of the frame thickness (z) [L]
const DEFAULT_NGRID    = 128                  # number of voxels in (x, y)
const DEFAULT_NGRID_Z  = 8                    # number of voxels in z
const DEFAULT_DEL_TSIM = 0.01                 # timestep [T]
const DEFAULT_TSIM     = 5.                   # total simulation time [T]

## Load modules
include("Param.jl")
include("Halo.jl")
include("Orbit.jl")
include("Bubble.jl")
include("OrbitAveraging.jl")
include("Sampling.jl")
include("IO.jl")
include("LinearRegression.jl")
include("LocalDeflection.jl")
include("LocalDiffusion.jl")

