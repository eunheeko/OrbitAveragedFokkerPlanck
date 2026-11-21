using NaNMath

include("Param.jl")

#=
 NFW potential
=#
function Ψ_NFW(r::Float64, G::Float64, rhalo::Float64, ρhalo::Float64, ϵ::Float64)::Float64
    
    ## softened r to avoid divergence near r = 0
    rsoft = sqrt(r^2 + ϵ^2)
    
    ## constant
    A = 4 * π * G * ρhalo * rhalo^3 
    
    return -A / rsoft * NaNMath.log(1 + rsoft / rhalo)
end

#=
1st derivative of NFW potential (consistent with softend potential)
=#
function dΨ_NFW(r::Float64, G::Float64, rhalo::Float64, ρhalo::Float64, ϵ::Float64)::Float64
    
    ## softened r to avoid divergence near r = 0
    rsoft = sqrt(r^2 + ϵ^2)

    ## contants
    A = 4 * π * G * ρhalo * rhalo^3 
    B = NaNMath.log(1 + rsoft/rhalo) - rsoft / (rhalo + rsoft)
    
    return A * r / rsoft^3 * B
end

#=
2nd derivative of NFW potential (consistent with softend potential)
=#
function d2Ψ_NFW(r::Float64, G::Float64, rhalo::Float64, ρhalo::Float64, ϵ::Float64)::Float64
    
    ## softened r to avoid divergence near r = 0
    rsoft = sqrt(r^2 + ϵ^2)

    ## contants
    A = 4 * π * G * ρhalo * rhalo^3 
    B = NaNMath.log(1 + rsoft/rhalo) - rsoft / (rhalo + rsoft)
    C = r^2 * rsoft^2 / (rhalo + rsoft)^2
    
    return A * (
        B * (rsoft^2 - 3 * r^2) - C
    ) / rsoft^5
end

#=
 NFW potential in cartesian position vector
=#
function Ψ_NFW_cartesian(pos::Vector, G::Float64, rhalo::Float64, ρhalo::Float64, ϵ::Float64)
    r = hypot(pos...)
    return Ψ_NFW(r, G, rhalo, ρhalo, ϵ)
end

#=
1st derivative of NFW potential in cartesian position vector (consistent with softend potential)
=#
function dΨ_NFW_cartesian(pos::Vector, G::Float64, rhalo::Float64, ρhalo::Float64, ϵ::Float64)
    
    x, y, z = pos
    r = hypot(pos...)
    
    dphi = dΨ_NFW(r, G, rhalo, ρhalo, ϵ)
    
    return [dphi * x/r, dphi * y/r, dphi * z/r]
end

## 1st derivative of NFW potential
function dΨ_NFW_IC(r::Float64, G::Float64, rhalo::Float64, ρhalo::Float64, ϵ::Float64)::Float64
    
    ## softened r to avoid divergence near r = 0
    rsoft = sqrt(r^2 + ϵ^2)
    
    return 4 * π * G * ρhalo * rhalo^3 / rsoft^2 * NaNMath.log(1+rsoft/rhalo) - 4 * π * G * ρhalo * rhalo^2 / rsoft / (1 + rsoft/rhalo) 
end

## 2nd derivative of NFW potential
function d2Ψ_NFW_IC(r::Float64, G::Float64, rhalo::Float64, ρhalo::Float64, ϵ::Float64)::Float64
    
    ## softened r to avoid divergence near r = 0
    rsoft = sqrt(r^2 + ϵ^2)
    
    return -8 * π * G * ρhalo * rhalo^3 / rsoft^3 * NaNMath.log(1+rsoft/rhalo) + 8 * π * G * ρhalo * rhalo^2 / rsoft^2 / (1 + rsoft/rhalo) + 4 * π * G * ρhalo * rhalo /  rsoft / (1 + rsoft/rhalo)^2
end

#=
# 1st derivative of NFW potential in cartesian position vector (consistent with softend potential)
# =#
# function dΨ_NFWp_cartesian(pos::Vector, G::Float64, rhalo::Float64, ρhalo::Float64, ϵ::Float64, FSIZE::Float64, FSIZE_Z::Float64)
    
#     xp, yp, zp = pos
#     rp = hypot(pos...)

#     r = sqrt((xp - FSIZE)^2 + (yp - FSIZE)^2 + (zp - FSIZE_Z)^2)
    
#     dphi = dΨ_NFW(r, G, rhalo, ρhalo, ϵ)
    
#     return [dphi * xp/rp, dphi * yp/rp, dphi * zp/rp]
# end