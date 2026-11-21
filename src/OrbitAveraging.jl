# include("Param.jl")
import OrbitalElements as OE

#=============================================
Numerical orbit-averaging
=============================================#

#=
PDF in the anomaly u (Eq. A6 in Petersen+24)
=#
function pu(u::Float64, sma::Float64, ecc::Float64, Ωr::Float64, Ψ::Function, dΨ::Function, d2Ψ::Function, params::OE.OrbitalParameters) 
    return Ωr * OE.ΘAE(Ψ, dΨ, d2Ψ, u, sma, ecc, params) / pi # Ωr/π 1/vᵣ dr/du 
end

#=
r of u
=#
r_of_u(u::Float64, sma::Float64, ecc::Float64) = OE.ru(u,sma,ecc)


#=
 v of r (total velocity |v|)
=#
function v_of_ru(u::Float64, sma::Float64, ecc::Float64, E::Float64, Ψ::Function)
    r = r_of_u(u, sma, ecc)
    return sqrt(2*(E - Ψ(r)))
end


#=
function to compute the orbit-averaged diffusion coefficients
=#
function orbit_averaging(local_dc::Vector, Σ::Vector, vec_u::Vector, sma::Float64, ecc::Float64, Ωr::Float64 ,Ψ::Function, dΨ::Function, d2Ψ::Function, params::OE.OrbitalParameters)
    
    ## trapezoidal rule
    int_du = pu.(vec_u, sma, ecc, Ωr, Ψ, dΨ, d2Ψ, Ref(params)) .* local_dc .* Σ # Ωr/π 1/vᵣ dr/du × d × Σ(r[u])
    
    ## orbit-averaged diffusion coefficient
    global_dc = sum(
        0.5 * (int_du[1:end-1] .+ int_du[2:end]) .* diff(vec_u)
    )

    return global_dc 
end