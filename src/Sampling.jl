import OrbitalElements as OE

#=
sample azimuthal angle φ ∈ [0, 2π]
=#
function mcsample_φ()
    return 2*pi*rand()
end

#=
sample polar angle θ ∈ [0, π]
=#
function mcsample_θ()
    # return pi*rand()
    return acos(2*rand() - 1)
end

#=
sample radius
=#
function mcsample_r(u::Float64, sma::Float64, ecc::Float64)
    return  OE.ru(u,sma,ecc)
end

#=
wrapper function to sample initial condition on the orbital plane.
=#
function mcsample_2d(E::Float64, L::Float64, rp::Float64, ra::Float64, sma::Float64, ecc::Float64,Ψ::F0,dΨ::F1,d2Ψ::F2,params::OE.OrbitalParameters) where {F0 <: Function, F1 <: Function, F2 <: Function}

    ## load frequencies
    Ωr, _ = OE.ComputeFrequenciesAE(Ψ,dΨ,d2Ψ,sma,ecc,params)

    ## PDF in the anomaly u
    function pu(u::Float64) 
        return Ωr * OE.ΘAE(Ψ,dΨ,d2Ψ,u,sma,ecc,params) / pi
    end
    maxpu = max(pu(1.0),pu(-1.0))

    ## sample azimuthal angle
    φ = mcsample_φ()

    ## sample r and velocities
    if ra == rp # circular orbit
        r = ra
        vr = 0.
        vφ = L/r^2
        return r, φ, vr, vφ
    else
        ## pre-allocation
        u = 0.
        while true
            # uniform between -1 and 1
            u = 2*rand()-1

            # Acceptance test
            utest = rand()
            if utest < pu(u) / maxpu # Reject condition with uniform control PDF
                break
            end
        end

        ## sample r
        r = mcsample_r(u, sma, ecc)

        ## sample random sign
        sign_vr = rand((-1.0, 1.0))
        sign_vφ = rand((-1.0, 1.0))
        
        ## sample velocity
        vr = sign_vr * OE.Vrad(Ψ,dΨ,u,sma,ecc,params)
        vφ = sign_vφ * L/r^2
        
        return r, φ, vr, vφ
    end
end