using LinearAlgebra

#=
cross product
=#
function cal_cross_product(vec_in1::AbstractMatrix, vec_in2::AbstractMatrix)
    
    ## initialize a table to save the final outer product
    vec_out = zeros(axes(vec_in1))

    @inbounds for iax in axes(vec_out, 1)
        vec_out[iax, :] .= cross(vec_in1[iax, :], vec_in2[iax, :])
    end

    return vec_out
end

#=
dot product
=#
function cal_dot_product(vec_in1::AbstractMatrix, vec_in2::AbstractMatrix)
    
    ## initialize a table to save the final dot product
    vec_out = zeros(axes(vec_in1, 1))
    
    @inbounds for iax in axes(vec_out, 1)
        vec_out[iax] = dot(vec_in1[iax, :], vec_in2[iax, :])
    end

    return vec_out
end

#=
calculate norm
=#
function cal_norm(vec_in::AbstractMatrix; ϵnorm::Float64=1e-15)

    ## initialize a table to save the final norm
    norm_out = zeros(axes(vec_in, 1))

    @inbounds for iax in axes(vec_in, 1)
        norm_in = norm(vec_in[iax, :])
        # n_in    = max(norm_in, ϵnorm) # in case of small norm
        n_in = norm_in 
        
        norm_out[iax] = n_in
    end

    return norm_out
end

#=
calculate unit vector
=#
function cal_unit_vector(vec_in::AbstractMatrix; ϵnorm::Float64=1e-15)

    ## initialize a table to save the final unit vector
    vec_unit_out = zeros(axes(vec_in))
    norm_in      = cal_norm(vec_in; ϵnorm = ϵnorm)

    @inbounds for iax in axes(vec_in, 1)
        n_in = norm_in[iax]
        vec_unit_out[iax, :] .= vec_in[iax, :] ./ n_in
    end

    return vec_unit_out
end

#=
calculate force fluctuations at a given position at a given time
=#
function cal_δF(vec_r::AbstractMatrix, cal_δF_local!::Function, times_save::AbstractVector)

    ## pre-allocation for force fluctuations 
    vec_δF = zeros(axes(vec_r)) ## force fluctuations
    
    ## pre-allocation of vectors for intermediate calculations
    _vec_pos_loc = [0., 0., 0.] # position vector 𝐫 at t
    _vec_δF_loc  = [0., 0., 0.] # forces vector 𝐅 at t

    ## calculate force fluctuations at t
    @inbounds for iax in axes(vec_r, 1)
        
        ## load time t
        t = times_save[iax]

        ## update position and velocity of a test particle
        _vec_pos_loc .= vec_r[iax, :]

        ## update force_bubble at (r, t)
        cal_δF_local!(_vec_pos_loc, _vec_δF_loc, t)

        #@ save SB force at t [δF]
        vec_δF[iax, :] .= _vec_δF_loc
    end

    return vec_δF

end


#=
calculate local deflection in (E, L)
=#
function local_orbit_change(vec_r::AbstractMatrix, vec_v::AbstractMatrix, vec_δv::AbstractMatrix, type_orbit::String)

    ## local changes in E
    δE      = cal_dot_product(vec_v, vec_δv) .+ 0.5 .* cal_dot_product(vec_δv, vec_δv) # δE = 𝐯 · δ𝐯 + 0.5 |δ𝐯|² + 𝑂(δ³)
    δEδE    = cal_dot_product(vec_v, vec_δv).^2                                        # (δE)² = (𝐯 · δ𝐯)² + 𝑂(δ³)
    
    ## local changes in L
    if type_orbit == "1D-radial" # radial orbit: no L direction

        vec_δL  = cal_cross_product(vec_r, vec_δv)  # angular momentum deflection vector δ𝐋 = 𝐫 × δ𝐯
        
        δL      = cal_norm(vec_δL)                  # |δ𝐋| + 𝑂(δ³)
        δLδL    = cal_dot_product(vec_δL, vec_δL)   # (δL)² = |δ𝐋|² + 𝑂(δ³)
        δEδL    = δE .* δL

        return δE, δEδE, δL, δLδL, δEδL
    else # 2D orbit

        vec_L   = cal_cross_product(vec_r, vec_v)  # angular momentum vector 𝐋 = 𝐫 × 𝐯
        vec_δL  = cal_cross_product(vec_r, vec_δv) # angular momentum deflection vector δ𝐋 = 𝐫 × δ𝐯
        norm_L  = cal_norm(vec_L)                  # |L|: norm of angular momentum vector

        if norm_L .>= 1e-10 # non-radial orbit
            norm_LδL  = cal_dot_product(vec_L, vec_δL)
            norm_δLδL = cal_dot_product(vec_δL, vec_δL)

            δL    = norm_LδL ./ norm_L .+ 0.5 .* norm_δLδL.^2 ./ norm_L .- 0.5 .* (norm_LδL.^2) ./ (norm_L .^3)
            δLδL  = (norm_LδL ./ norm_L).^2
            δEδL  = δE .* δL

            return δE, δEδE, δL, δLδL, δEδL
        else # radial orbit limit
            
            δL      = cal_norm(vec_δL)                  # |δ𝐋| + 𝑂(δ³)
            δLδL    = cal_dot_product(vec_δL, vec_δL)   # (δL)² = |δ𝐋|² + 𝑂(δ³)
            δEδL    = δE .* δL

        return δE, δEδE, δL, δLδL, δEδL
        end

    end

end

# function from_ref_to_loc(vec_ref::AbstractMatrix, vec_basis_loc1::AbstractMatrix, vec_basis_loc2::AbstractMatrix, vec_basis_loc3::AbstractMatrix)
    
#     vec_loc = zeros(axes(vec_ref))

#     @inbounds for iax in axes(vec_loc, 1)

#         mat_R = hcat(vec_basis_loc1[iax, :], vec_basis_loc2[iax, :], vec_basis_loc3[iax, :])' # z = L
#         # mat_R = hcat(vec_basis_loc1[iax, :], vec_basis_loc2[iax, :], vec_basis_loc3[iax, :])

#         vec_loc[iax, :] .= mat_R * vec_ref[iax, :]
#         # vec_loc[iax, 1] = dot(vec_ref[iax, :], vec_basis_loc1[iax, :])
#         # vec_loc[iax, 2] = dot(vec_ref[iax, :], vec_basis_loc2[iax, :])
#         # vec_loc[iax, 3] = dot(vec_ref[iax, :], vec_basis_loc3[iax, :])
#     end

#     return vec_loc
# end


# ## local comoving frame basis (v', L', w')
# function cal_local_basis(vec_r::AbstractMatrix, vec_v::AbstractMatrix, L₀::Float64; ϵnorm  = ϵnorm)

#     ## v' ∥ v
#     vec_vp_hat = cal_unit_vector(vec_v; ϵnorm  = ϵnorm)

#     ## L' ∥ z
#     vec_L = zero(vec_r)

#     # radial orbit
#     if isnan(L₀)
#         vec_L[:, 3] .= 1.0 
#     else # 2D orbit
#         vec_L .= cal_cross_product(vec_r, vec_v)
#     end
#     vec_Lp_hat = cal_unit_vector(vec_L)
    
#     ## w' = (v × L)/|v×L|
#     vec_wp = cal_cross_product(vec_vp_hat, vec_Lp_hat)
#     vec_wp_hat = cal_unit_vector(vec_wp; ϵnorm  = ϵnorm)

#     return vec_vp_hat, vec_Lp_hat, vec_wp_hat

# end

# function cal_local_basis(vec_r::AbstractMatrix, vec_v::AbstractMatrix; ϵnorm  = ϵnorm)

#     ## L vector
#     vec_L = zeros(vec_r)

#     # radial orbit
#     if cal_dot_product(vec_r, vec_v) <= ϵnorm
#         vec_L[:, 3] .= 1.0 
#     else
#         vec_L .= cross(vec_r, vec_v)
#     end

#     ##  build local frame (r, t, L)
#     vec_rp_hat = cal_unit_vector(vec_r; ϵnorm  = ϵnorm)

#     vec_Lp_hat = cal_unit_vector(vec_L; ϵnorm  = ϵnorm)

#     vec_tp = cal_cross_product(vec_L_hat, vec_r_hat)
#     vec_tp_hat = cal_unit_vector(vec_tp; ϵnorm  = ϵnorm)

#     return vec_rp_hat, vec_tp_hat, vec_Lp_hat
    	
# end

# function local_change(vec_r::AbstractMatrix, vec_v::AbstractMatrix, vec_δv::AbstractMatrix, L₀::Float64)

#     vec_vp_hat, vec_Lp_hat, vec_wp_hat = cal_local_basis(vec_r, vec_v, L₀)
#     vec_L = cal_cross_product(vec_r, vec_v)
#     vec_δL = cal_cross_product(vec_r, vec_δv)

#     # velocity vector in local comoving frame (v', L', w')
#     # vec_r_loc = from_ref_to_loc(vec_r, vec_vp_hat, vec_Lp_hat, vec_wp_hat)
#     vec_v_loc = from_ref_to_loc(vec_v, vec_vp_hat, vec_Lp_hat, vec_wp_hat)
#     vec_δv_loc = from_ref_to_loc(vec_δv, vec_vp_hat, vec_Lp_hat, vec_wp_hat)
#     vec_L_loc = from_ref_to_loc(vec_L, vec_vp_hat, vec_Lp_hat, vec_wp_hat)
#     vec_δL_loc = from_ref_to_loc(vec_δL, vec_vp_hat, vec_Lp_hat, vec_wp_hat)

#     δE = cal_δE(vec_v_loc, vec_δv_loc)
#     δEδE = cal_δEδE(vec_v_loc, vec_δv_loc)

#     if isnan(L₀) # radial orbit
#         δL = fill(NaN, size(vec_r, 1))
#         δLδL = cal_dot_product(vec_δL_loc, vec_δL_loc)
#         δEδL = fill(NaN, size(vec_r, 1))
#         return δE, δEδE, δL, δLδL, δEδL
#     else # 2D orbit
#         norm_L_loc = cal_norm(vec_L_loc) # |L|: norm of angular momentum vector
#         L_dot_δL_loc = cal_dot_product(vec_δL_loc, vec_L_loc) # angular momentum deflection δL = 𝐋̂ · δ𝐋
#         δL_dot_δL_loc = cal_dot_product(vec_δL_loc, vec_δL_loc) # angular momentum deflection δL = 𝐋̂ · δ𝐋
        
#         δL = L_dot_δL_loc ./ norm_L_loc .+ 0.5 .* (δL_dot_δL_loc .^2) ./ norm_L_loc .- 0.5 .* (L_dot_δL_loc .^2) ./ (norm_L_loc .^3)
#         δLδL = (L_dot_δL_loc ./ norm_L_loc).^2
#         δEδL = δE .* δL

#         return δE, δEδE, δL, δLδL, δEδL
#     end

# end

