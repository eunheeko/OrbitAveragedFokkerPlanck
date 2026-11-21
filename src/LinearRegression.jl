using StatsBase
using GLM
using DataFrames

#=
linar regression
=#
function lr(df)

    ## linear regression model
    model = lm(@formula(y ~ x), df) 

    ## coefficients and covariance
    mu    = coef(model) # fitted intercept and slope
    sigma = vcov(model) # covariance

    return mu, sigma
end

#=
Bootstrapping function
=#
function bootstrap_lr(tab_time, tab_delta)

    df = DataFrame(x = tab_time, y= tab_delta)
    
    ## linear regression
    mu, sigma = lr(df)

    ## random sample folowing N(mu, sigma)
    d = MvNormal(mu, sigma)
    
    return rand(d, 1)[2] # newly sampled slope
    
    # try 
    #     # linear regression
    #     mu, sigma = lr(df)
    #     # random sample folowing N(mu, sigma)
    #     d = MvNormal(mu, sigma)
    #     return rand(d, 1)[2] # newly sampled slope
    # catch e   
    #    if e isa PosDefException
    #         return NaN   # or NaN
    #     else
    #         rethrow(e)
    #     end
    # end
end