# import necessary packages 
#import ApproxFun

# export the tools provided in this file 
# prior-normalizing KR map s: R₊×R → R², [x,θ] ↦ [u,τ] 
export priorNormalizing_KR 
export priorNormalizing_KR_sθ
export priorNormalizing_KR_sx
# inverse of the prior-normalizing KR map, T: R² → R₊×R, [u,τ] ↦ [x,θ]
export priorNormalizing_KR_inv
export priorNormalizing_KR_inv_tτ 
export priorNormalizing_KR_inv_tu
export priorNormalizing_KR_inv_tu_fast


# Define the inverse of the prior-normalizing KR map, t: R² → R₊×R, [u,τ] ↦ [x,θ] 
# The last argument, Φ, is the function (approximation) for gammainvccdf( β, 1, 0.5*erfc( τ/sqrt(2) ) ) 
function priorNormalizing_KR_inv( u, τ; r, β, ϑ, Φ::Function )  
    θ = priorNormalizing_KR_inv_tτ( τ; r, β, ϑ, Φ ) 
    x = priorNormalizing_KR_inv_tu( u, τ; r, β, ϑ, Φ )
    return x, θ
end


# Define the first component of the inverse KR map, t^τ
function priorNormalizing_KR_inv_tτ_aux( τ; r, β, ϑ, Φ::Function ) 
    # Compute quantile of gamma distribution
    #gamma_quantile = gammainvccdf( β, 1, 0.5*erfc( τ/sqrt(2) ) ) 
    gamma_quantile = abs.(Φ.(τ)) 

    # Translate from gamma to generalized gamma distribution 
    θ = ϑ * gamma_quantile.^( 1/r ) # transform into GG(r,ϑ,β) 

    return θ
end


# Define the first component of the inverse KR map, t^τ
function priorNormalizing_KR_inv_tτ( τ; r, β, ϑ, Φ::Function ) 
    B = 5
    if abs(τ) <= B 
        θ = priorNormalizing_KR_inv_tτ_aux( τ; r, β, ϑ, Φ )
    elseif τ < -B 
        # Value of t^τ at τ = R 
        θ_B = priorNormalizing_KR_inv_tτ_aux( -B; r, β, ϑ, Φ ) 
        # Value of t^τ at τ = R-δ 
        δ = .1*B # scaled FD step size
        θ_δ = priorNormalizing_KR_inv_tτ_aux( δ-B; r, β, ϑ, Φ )
        # Compute (t^τ)'(R) using an FD approximation 
        θ_deriv = ( θ_δ - θ_B ) / δ
        # compute the value of the linear extension 
        θ = θ_B + θ_deriv.*( τ + B )
    elseif τ > B 
        # Value of t^τ at τ = R 
        θ_B = priorNormalizing_KR_inv_tτ_aux( B; r, β, ϑ, Φ ) 
        # Value of t^τ at τ = R-δ 
        δ = .1*B # scaled FD step size
        θ_δ = priorNormalizing_KR_inv_tτ_aux( B-δ; r, β, ϑ, Φ )
        # Compute (t^τ)'(R) using an FD approximation 
        θ_deriv = ( θ_B - θ_δ ) / δ
        # compute the value of the linear extension 
        θ = θ_B + θ_deriv.*( τ - B )
    end

    # Make sure that output is positive 
    if θ > 0 
        return θ 
    else 
        return 0 
    end
end


# Define the second component of the inverse KR map, t^{u|τ}
function priorNormalizing_KR_inv_tu( u, τ; r, β, ϑ, Φ::Function ) 
    θ = priorNormalizing_KR_inv_tτ( τ; r, β, ϑ, Φ )
    x = u .* sqrt.(θ)
    return x
end


# Fast version of the inverse KR map, t: R² → R₊×R, [u,τ] ↦ [x,θ] 
# This version assumes that τ is ensured to be in the domain [-5,5]
# Define the second component of the inverse KR map, t^{u|τ}
function priorNormalizing_KR_inv_tu_fast( u, τ; r, β, ϑ, Φ::Function ) 
    θ = priorNormalizing_KR_inv_tτ_aux( τ; r, β, ϑ, Φ )
    x = u .* sqrt.(θ)
    return x
end



# Define the prior-normalizing KR map s: R₊×R → R², [x,θ] ↦ [u,τ] 
function priorNormalizing_KR( x, θ; r, β, ϑ )  
    u = priorNormalizing_KR_sx( x, θ; r, β, ϑ )
    τ = priorNormalizing_KR_sθ( θ; r, β, ϑ ) 
    return u, τ
end


# Define the first component of the KR map, s^θ
function priorNormalizing_KR_sθ( θ; r, β, ϑ ) 
    distr_ref = Normal() # set-up the reference standard normal distribution
    distr_target = Gamma( β, 1 ) # set-up the target gamma distribution

    if θ <= 0        
        τ = 0 
    else
        GG_CDF_aux = cdf( distr_target, ( θ/ϑ ).^r ) # compute the CDF of GG(r,β,ϑ) 
        if r > 0 
            GG_CDF = GG_CDF_aux
        else 
            GG_CDF = 1 - GG_CDF_aux
        end
        τ = quantile( distr_ref, GG_CDF ) # compute the quantile of the standard normal distribution    
    end

    return τ
end


# Define the second component of the KR map, s^{x|θ}
function priorNormalizing_KR_sx( x, θ; r, β, ϑ ) 
    u = x/sqrt(θ)
end