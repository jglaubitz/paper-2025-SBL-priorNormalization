# import necessary packages 
import Base.rand, StatsBase.params  # to extend the function "rand"  

# export the tools provided in this file 
export GeneralizedGamma 


# Define a new type for the Generalized Gamma distribution
struct GeneralizedGamma{T<:Real} <: ContinuousUnivariateDistribution
    # Define the parameters 
    r::T # power parameter 
    β::T # shape parameter
    ϑ::T # scale parameter
    # rβ::T 

    # use gamma distribution to compute for CDF and quantile method below.  
    gamma_dist::Distributions.Gamma 

    function GeneralizedGamma{T}(r::T, β::T, ϑ::T; check_args=true) where {T<:Real}
        # Check if appropriate parameters are used 
        if check_args && r==0
            throw("Parameter r must be non-zero.")
        elseif check_args && (β<=0 || ϑ<=0)
            throw("Parameters β and ϑ must be positive.")       
        end
        new{T}(r, β, ϑ, Gamma{T}(β,1) )
    end
end


# Constructor functions for implicitly supplied type 
# Constructor for no type and params Float64
function GeneralizedGamma(r::Float64, β::Float64, ϑ::Float64; check_args = true)
    return GeneralizedGamma{Float64}(r, β, ϑ; check_args = check_args)
end


# Constructor for real params - use promote to make all params the same type 
GeneralizedGamma(r::Real, β::Real, ϑ::Real) = GeneralizedGamma(promote(r, β, ϑ)... ) 
GeneralizedGamma(r::Integer, β::Integer, ϑ::Integer) = GeneralizedGamma(float(r), float(β), float(ϑ) )

# Implement the conversion logic from GeneralizedGamma to WikiGeneralizedGamma parameters
function show(io::IO, dist::GeneralizedGamma)
    print(io, "GeneralizedGamma (r=$(dist.r), β=$(dist.β), ϑ=$(dist.ϑ))")
end

Distributions.params(dist::GeneralizedGamma{T}) where {T <: Real} = (dist.r, dist.β, dist.ϑ)

function Distributions.logpdf(dist::GeneralizedGamma{T}, x::Real) where {T<:Real}
    (r, β, ϑ) = params(dist) # get parameters 
    if x <= 0 
        return zero(T) # equivalent of zero for type T 
    else 
        return log(abs(r)) - r*β*log(ϑ) - loggamma(β) + (r*β-1)*log(x) - (x/ϑ)^r
    end
end

# # Utilize the WikiGeneralizedGamma under the hood
# function wikiGG_from_generalizedGG(dist::GeneralizedGamma) 
#     a = dist.ϑ
#     d = dist.r * dist.β
#     p = dist.r
#     WikiGeneralizedGamma(a, d, p)
# end

# Define required methods by forwarding to the GeneralizedGamma equivalent
## We need to define the following eight methods to set-up our new distribution
# 1) rand(::AbstractRNG, dist::UnivariateDistribution)
# 2) sampler(dist::Distribution)
# 3) logpdf(dist::UnivariateDistribution, x::Real)
# 4) cdf(dist::UnivariateDistribution, x::Real)
# 5) quantile(dist::UnivariateDistribution, q::Real)
# 6) minimum(dist::UnivariateDistribution)
# 7) maximum(dist::UnivariateDistribution)
# 8) insupport(dist::UnivariateDistribution, x::Real)



# 4) cdf(dist::UnivariateDistribution, x::Real)
# CDF for the GG distribution 
#function Distributions.cdf(dist::GeneralizedGamma{T}, x::Real) where {T<:Real} 
#    (r, β, ϑ) = params(dist) # get parameters 
#    return T(Distributions.cdf( dist.gamma_dist, T((x/ϑ)^(r*β)) ) )
#end 
function Distributions.cdf(dist::GeneralizedGamma{T}, x::Real) where {T<:Real} 
    (r, β, ϑ) = params(dist) # get parameters 
    if r > 0 
        return T(Distributions.cdf(dist.gamma_dist, T((x/ϑ)^(r*β)))) 
    else 
        return 1 - T(Distributions.cdf(dist.gamma_dist, T((x/ϑ)^(r*β))))
    end
end


# 5) quantile(dist::UnivariateDistribution, q::Real)
# Quantile function for the Wiki-GG distribution
#function Distributions.quantile(dist::GeneralizedGamma{T}, x::Real) where {T<:Real} 
#    (r, β, ϑ) = params(dist) # get parameters 
#    return T( ϑ*Distributions.quantile( dist.gamma_dist, T(x))^(1/(r*β)) )
#end
function Distributions.quantile(dist::GeneralizedGamma{T}, x::Real) where {T<:Real} 
    (r, β, ϑ) = params(dist) # get parameters 
    if r > 0 
        return T( ϑ*Distributions.quantile( dist.gamma_dist, T(x))^(1/(r*β)) )
    else 
        return T( ϑ*Distributions.quantile( dist.gamma_dist, T(1-x))^(1/(r*β)) )
    end
end


# 6) minimum(dist::UnivariateDistribution) 
# Define the minimum value of the distribution 
function Base.minimum(dist::GeneralizedGamma) 
    return(0) 
end


# 7) maximum(dist::UnivariateDistribution) 
# Define the maximum value of the distribution  
function Base.maximum(dist::GeneralizedGamma) 
    return(Inf) 
end


# 8) insupport(dist::UnivariateDistribution, x::Real) 
# Check if a value is within the support of the distribution 
function Distributions.insupport(dist::GeneralizedGamma) 
    insupport(dist::GeneralizedGamma, x::Real) = zero(x) < x
end


## Additional recommended statistics functions 

# A) Define the mean of the distribution 
function Statistics.mean(dist::GeneralizedGamma) 
    (r, β, ϑ) = params(dist) # get parameters 
    return ϑ*gamma( (r*β + 1)/r ) / gamma( β )
end 


# B) Define the mode of the distribution 
function StatsBase.mode(dist::GeneralizedGamma{T}) where {T<:Real} 
   (r, β, ϑ) = params(dist) # get parameters 
    if r*β > 1 
        return T( ϑ*( (r*β-1)/r )^(1/r) ) 
    else 
        return T(0)
    end
end 


# C) Define the variance of the distribution 
function Statistics.var(dist::GeneralizedGamma) 
    (r, β, ϑ) = params(dist) # get parameters 
    return ϑ^2 * ( gamma((r*β+2)/r) / gamma(β) - ( gamma((r*β+1)/r) / gamma(β) )^2 )
end 



# Utilize the WikiGeneralizedGamma under the hood
function wikiGG_from_generalizedGG(dist::GeneralizedGamma) 
    a = dist.ϑ
    d = dist.r * dist.β
    p = dist.r
    WikiGeneralizedGamma(a, d, p)
end


# # 1) rand(::AbstractRNG, dist::UnivariateDistribution)
# Base.rand(rng::AbstractRNG, dist::GeneralizedGamma) = rand(rng, wikiGG_from_generalizedGG(dist))

# # 2) sampler(dist::Distribution) 
# Distributions.sampler(rng::AbstractRNG, dist::GeneralizedGamma) = Base.rand(rng::AbstractRNG, dist::GeneralizedGamma)

# # 3) logpdf(dist::UnivariateDistribution, x::Real) 
# Distributions.logpdf(dist::GeneralizedGamma, x::Real) = logpdf(wikiGG_from_generalizedGG(dist), x)

# # 4) cdf(dist::UnivariateDistribution, x::Real) 
# Distributions.cdf(dist::GeneralizedGamma, x::Real) = cdf(wikiGG_from_generalizedGG(dist), x)

# # 5) quantile(dist::UnivariateDistribution, q::Real) 
# Distributions.quantile(dist::GeneralizedGamma, q::Real) = quantile(wikiGG_from_generalizedGG(dist), q)

# # 6) minimum(dist::UnivariateDistribution) 
# Base.minimum(dist::GeneralizedGamma) = minimum(wikiGG_from_generalizedGG(dist))

# # 7) maximum(dist::UnivariateDistribution) 
# Base.maximum(dist::GeneralizedGamma) = maximum(wikiGG_from_generalizedGG(dist))

# # 8) insupport(dist::UnivariateDistribution, x::Real) 
# Distributions.insupport(dist::GeneralizedGamma, x::Real) = insupport(wikiGG_from_generalizedGG(dist), x)


# ## Additional recommended statistics functions 

# # A) Define the mean of the distribution 
# Statistics.mean(dist::GeneralizedGamma) = mean(wikiGG_from_generalizedGG(dist))

# # B) Define the mode of the distribution 
# StatsBase.mode(dist::GeneralizedGamma) = mode(wikiGG_from_generalizedGG(dist))

# # C) Define the variance of the distribution
# Statistics.var(dist::GeneralizedGamma) = var(wikiGG_from_generalizedGG(dist))


