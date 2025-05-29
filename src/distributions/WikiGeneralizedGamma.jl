# import necessary packages 
import Base.rand, StatsBase.params  # to extend the function "rand"  
import Random, Distributions, Statistics, StatsBase

# export the tools provided in this file 
export WikiGeneralizedGamma 


# Define a new type for the GG distribution (Wikipedia parametrization)
struct WikiGeneralizedGamma{T<:Real} <: ContinuousUnivariateDistribution 
    # Define the parameters 
    a::T # scale parameter 
    d::T # shape parameter 
    p::T # power parameter 

    # use gamma distribution to compute for CDF and quantile method below.  
    gamma_dist::Distributions.Gamma 

    # inner constructor function to instantiate new WikiGeneralizedGamma objects 
    function WikiGeneralizedGamma{T}(a::T, d::T, p::T; check_args = true) where {T<:Real}
        return new{T}(a, d, p, Gamma{T}(d/p,1) )
    end

end


# Constructor functions for implicitly supplied type 
# Constructor for no type and params Float64
function WikiGeneralizedGamma(a::Float64, d::Float64, p::Float64; check_args = true)
    return WikiGeneralizedGamma{Float64}(a, d, p; check_args = check_args)
end


# Constructor for real params - use promote ti make all params the same type 
WikiGeneralizedGamma(a::Real, d::Real, p::Real) = WikiGeneralizedGamma(promote(a, d, p)... ) 
WikiGeneralizedGamma(a::Integer, d::Integer, p::Integer) = WikiGeneralizedGamma(float(a), float(d), float(p) )


## We need to define the following eight methods to set-up our new distribution
# 1) rand(::AbstractRNG, dist::UnivariateDistribution)
# 2) sampler(dist::Distribution)
# 3) logpdf(dist::UnivariateDistribution, x::Real)
# 4) cdf(dist::UnivariateDistribution, x::Real)
# 5) quantile(dist::UnivariateDistribution, q::Real)
# 6) minimum(dist::UnivariateDistribution)
# 7) maximum(dist::UnivariateDistribution)
# 8) insupport(dist::UnivariateDistribution, x::Real)


# 0) Helper function 
StatsBase.params(dist::WikiGeneralizedGamma) = (dist.a, dist.d, dist.p)


# 1) rand(::AbstractRNG, dist::UnivariateDistribution)
# Constructor for the Wiki-GG distribution
function Base.rand(rng::AbstractRNG, dist::WikiGeneralizedGamma) 
    return quantile(dist, rand())
end


# 2) sampler(dist::Distribution)
Distributions.sampler(rng::AbstractRNG, dist::WikiGeneralizedGamma) = Base.rand(rng::AbstractRNG, dist::WikiGeneralizedGamma)


# 3) logpdf(dist::UnivariateDistribution, x::Real)
# Log-PDF for the Wiki-GG distribution 
function Distributions.logpdf(dist::WikiGeneralizedGamma{T}, x::Real) where {T<:Real}
    (a, d, p) = params(dist) # get parameters 
    if x <= 0 
        return zero(T) # equivalent of zero fpr type T 
    else 
        return log(p) - d*log(a) - loggamma(d/p) + (d-1)*log(x) - (x/a)^p

    end
end


# 4) cdf(dist::UnivariateDistribution, x::Real)
# CDF for the Wiki-GG distribution 
function Distributions.cdf(dist::WikiGeneralizedGamma{T}, x::Real) where {T<:Real} 
    (a, d, p) = params(dist) # get parameters 
    return T(Distributions.cdf( dist.gamma_dist, T((x/a)^d) ) )
end 


# 5) quantile(dist::UnivariateDistribution, q::Real)
# Quantile function for the Wiki-GG distribution
function Distributions.quantile(dist::WikiGeneralizedGamma{T}, x::Real) where {T<:Real} 
    (a, d, p) = params(dist) # get parameters 
    return T( a*Distributions.quantile( dist.gamma_dist, T(x))^(1/d) )
end


# 6) minimum(dist::UnivariateDistribution) 
# Define the minimum value of the distribution 
function Base.minimum(dist::WikiGeneralizedGamma) 
    return(0) 
end


# 7) maximum(dist::UnivariateDistribution) 
# Define the maximum value of the distribution  
function Base.maximum(dist::WikiGeneralizedGamma) 
    return(Inf) 
end


# 8) insupport(dist::UnivariateDistribution, x::Real) 
# Check if a value is within the support of the distribution 
function Distributions.insupport(dist::WikiGeneralizedGamma) 
    insupport(dist::WikiGeneralizedGamma, x::Real) = zero(x) < x
end


## Additional recommended statistics functions 

# A) Define the mean of the distribution 
function Statistics.mean(dist::WikiGeneralizedGamma) 
    (a, d, p) = params(dist) # get parameters 
    return a*gamma( (d + 1)/p ) / gamma( d/p )
end 


# B) Define the mode of the distribution 
function StatsBase.mode(dist::WikiGeneralizedGamma{T}) where {T<:Real} 
    (a, d, p) = params(dist) # get parameters 
    if d > 1 
        return T( a*( (d-1)/p )^(1/p) ) 
    else 
        return T(0)
    end
end 


# C) Define the variance of the distribution 
function Statistics.var(dist::WikiGeneralizedGamma) 
    (a, d, p) = params(dist) # get parameters 
    return a^2 * ( gamma((d+2)/p) / gamma(d/p) - ( gamma((d+1)/p) / gamma(d/p) )^2 )
end 