module PriorNormalization

# Import packages. 
using KernelDensity # for estimating densities from samples 
using Random 
using SpecialFunctions # for loggamma function 
using StatsBase # for defining customized distributions 
using StatsFuns # for defining customized distributions 
using StatsPlots # for plotting 
using Turing # for setting up the model and sampling 
using Optim # for ML and MAP estimation

# Export specific functionalities you want to be available outside the module

## Generalized gamma distribution 
include("distributions/WikiGeneralizedGamma.jl")
include("distributions/GeneralizedGamma.jl")

## Transforms 
# Proposed prior-normalizing transform 
include("transforms/transportMap_priorNormalizing.jl") 
# Transform proposed by Calvetti and Somersalo in 2024 
include("transforms/transform_CalvettiSomersalo2024.jl")

end # module PriorNormalization