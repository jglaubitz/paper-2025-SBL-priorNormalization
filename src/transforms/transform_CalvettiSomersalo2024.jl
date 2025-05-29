# import necessary packages 
#import ApproxFun

# Transform proposed by Calvetti and Somersalo in their 2024 paper 
# "Computationally efficient sampling methods for sparsity promoting hierarchical Bayesian models"

# export the tools provided in this file 
export transform_CalvettiSomersalo2024 

# Define the transform T = [t_1,...,t_n] with t_i: R² → R₊×R, [u,τ] ↦ [x,θ] 
function transform_CalvettiSomersalo2024( v, ω; r, β, ϑ )  
    # Compute θ from v and ω
    θ = 2^(-1/r) * ϑ .* abs.(ω).^(2/r)
    
    # Compute x from v and ω
    x = v .* sqrt.(θ)
    
    # Return the results 
    return x, θ
end