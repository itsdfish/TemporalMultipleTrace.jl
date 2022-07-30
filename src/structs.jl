abstract type AbstractFMTP end 

"""
    FMTPModel <: AbstractFMTP   

Model object for formal multiple trace theory of temporal preparation.

# Fields

- `τs`: a vector containing the time of at maximum fire rate for each timing cell``
- `κ`: temporal smearing parameter
- `c`: time persistance parameter 
- `λ`: decay parameter
- `u_act_weights`: a dictionary of unique activation weights 
- `u_inhib_weights`: a dictionary of unique inhabition weights

# References

Salet, J. M., Kruijne, W., van Rijn, H., Los, S. A., & Meeter, M. (2022). FMTP: A unifying computational framework of temporal preparation across time scales. 
Psychological Review.
"""
@concrete mutable struct FMTPModel <: AbstractFMTP
    τs
    κ
    c 
    λ
    u_act_weights
    u_inhib_weights    
end

"""
    FMTPModel(;
        τs,
        κ,
        c, 
        λ, 
        u_act_weights,
        u_inhib_weights)

Constructor for the FMPTModel object. 

# Arguments

- `τs`: a vector containing the time of at maximum fire rate for each timing cell``
- `κ`: temporal smearing parameter
- `c`: time persistance parameter 
- `λ`: decay parameter
- `u_act_weights`: a dictionary of unique activation weights 
- `u_inhib_weights`: a dictionary of unique inhabition weights
"""

function FMTPModel(;
    τs,
    κ,
    c, 
    λ, 
    u_act_weights,
    u_inhib_weights)
    return FMTPModel(τs, κ, c, λ, u_act_weights, u_inhib_weights)
end

abstract type AbstractTask end

@concrete mutable struct Task
    dist
    n_trials 
end