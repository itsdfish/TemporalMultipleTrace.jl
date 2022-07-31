abstract type AbstractFMTP end 

"""
    FMTPModel <: AbstractFMTP   

Model object for formal multiple trace theory of temporal preparation.

# Fields

- `τs`: a vector containing the time of at maximum fire rate for each timing cell``
- `κ=4`: temporal smearing parameter
- `c=1e-4`: time persistance parameter 
- `λ=2.81`: decay parameter
- `act_ω`: a dictionary of unique activation weights 
- `inhib_ω`: a dictionary of unique inhibition weights

# References

Salet, J. M., Kruijne, W., van Rijn, H., Los, S. A., & Meeter, M. (2022). FMTP: A unifying computational framework of temporal preparation across time scales. 
Psychological Review.
"""
@concrete mutable struct FMTPModel <: AbstractFMTP
    τs
    κ
    c 
    λ
    act_ω    
    inhib_ω    
end

"""
    FMTPModel(;
        τs,
        κ,
        c, 
        λ, 
        act_ω,
        inhib_ω)

Constructor for the FMPTModel object. 

# Arguments

- `τs`: a vector containing the time of at maximum fire rate for each timing cell``
- `κ`: temporal smearing parameter
- `c`: time persistance parameter 
- `λ`: decay parameter
- `act_ω`: a dictionary of unique activation weights 
- `inhib_ω`: a dictionary of unique inhibition weights
"""
function FMTPModel(;
    τs,
    κ =4,
    c = 1e-4, 
    λ = 2.81, 
    act_ω,
    inhib_ω)
    return FMTPModel(τs, κ, c, λ, act_ω, inhib_ω)
end