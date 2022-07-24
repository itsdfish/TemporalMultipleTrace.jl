abstract type AbstractFMTP end 

"""
    FMTPModel

Model object for formal multiple trace theory of temporal preparation.

# Fields

- `τs`: a vector containing the time of at maximum fire rate for each timing cell``
- `τ_min`: minium τ
- `τ_max`: maximum τ
- `Δt`: time step
- `κ`: temporal smearing parameter
- `c`: time persistance parameter 
- `λ`: decay parameter
- `n_cells`: number of temporal cells

# References

Salet, J. M., Kruijne, W., van Rijn, H., Los, S. A., & Meeter, M. (2022). FMTP: A unifying computational framework of temporal preparation across time scales. Psychological Review.
"""
@concrete mutable struct FMTPModel <: AbstractFMTP
    τs
    τ_min
    τ_max
    Δt 
    κ
    c 
    λ
    n_cells
end

"""
    FMTPModel(;
        τs,
        τ_min,
        τ_max,
        Δt, 
        κ,
        c, 
        λ,
        n_cells)

Constructor for the FMPTModel object. 

# Arguments

- `τs`: a vector containing the time of at maximum fire rate for each timing cell``
- `τ_min`: minium τ
- `τ_max`: maximum τ
- `Δt`: time step
- `κ`: temporal smearing parameter
- `c`: time persistance parameter 
- `λ`: decay parameter
- `n_cells`: number of temporal cells
"""

function FMTPModel(;
    τs,
    τ_min,
    τ_max,
    Δt, 
    κ,
    c, 
    λ,
    n_cells)
    return FMTPModel(τs, τ_min, τ_max, Δt, κ, c, λ, n_cells)
end

abstract type AbstractTask end

@concrete mutable struct Task
    dist
    n_trials 
end