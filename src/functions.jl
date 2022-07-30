"""
    cell_activation_func(t, τ, κ)

Computes the activation of a timing cell with peak firing rate at time τ.

# Arguments

- `t`: time 
- `τ`: time at peak firing rate 
- `κ`: temporal smearing parameter 
"""
function cell_activation_func(t, τ, κ)
    x1 = (1 / t)
    x2 = (κ^(κ+1)) / factorial(κ)
    x3 = (-t / -τ)^(κ+1)
    x4 = exp(-κ * (-t / -τ))
    return  x1 * x2 * x3 * x4 
end

"""
    compute_weight(tmin, tmax, τ, κ)

Computes the weight of the time unit by integrating the activity of time cell from time tmin to tmax.
Typically, tmin = 0 and tmax = fp - .050 for inhibition weights, and tmin = fp and tmax = fp + 0.30 for 
activation weights. 

# Arguments

- `tmin`: lower bound of integration over t
- `tmax`: upper bound of integration over t 
- `τ`: time at peak firing rate 
- `κ`: temporal smearing parameter 
"""
function compute_weight(tmin, tmax, τ, κ; kwargs...)::Float64
    f(t) = cell_activation_func(t[1], τ, κ)
    v,_ = hcubature(f, [tmin], [tmax]; kwargs...)
    return v
end

"""
    decay_func(i, r, c)

Computes decay of the trace computed i trials ago. 

# Arguments

- `i`: number of trial ago that the trace was formed
- `λ`: decay rate parameter
- `c`: memory persistance parameter
"""
decay_func(i, λ, c) = i^-λ + c

"""
    precompute_weights(τs, fps)

Precomputes dictionaries of activation and inhibition weights for a set of timing cells at each 
unique for period. Each key of the dictionary coresponds a forepriod, and each value corresponds 
to a vector of weights (one per timing cell). Weights are passed to the model object where they are 
cashed and referenced for efficiency. 

# Arguments

- `τs`: a vector time at peak firing rates for the timing cells
- `κ`: temporal smearing parameter 
- `ufps`: a vector of unique forepriods
"""
function precompute_weights(τs, κ, ufps)
    # activation weights 
    aws = Dict(fp => compute_weight.(fp, fp + .30, τs, κ) for fp in ufps)
    # inhibition weights 
    iws = Dict(fp => compute_weight.(0, fp - .05, τs, κ) for fp in ufps)
    return aws, iws
end
