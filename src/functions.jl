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

function inner_weight_func(λ, c, tidx, τidx, weights, fp)
    return  weights[fp][τidx] * decay_func(tidx, λ, c)
end

function trace_weight_func(λ, c, τidx, weights, fps)
    val = 0.0
    n_fp = length(fps)
    for i in 1:(n_fp - 1)
        tidx = n_fp - i + 1
        val += inner_weight_func(λ, c, tidx, τidx, weights, fps[i])
    end
    return val
end

function total_weight_func(λ, c, τs, κ, t, weights, fps)
    val = 0.0
    for i in 1:length(τs)
        val += cell_activation_func(t, τs[i], κ) * trace_weight_func(λ, c, i, weights, fps)
    end
    return val
end

"""
    motor_prep_func(model::AbstractFMTP, t, fps)

Computes motor preparation at time `t` given a vector of foreperiods.

# Arguments

- `model`: `AbstractFMTP` model object 
- `fps`: a vector of forperiods
- `t`: time at stimulus onset
"""
function motor_prep_func(model::AbstractFMTP, t, fps)
    (;λ,c,τs,κ,act_ω,inhib_ω) = model
    return motor_prep_func(λ, c, τs, κ, act_ω, inhib_ω, t, fps)
end

function motor_prep_func(λ, c, τs, κ, act_ω, inhib_ω, t, fps)
    act = total_weight_func(λ, c, τs, κ, t, act_ω, fps)
    inhib = total_weight_func(λ, c, τs, κ, t, inhib_ω, fps)
    return inhib / act 
end

function motor_prep_func(λ, c, τs, κ, act_ω, inhib_ω, fps)
    return motor_prep_func(λ, c, τs, κ, act_ω, inhib_ω, fps[end], fps)
end

function motor_prep_func(model::AbstractFMTP, fps)
    return motor_prep_func(model::AbstractFMTP, fps[end], fps)
end

"""
    motor_preps_func(model::AbstractFMTP, fps)

Computes motor preparation for a vector of foreperiods.

# Arguments

- `model`: `AbstractFMTP` model object 
- `fps`: a vector of forperiods
"""
function motor_preps_func(model::AbstractFMTP, fps)
    (;λ,c,τs,κ,act_ω,inhib_ω) = model
    func(t) = motor_prep_func(λ, c, τs, κ, act_ω, inhib_ω, @view fps[1:t])
    return map(func, 2:length(fps))
end