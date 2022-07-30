"""
    compute_activation(t, τ, κ)

Computes the activation of a timing cell with peak firing rate at time τ.

# Arguments

- `t`: time 
- `τ`: time at peak firing rate 
- `κ`: temporal smearing parameter 
"""
function compute_activation(t, τ, κ)
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
    f(t) = compute_activation(t[1], τ, κ)
    v,_ = hcubature(f, [tmin], [tmax]; kwargs...)
    return v
end
