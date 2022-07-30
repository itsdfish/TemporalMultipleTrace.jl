###################################################################################################
#                                 load packages
###################################################################################################
cd(@__DIR__)
using Pkg
Pkg.activate("")
using Revise, TemporalMultipleTrace

function inner_weight_func(λ, c, tidx, τidx, weights, fp)
    return  weights[fp][τidx] * decay_func(tidx, λ, c)
end

function trace_weight_func(λ, c, τidx, weights, fps)
    val = 0.0
    for i in 1:length(fps)
        val += inner_weight_func(λ, c, i, τidx, weights, fps[i])
    end
    return val
end

function total_weight_func(λ, c, τs, κ, t, weights, fps)
    val = 0.0
    for (i,τ) in enumerate(τs)
        val += cell_activation_func(t, τ, κ) * trace_weight_func(λ, c, i, weights, fps)
    end
    return val
end

τs = range(.05, 5, length=5)
κ = 4
ufps = [.5,.6]
λ = 2.81
c = 1e-4

uaw,uiw = precompute_weights(τs, κ, ufps)

tidx = 2
τidx = 1
fps = [.6,.5]
τidx = 1
t = .5

total_activation = total_weight_func(λ, c, τs, κ, t, uaw, fps)
total_inhibition = total_weight_func(λ, c, τs, κ, t, uiw, fps)

total_inhibition / total_activation 

