###################################################################################################
#                                 load packages
###################################################################################################
cd(@__DIR__)
using Pkg
Pkg.activate("")
using Revise, TemporalMultipleTrace

fps = repeat([.2,2.0,.4], inner = 166)

τs = range(.05, 5, length=50)
κ = 4
ufps = unique(fps)
λ = 2.81
c = 1e-4

act_ω,inhib_ω = precompute_weights(τs, κ, ufps)

model = FMTPModel(;τs, κ, λ, c, act_ω, inhib_ω)

prep = motor_preps_func(model, fps)

@time motor_preps_func(model, fps);