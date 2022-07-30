###################################################################################################
#                                 load packages
###################################################################################################
cd(@__DIR__)
using Pkg
Pkg.activate("")
using Revise, TemporalMultipleTrace

τs = range(.05, 5, length=5)
κ = 4
ufps = [.5,.6]
λ = 2.81
c = 1e-4

act_ω,inhib_ω = precompute_weights(τs, κ, ufps)

model = FMTPModel(;τs, κ, λ, c, act_ω, inhib_ω)

fps = [.6,.5,.5,.5,.5,.6]
t = .6

motor_prep_func(model, t, fps)

