###################################################################################################
#                                 load packages
###################################################################################################
cd(@__DIR__)
using Pkg
Pkg.activate("")
using Revise, TemporalMultipleTrace

# sample foreperiods from a discrete uniform distribution
fps = rand([.2,4,.4], 20)
# unique foreperiods
ufps = unique(fps)
# a vector of time points corresponding to maximum fire rate of each time cell 
τs = range(.05, 5, length=50)
# temporal smearing parameter (controls width of activation function)
κ = 4
# trace decay rate 
λ = 2.81
# memory persistance (lower bound asymptote)
c = 1e-4

# cache activation and inhibition weights 
act_ω,inhib_ω = precompute_weights(τs, κ, ufps)
# model object
model = FMTPModel(;τs, κ, λ, c, act_ω, inhib_ω)
# motor preparation values (ratio of inhabition to activation)
prep = motor_preps_func(model, fps)

