# TemporalMultipleTrace

This package provides a Julia implementation of the formalized multiple trace theory of temporal preparation (fMTP) as described in fMTP: A Unifying Computational Framework of Temporal
Preparation Across Time Scales. The model is basically a neural network consisting of layer of time sensitive nodes which are connected to an output layer consisting of an inhibition node and an activation node. Mean motor response time is a linear function of the ratio of inhibition to activation. The weights between time nodes and motor nodes are based on memory traces containing activation and inhibition of previous experiences. The influence of old memory traces decay across time. These dynamics allow the model to produce foreperiod effects, sequential effects, and distribution effects. 

# Installation 

To install the package, enter the following in to the REPL:

```julia 
] add https://github.com/itsdfish/TemporalMultipleTrace.jl
```

# Basic Usage

The code block below provides an example of basic usage. 
```julia 
using TemporalMultipleTrace

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
# memory persistance (decay asymptote)
c = 1e-4

# cache activation and inhibition weights 
act_ω,inhib_ω = precompute_weights(τs, κ, ufps)
# model object
model = FMTPModel(;τs, κ, λ, c, act_ω, inhib_ω)
# motor preparation values (ratio of inhabition to activation)
prep = motor_preps_func(model, fps)
```

# References

Salet, J. M., Kruijne, W., van Rijn, H., Los, S. A., & Meeter, M. (2022). FMTP: A unifying computational framework of temporal preparation across time scales. Psychological Review.
