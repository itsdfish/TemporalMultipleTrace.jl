module TemporalMultipleTrace
    using ConcreteStructs, HCubature
    export FMTPModel, cell_activation_func, compute_weight
    export precompute_weights, decay_func


    include("structs.jl")
    include("functions.jl")

# Write your package code here.

end
