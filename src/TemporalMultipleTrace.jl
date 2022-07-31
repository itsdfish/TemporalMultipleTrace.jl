module TemporalMultipleTrace
    using ConcreteStructs
    using HCubature
    
    export FMTPModel
    export cell_activation_func
    export compute_weight
    export precompute_weights
    export decay_func
    export motor_prep_func
    export motor_preps_func

    include("structs.jl")
    include("functions.jl")
end
