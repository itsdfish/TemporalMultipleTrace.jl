module TemporalMultipleTrace
    using ConcreteStructs, HCubature
    export FMTPModel, cell_activation_func, compute_weight
    export precompute_weights, decay_func, motor_prep_func
    export motor_preps_func

    include("structs.jl")
    include("functions.jl")
end
