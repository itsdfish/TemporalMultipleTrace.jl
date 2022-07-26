using SafeTestsets

@safetestset "cell_activation_func" begin
    using Test, TemporalMultipleTrace
    t = range(.1, 5, length=10)
    τ = 2.0
    κ = 4

    # taken from Python
    # t  = np.linspace(.1, 5, 10)
    # tau_star = 2.0
    # k = 4
    # fact = np.math.factorial
    # A = ( (1/t) * (k ** (k+1))/fact(k) * (-t/-tau_star)**(k+1) * 
    #                 np.exp(-k*(-t/-tau_star)) )
    # A
    true_vals = [1.09164100e-04, 6.33757834e-02, 2.47085945e-01, 3.75762026e-01,
    3.77165103e-01, 2.99193731e-01, 2.03934030e-01, 1.25023429e-01,
    7.08760990e-02, 3.78332748e-02]

    vals = cell_activation_func.(t, τ, κ)
    @test true_vals ≈ vals

    t = range(.1, 5, length=10)
    τ = 4.0
    κ = 3

    # taken from Python
    true_vals = [4.89239729e-05, 8.70444842e-03, 3.63305268e-02, 7.48440437e-02,
    1.12903685e-01, 1.42761391e-01, 1.61101874e-01, 1.67905060e-01,
    1.65018258e-01, 1.55024204e-01]

    vals = cell_activation_func.(t, τ, κ)
    @test true_vals ≈ vals

    t = range(1, 3, length=200)
    τ = 2.0
    κ = 4

    vals = cell_activation_func.(t, τ, κ)
    _,mi = findmax(vals)
    @test t[mi] ≈ τ atol = 1e-2


    t = range(3, 5, length=200)
    τ = 4.0
    κ = 2

    vals = cell_activation_func.(t, τ, κ)
    _,mi = findmax(vals)
    @test t[mi] ≈ τ atol = 1e-2
end

@safetestset "compute_weight" begin
    using Test, TemporalMultipleTrace

    # tests based on python code
    # Temporal smear 
    # k = 4 

    # # Forgetting curve    
    # r = -2.81 # rate of forgetting
    # c = 1e-4 # memory persistence
    # N = 5 # number of time cells

    # # Initialize fMTP class
    # fmtp = fMTP(r, c, k, N=N)

    # # Define FPs
    # FP = np.array([.5,])

    # # Run experiment using object "exp" and "fmtp"
    # #state_discr, state_con = exp.run_exp(fmtp) 
    # fmtp.trace_formation(FP)
    # # weights at tau = .05
    # fmtp.W[0]
    # # weights at tau = 1.2875
    # fmtp.W[1,0,:]

    fp = .5
    τ = .05
    κ = 4

    # inhabition weight 
    iw = compute_weight(0, fp - .05, τ, κ)
    @test iw ≈ 1.0

    # activation weight 
    aw = compute_weight(fp, fp +.3, τ, κ)
    @test aw ≈ 5.2039024e-13 atol = 1e-5

    fp = .5
    τ = 1.2875
    κ = 4

    # inhabition weight 
    iw = compute_weight(0, fp -.05, τ, κ)
    @test iw ≈ 0.01411575 atol = 1e-3

    # activation weight 
    aw = compute_weight(fp, fp +.3, τ, κ)
    @test aw ≈ 0.08555665 atol = 1e-3
end

@safetestset "motor_prep_func" begin
    using Test, TemporalMultipleTrace

    # # Temporal smear 
    # k = 4 

    # # Forgetting curve    
    # r = -2.81 # rate of forgetting
    # c = 1e-4 # memory persistence
    # N = 5 # number of time cells
    # dt = .00001 # increase precision

    # # Initialize fMTP class
    # fmtp = fMTP(r, c, k, N=N, dt=dt)

    # # Define FPs
    # FP = np.array([.6,.5,.5,.7,.6])                        

    # # Run experiment using object "exp" and "fmtp"
    # #state_discr, state_con = exp.run_exp(fmtp) 
    # fmtp.trace_formation(FP)
    # A,I = fmtp.trace_expression(FP)
    # fmtp.prep

    # based on the Python code above
    true_values = [NaN, .27214065, 0.19634961, 0.17796459, 0.34751846]

    τs = range(.05, 5, length=5)
    κ = 4
    ufps = [.7,.5,.6]
    λ = 2.81
    c = 1e-4

    act_ω,inhib_ω = precompute_weights(τs, κ, ufps)

    model = FMTPModel(;τs, κ, λ, c, act_ω, inhib_ω)

    fps = [.6,]
    t = .5
    push!(fps, t)
    prep = motor_prep_func(model, t, fps)
    @test prep ≈ true_values[2] atol = 1e-4

    t = .5
    push!(fps, t)
    prep = motor_prep_func(model, t, fps)
    @test prep ≈ true_values[3] atol = 1e-4

    t = .7
    push!(fps, t)
    prep = motor_prep_func(model, t, fps)
    @test prep ≈ true_values[4] atol = 1e-4

    t = .6
    push!(fps, t)
    prep = motor_prep_func(model, t, fps)
    @test prep ≈ true_values[5] atol = 1e-4
end

@safetestset "motor_preps_func" begin
    using Test, TemporalMultipleTrace

    # based on the Python code above
    true_values = [.27214065, 0.19634961, 0.17796459, 0.34751846]

    τs = range(.05, 5, length=5)
    κ = 4
    ufps = [.7,.5,.6]
    λ = 2.81
    c = 1e-4

    act_ω,inhib_ω = precompute_weights(τs, κ, ufps)

    model = FMTPModel(;τs, κ, λ, c, act_ω, inhib_ω)

    fps = [.6,.5,.5,.7,.6]
    prep = motor_preps_func(model, fps)
    @test prep ≈ true_values atol = 1e-4

    # based on the Python code above
    true_values = [0.16330973, 2.37270496, 2.73363641, 1.18045448,
    6.31565303, 3.27070813]
    fps = [.5,1.5,1.5,.7,2.3,.3,.5]

    τs = range(.05, 5, length=5)
    κ = 4
    ufps = unique(fps)
    λ = 2.81
    c = 1e-4

    act_ω,inhib_ω = precompute_weights(τs, κ, ufps)

    model = FMTPModel(;τs, κ, λ, c, act_ω, inhib_ω)

    prep = motor_preps_func(model, fps)
    @test prep ≈ true_values atol = 1e-4


    # based on the Python code above
    true_values = [0.25434547,  2.75502982,  8.17739472,  1.14787027,
    30.69833091,  1.57093889]
    fps = [.5,1.5,1.5,.7,2.3,.3,.5]

    τs = range(.05, 5, length=50)
    κ = 4
    ufps = unique(fps)
    λ = 2.81
    c = 1e-4

    act_ω,inhib_ω = precompute_weights(τs, κ, ufps)

    model = FMTPModel(;τs, κ, λ, c, act_ω, inhib_ω)

    prep = motor_preps_func(model, fps)
    @test prep ≈ true_values atol = 1e-3


    # based on the Python code above
    true_values = [4.93008748, 88.97009508,  4.33771663,  9.12580015,
    17.82549668]
    fps = [2.4,3.0,.5,4.0,3.2,1.0]

    τs = range(.05, 5, length=50)
    κ = 3
    ufps = unique(fps)
    λ = 1.5
    c = 1e-4

    act_ω,inhib_ω = precompute_weights(τs, κ, ufps)

    model = FMTPModel(;τs, κ, λ, c, act_ω, inhib_ω)

    prep = motor_preps_func(model, fps)
    @test prep ≈ true_values atol = 1e-3
end