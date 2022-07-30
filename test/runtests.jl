using SafeTestsets

@safetestset "compute_activation" begin
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