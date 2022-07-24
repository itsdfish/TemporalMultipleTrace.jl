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

    vals = compute_activation.(t, τ, κ)
    @test true_vals ≈ vals

    t = range(.1, 5, length=10)
    τ = 4.0
    κ = 3

    # taken from Python
    true_vals = [4.89239729e-05, 8.70444842e-03, 3.63305268e-02, 7.48440437e-02,
    1.12903685e-01, 1.42761391e-01, 1.61101874e-01, 1.67905060e-01,
    1.65018258e-01, 1.55024204e-01]

    vals = compute_activation.(t, τ, κ)
    @test true_vals ≈ vals
end
