using Revise
using .ManipulatingDistributions
using  PiecewiseOrthogonalPolynomials, Test #, ManipulatingDistributions
includet("../src/legendre.jl")



@testset "FFT" begin
    @testset "sin and cos convolution" begin
        n = 1000
        h = fft_conv(sin, cos, n)
        θ = range(0, 2π, length=2n+1)[1:end-1]
        @test h ≈ π*sin.(θ) # π*sin(θ) == ∫_0^2π sin(t)cos(θ-t) dt
    end
end

@testset "Legendre series" begin
    f₁ = x -> 0.2sin(3x) + 7exp(0.4x)
    h₁ = x-> (legendreseries(f₁, N=100))[x]
    xx = range(-1, 1, length=10)
    @test f₁.(xx) ≈ h₁.(xx)

    # @test legendreseries(f, N=10) ≈ [3; 2; zeros(8)]
    #@test legendreseries(f).[x]    
end

@testset "gamma left" begin
    # f₁ = x -> x^2
    # g₁ = x -> x + 1
    # γ_left_true = [8/15, 2/9, 2/105, 0, 4/945]
    # g_lam = k -> gammaleft(k, f, g, N=10)
    # γ_exp = g_lam.(0:4)
    # @test γ_left_true ≈ γ_exp
    @test gammaleft(0, x->1, x->1, N=10) ≈ 1
    f = x -> x^2
    g = x -> x + 1
    γ_left_true = [4/15, 6/18, 10/210, 0, 36/1890]
    g_lam = k -> gammaleft(k, f, g, N=10)
    γ_exp = g_lam.(0:4)
    @test γ_left_true ≈ γ_exp
    
end

@testset "gamma right" begin
    @testset "same function" begin
        @test gammaright(0, x->1, x->1, N=10) ≈ 1
    
    #     f = x -> x^2
    #     g = x -> x^2
    #     γ_right_true = [1/9, -1/21, 8/63, -4/27, -4/105, -4/945]
    #     g_lam_r = k -> gammaright(k, f ,g, N=10)
    #     γ_exp = g_lam_r.(0:5)
    #     # g_lam = k -> gammaright(k, f, g, N=10)
    #     # γ_exp = g_lam.(0:4)
    #     @test γ_right_true ≈ γ_exp
     end

    @testset "different function" begin
        f = x -> x^4
        g = x -> 1
        γ_right_true = [1/5, -3/35, 0, -4/45, 0]
        g_lam = k -> gammaright(k, f, g, N=10)
        γ_exp = g_lam.(0:5)
        println(γ_exp)
        println(γ_right_true)
        @test γ_right_true ≈ γ_exp
        
        # f = x -> x^2
        # g = x -> x + 1
        # γ_right_true = [2/5, -1/15, -1/21, -4/15, -2/105]
        # g_lam = k -> gammaright(k, f, g, N=10)
        # γ_exp = g_lam.(0:4)
        # @test γ_right_true ≈ γ_exp
    end
end





    # @test legendreseries(f, N=10) ≈ [3; 2; zeros(8)]
    #@test legendreseries(f).[x]    


# @testset "Legendre" begin
#     n = 10
#     h = legendre_conv(one, one, [0,1], [0,1]; N=100)

#     @test_broken iszero(h(-0.5))
#     @test_broken h(0.1) ≈ 0.1 # x
#     @test_broken h(1.1) ≈ 1.1 # 2-x
#     @test iszero(h(2.5))


#     P = ContinuousPolynomial{0}([0, 1, 2]) # Piecewise Legendre
#     x = axes(P,1)
#     c = P \ broadcast(x -> if 0 ≤ x ≤ 1
#                             x
#                         elseif 1 ≤ x ≤ 2
#                             2-x
#                         else
#                             zero(x)
#                         end, x)
#     h̃ = P * c

#     @test_broken h == h̃

#     @testset "examples of ContinuousPolynomial" begin
#         Q = ContinuousPolynomial{1}([0,1,2])
        
#         @test Q[0.5,1] == 1-0.5
#         @test Q[1.5,1] == 0

#         @test Q[0.5,4] == 4 * 0.5*(1-0.5)

#         Q = ContinuousPolynomial{1}([-1,1])
#         x = 0.5
#         @test Q[x,1] == (1 - x)/2
#         @test Q[x,2] == (1 + x)/2
#         for k = 3:10
#             @test Q[x,k] ≈ (1-x^2) * jacobip(k-3, 1, 1, x) # Degree k-1
#         end

#         # (1-x^2) * P_n^(1,1)(x) are orthogonal wrt
#         # <f,g>_1 := <f',g'> == ∫_{-1}^1 f'(x) g'(x) dx

#         P = ContinuousPolynomial{0}([0,1,2])
#         c = [randn(10); zeros(∞)]
#         g = Q * c
#         c̃ = (P \ Q) * c
#         g̃ = P * c̃
#         @test g[0.1] ≈ g̃[0.1]
#     end
# end