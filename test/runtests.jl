using ManipulatingDistributions, PiecewiseOrthogonalPolynomials, Test

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