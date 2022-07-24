using Revise, ManipulatingDistributions, PiecewiseOrthogonalPolynomials, Test
using Distributions
using ClassicalOrthogonalPolynomials

import ManipulatingDistributions: gammaleft, gammaright, legendreseries,
                                  gammaleft_matrix, gammaright_matrix


@testset "FFT" begin
    @testset "sin and cos convolution" begin
        n = 1000
        h = fft_conv(sin, cos, n)
        θ = range(0, 2π, length=2n+1)[1:end-1]
        @test h ≈ π*sin.(θ) # π*sin(θ) == ∫_0^2π sin(t)cos(θ-t) dt
    end

    @testset "normal distribution" begin
        n=10
        lower= 0
        upper = 2π
        d₁ = Normal(2, 0.5)
        d₂ = Normal(1.5, 0.5)
        # d₁ = Normal(0, 1)
        # d₂ = Normal(0, 1)

        #f_norm = Normal(0, sqrt(2))
        f_norm = Normal(3.5, 0.5*sqrt(2))
        # y₁ = pdf.(d, x₁)
        # y₂ = pdf.(d, x₂)
        # scatter(x, y)

        norm₁ = x -> pdf(d₁, x)
        norm₂ = y -> pdf(d₂, y)

        # fourierconv(norm₁, norm₁, n)
        z = range(lower, upper, length=2n+1)[1:end-1]

        data = abs.(fft_conv(norm₁, norm₂, n, lower=lower, upper=upper))
        

        #@test data ≈ pdf.(f_norm, z)
        @test isapprox(data,  pdf.(f_norm, z), atol=0.001)
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


@testset "gamma left matrix" begin
    f = x -> x^2
    g = x -> x + 1
    γ_left_true = [4/15, 6/18, 10/210, 0, 36/1890]
    γ_exp = gammaleft_matrix(f, g; α_s=2, β_s=1)
    @test γ_left_true ≈ γ_exp

    f = x -> (2x+1)^3
    g = x -> x^2 + 5x - 9
    γ_left_true = [-(365/21), -(292/7), -(2381/63), -(844/45), -(1964/385), 64/315, 32/3465]
    γ_exp = gammaleft_matrix(f, g; α_s=3, β_s=2)
    @test γ_left_true ≈ γ_exp


end

@testset "gamma right matrix" begin
    @testset "same function" begin
        @test gammaright_matrix(x->1, x->1, α_s=0, β_s=0) ≈ [1]
        f = x -> x^2
        g = x -> x^2
        γ_right_true = [1/9, -1/21, 8/63, -4/27, -4/105, -4/945]
        γ_exp = gammaright_matrix(f, g, α_s=2, β_s=2)
        @test γ_right_true ≈ γ_exp
    end

    @testset "different function" begin
        # f₁ = x -> 1
        # g₁ = x -> x^3
        # γ_right_true = [1/5, 0, -1/7, 0, -2/35]
        # γ_exp = gammaright_matrix(f₁, g₁, α_s=0, β_s=3)
        # @test γ_right_true ≈ γ_exp

        # f₁ = x -> x+1
        # g₁ = x -> x^3
        # γ_right_true = [1/5, 8/35, -2/7, -1/45, -4/35, -2/315]
        # γ_exp = gammaright_matrix(f₁, g₁, α_s=1, β_s=3)
        # @test γ_right_true ≈ γ_exp

        # f₁ = x -> (2x+1)^3
        # g₁ = x -> x^3 + 3
        # γ_right_true = [186/7, -743/105, -377/35, -862/165, -(1194/385),
        #                  -(484/1365), -(16/385), -(32/15015)]  
        # γ_exp = gammaright_matrix(f₁, g₁, α_s=3, β_s=3)
        # @test γ_right_true ≈ γ_exp

        f₁ = x -> (3x +2)^5
        g₁ = x -> (x+4)^4
        γ_right_true = [39434387/165, 7095810/77, -(4994755/39), -(1110657260/9009), -(295610936/5005), -(5556808/315),
                        -(11500480/3927), -(12160/77), -(9600/1729), -(576/5005), -(1728/1616615)]
        γ_exp = gammaright_matrix(f₁, g₁, α_s=5, β_s=4)
        @test γ_right_true ≈ γ_exp

        f₂ = x -> x^2
        g₂ = x -> x + 1
        γ_right_true = [2/5, -1/15, -1/21, -4/15, -2/105]
        γ_exp = gammaright_matrix(f₂, g₂, α_s=2, β_s=1)
        @test γ_right_true ≈ γ_exp
    end
end


# @testset "legendre convolution" begin
#     @testset "same interval" begin
#         n = 1000
#         h = legendre_conv(sin, cos, [-1, 1], [-1, 1], N = 3)
#         θ = range(-1, 1, length=2n+1)[1:end-1]
#         @test h.(θ) ≈ π*sin.(θ)
#     end

#     @testset "general interval" begin
        
#     end
# end



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