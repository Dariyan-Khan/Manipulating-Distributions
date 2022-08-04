using Revise, ManipulatingDistributions, PiecewiseOrthogonalPolynomials, Test
using Distributions, LinearAlgebra
using ClassicalOrthogonalPolynomials

import ManipulatingDistributions: gammaleft, gammaright, legendreseries,
                                  gammaleft_matrix, gammaright_matrix, bleft_matrix,
                                  legendrecoeff, bright_matrix, legendre_same_length,
                                  conv_unit, left_conv_unit, right_conv_unit


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

# @testset "gamma left" begin
#     # f₁ = x -> x^2
#     # g₁ = x -> x + 1
#     # γ_left_true = [8/15, 2/9, 2/105, 0, 4/945]
#     # g_lam = k -> gammaleft(k, f, g, N=10)
#     # γ_exp = g_lam.(0:4)
#     # @test γ_left_true ≈ γ_exp
#     @test gammaleft(0, x->1, x->1, N=10) ≈ 1
#     f = x -> x^2
#     g = x -> x + 1
#     γ_left_true = [4/15, 6/18, 10/210, 0, 36/1890]
#     g_lam = k -> gammaleft(k, f, g, N=10)
#     γ_exp = g_lam.(0:4)
#     @test γ_left_true ≈ γ_exp

#     f = x -> x^4 + 2x + 5
#     g = x -> x^5 + 3x^4 + x^2 +5x + 7
#     γ_left_true = [52116/1925, 9674/275, 157666/15015, 525544/225225, 
#                         23396/75075, 1192/2925, 3776/98175, 3616/765765,
#                         -(64/611325), -(128/3828825), 64/14549535]
#     g_lam = k -> gammaleft(k, f, g, N=20)
#     γ_exp = g_lam.(0:10)
#     @test γ_left_true ≈ γ_exp
# end

@testset "Bleft matrix" begin
    f = x -> x^4 + 2x + 5
    g = x -> x^5 + 3x^4 + x^2 +5x + 7
    α_s = 5
    β_s = 6
    α = legendrecoeff(f, N=α_s + β_s + 4)
    B_left = bleft_matrix(α, α_s=α_s, β_s=β_s)

    B_actual = [4.533333333333331       -1.6952380952380965     0.13333333333333203     -0.012698412698412829   9.853815308202007e-16   -0.0023088023088021557  2.9219050110271733e-16
    5.085714285714289       -0.6857142857142899     -0.6095238095238138     0.019047619047618057    -0.0034632034632013763  -0.006926406926405424   -0.0005328005328000812
    0.6666666666666605      1.015873015873023       -0.03809523809524781    -0.49639249639250127    0.01443001443001645     -0.012876012876009355   -0.002664002664001851
    0.08888888888889027     0.044444444444441795    0.6949494949495016      -0.008080808080822155   -0.38414918414918636    -0.005283605283600082   -0.008080808080806029
    7.8130018422289e-15     0.010389610389604412    0.025974025974029494    0.49390609390609663     0.023176823176807384    -0.33086913086912795    -0.00759240759240333
    0.025396825396818676    -0.02539682539681831    0.028327228327220254    -0.008302808302800046   0.4043956043956006      0.022466422466410393    -0.2883286647992487
    3.424948443775583e-16   0.002308802308801455    -0.006926406926404797   0.015007215007210402    -0.010966810966804278   0.3407520583991118      0.018538324420667576
    -3.9244851485288357e-16 5.663381600147883e-16   0.0005328005328003373   -0.002664002664002294   0.008540479128712473    -0.009598244892357378   0.29346092503986626
    2.019570004747791e-15   -4.839568725926056e-16  6.837529970273157e-17   0.00017760017760049072  -0.0012432012432007727  0.005267247372508485    -0.007900871058760939
    1.4038196738666708e-15  -7.216417072689524e-16  2.726241851317163e-16   -1.560369050955016e-16  7.312948489416502e-05   -0.0006581653640481252  0.003463203463203595
    1.9527592505349253e-15  7.309257366028152e-16   -7.395075446776022e-16  2.8857125022873125e-16  -1.98000590674196e-16   3.464028231839639e-05   -0.00038104310550083676
    4.2108383388693185e-15  -1.1353226150092589e-15 4.93603174207213e-16    -6.406997100216386e-16  6.208513366111186e-16   -1.296749896734384e-16  1.814490978562449e-05
    9.674825181169094e-15   -3.1197143083147296e-15 5.780920420057809e-16   -1.411011684477542e-16  3.162128953121657e-16   -1.8532112848269778e-16 -1.4689572467013252e-16]
    
    #@test isapprox(B_actual, B_left, atol=0.1) #B_actual ≈ B_left
    #print(opnorm(B_actual - B_left, Inf))
    @test B_actual ≈ B_left
end

@testset "Bright matrix" begin
    f = x -> x^4 + 2x + 5
    g = x -> x^5 + 3x^4 + x^2 +5x + 7
    α_s = 5
    β_s = 6
    α = legendrecoeff(f, N=α_s + β_s + 4)
    B_right = bright_matrix(α, α_s=α_s, β_s=β_s)

    B_actual = [5.866666666666664       1.6952380952380965      -0.13333333333333194    0.012698412698412905    -5.253034539442149e-16  0.002308802308801759    9.43130274373311e-18
    -5.085714285714289      0.9142857142857067      1.4095238095238056      -0.09523809523809461    0.003463203463205556    -0.006926406926405403   0.0005328005328004862
    -0.6666666666666589     -2.349206349206344      0.3428571428571413      1.0678210678210636      -0.04906204906205289    0.012876012876011202    -0.0026640026640022147
    -0.08888888888889065    -0.22222222222222       -1.494949494949489      0.16969696969697398     0.8285936285936324      -0.045687645687645044   0.008080808080807102
    -3.3359611554976743e-15 -0.01038961038961756    -0.08831168831169522    -1.0653346653346698     0.1270729270729235      0.6945054945054958      -0.03556443556443399
    -0.025396825396821983   -0.025396825396819787   -0.028327228327224424   -0.07179487179487097    -0.8488400488400496     0.09084249084248691     0.5960209724915602
    3.3806678584396337e-15  -0.002308802308801617   -0.006926406926405417   -0.015007215007212391   -0.051370851370849624   -0.7043884220354806     0.0670231729055224
    1.4561753125174758e-15  -3.2682066400649204e-16 -0.0005328005327999742  -0.002664002664002189   -0.00854047912871351    -0.03757027286438865    -0.6011532327321787
    -4.5485394934285445e-15 1.0921942913569594e-15  -3.187852511003292e-16  -0.00017760017760090299 -0.0012432012432005025  -0.005267247372509637   -0.028413691571584906
    1.0989688534422927e-14  -1.5965684228694373e-15 -1.4720058009811008e-15 6.433035137675858e-16   -7.312948489458878e-05  -0.0006581653640477087  -0.0034632034632035343
    -1.805225443781307e-14  -1.9697785069717904e-16 1.120560747160094e-15   -8.381316474181873e-16  1.3766532609290025e-16  -3.464028231899424e-05  -0.0003810431055009706
    9.746632897029716e-15   -1.3282795468361282e-15 -3.807470212129644e-16  5.662423145933649e-16   -1.0682860845934043e-15 1.9048498566772826e-16  -1.814490978571768e-05
    1.4962757238543455e-15  3.3743469328241673e-15  -5.649692219275219e-16  -7.786120458692245e-16  3.4042667070763786e-16  4.4200599281864556e-16  4.769893848909057e-16]
    #@test isapprox(B_actual, B_left, atol=0.1) #B_actual ≈ B_left
    #print(opnorm(B_actual - B_left, Inf))
    @test B_actual ≈ B_right
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

    f = x -> x^4 + 2x + 5
    g = x -> x^5 + 3x^4 + x^2 +5x + 7
    γ_left_true = [52116/1925, 9674/275, 157666/15015, 525544/225225, 
                   23396/75075, 1192/2925, 3776/98175, 3616/765765,
                   -(64/611325), -(128/3828825), 64/14549535, 0, 0]
    γ_exp = gammaleft_matrix(f, g; α_s=5, β_s=6)

    @test γ_left_true ≈ γ_exp

    γ_left_true = [-0.2726756433, -0.05775281530, 0.2859244576, 0.05956759562,
                   -0.01343585768, -0.001833723566, 0.0001882724827, 0.00001904201277,
                     0]
    γ_exp = gammaleft_matrix(sin, cos; α_s=3, β_s=4)

    @test isapprox(γ_left_true,  γ_exp, atol=0.001)

    
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
        f₁ = x -> 1
        g₁ = x -> x^3
        γ_right_true = [1/5, 0, -1/7, 0, -2/35]
        γ_exp = gammaright_matrix(f₁, g₁, α_s=0, β_s=3)
        @test γ_right_true ≈ γ_exp

        f₁ = x -> x+1
        g₁ = x -> x^3
        γ_right_true = [1/5, 8/35, -2/7, -1/45, -4/35, -2/315]
        γ_exp = gammaright_matrix(f₁, g₁, α_s=1, β_s=3)
        @test γ_right_true ≈ γ_exp

        f₁ = x -> (2x+1)^3
        g₁ = x -> x^3 + 3
        γ_right_true = [186/7, -743/105, -377/35, -862/165, -(1194/385),
                         -(484/1365), -(16/385), -(32/15015)]  
        γ_exp = gammaright_matrix(f₁, g₁, α_s=3, β_s=3)
        @test γ_right_true ≈ γ_exp

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

@testset "left convolution" begin
    h_series = left_conv_unit(sin, cos, α_s=10, β_s=10)
    h_lam = x -> h_series[x]
    n=1000
    θ = range(-2,0, length=2n+1)[1:end-1]
    h_actual = x -> 0.5 * (2+x) * sin(x)
    @test h_lam.(θ) ≈ h_actual.(θ)
end

@testset "right convolution" begin
    h_series = right_conv_unit(sin, cos, α_s=10, β_s=10)
    h_lam = x -> h_series[x]
    n=1000
    θ = range(0,2, length=2n+1)[1:end-1]
    h_actual = x -> -0.5 * (x-2) * sin(x)
    @test h_lam.(θ) ≈ h_actual.(θ)
end

@testset "[-1, 1] convolution" begin
    n = 1000
    h_lam = x -> conv_unit(x->sin(x), x->cos(x), x, α_s = 10, β_s = 11)
    θ = range(-2, 2, length=2n+1)[1:end-1]

    function h_trig(x::Real)
        if x in -2..0
            return 0.5 * (2+x) * sin(x)
        elseif x in 0..2
            return -0.5 * (x-2) * sin(x)
        else
            return 0
        end
    end

    @test h_lam.(θ) ≈ h_trig.(θ)
    
end


@testset "same interval convolution" begin
    n = 1000
    h = legendre_same_length(sin, cos, [-1, 1], [-1, 1], α_s = 10, β_s = 11)
    θ = range(-2, 2, length=2n+1)[1:end-1]

    function h_trig(x::Real)
        if x in -2..0
            return 0.5 * (2+x) * sin(x)
        elseif x in 0..2
            return -0.5 * (x-2) * sin(x)
        else
            return 0
        end
    end
    @test h.(θ) ≈ h_trig.(θ)
    
end

@testset "legendre convolution" begin
    @testset "same int polynomial" begin
        f = x -> x^2
        g = x -> x^3

        function h_mono(x::Real)
            if x in -2..0
                return x/5 + (x^2)/2 + (x^3)/3 + x^6/60
            elseif x in 0..2
                return x/5 - (x^2)/2 + (x^3)/3 - x^6/60
            else
                return 0
            end
        end

        n = 1000
        h = legendre_conv(f, g, [-1, 1], [-1, 1], α_s = 3, β_s = 4)
        θ = range(-1, 1, length=2n+1)[1:end-1]
        @test h.(θ) ≈ h_mono.(θ)
    end

    @testset "same interval" begin
        n = 1000
        h = legendre_conv(sin, cos, [-1, 1], [-1, 1], α_s = 10, β_s = 11)
        θ = range(-2, 2, length=2n+1)[1:end-1]

        function h_trig(x::Real)
            if x in -2..0
                return 0.5 * (2+x) * sin(x)
            elseif x in 0..2
                return -0.5 * (x-2) * sin(x)
            else
                return 0
            end
        end

        @test h.(θ) ≈ h_trig.(θ)
    end

    @testset "general intervals" begin
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