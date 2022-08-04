"""
Calculates the Fourier Transform of the mth degree Legendre Polynomial. (3.1)
"""
function legendreft(m) 
    f = x -> 2 * (-1im)^m * sphericalbesselj(m, x)
    return f
end


"""
Calculates the first N legendre coefficients in the Legendre expansion of
a function over a specified interval.
"""
function legendrecoeff(f; interval=-1 .. 1, N=∞)
    T = legendre(interval)
    x = axes(T, 1)
    N == ∞ ? T \ f.(x) : T[:,1:N] \ f.(x)
end


"""
Returns the Legendre series of a function truncated at N terms over a specified
interval. Output will be an anonymous function.
"""
function legendreseries(f; interval=-1..1, N=∞)
    T = legendre(interval)
    c = legendrecoeff(f, interval=interval, N=N)
    N == ∞ ? T * c : T[:, 1:N] * c
end


"""
 Calculates B_{k,n}^{left} from (4.5)
"""
function bleft_matrix(α; α_s=100, β_s=100)
    if β_s == 0
        return (α[1] - (α[2]/3))
    end

    B = zeros(Float64, α_s + β_s + 3, β_s+1)
    B[1, 1] = α[1] - α[2]/3

    for k in 2:(α_s + β_s + 3)
        k₋ = k-1
        B[k, 1] = (α[k₋]/(2k₋-1)) - (α[k+1]/(2k₋ +3))
    end
    
    B[1,2] = -B[2,1]/3

    for k in 2:(α_s + β_s + 2)
        k₋ = k - 1
        B[k, 2] = (B[k₋, 1]/(2k₋-1)) - (B[k+1, 1]/(2k₋+3)) - B[k, 1]
    end

    for n in 3: (β_s +1)
        B[1, n] = B[n, 1] * (-1)^(n-1) * (1/(2*(n-1)+1))
        B[2, n] = B[n, 2] * (-1)^(n-1+1) * (3/(2*(n-1)+1))
    end

    for n in 3: (β_s + 1)
        for k in n:(α_s + β_s + 2)
            k₋ = k-1
            n₋ = n-2
            B[k, n] = (-((2n₋+1)/(2k₋+3)) * B[k+1, n-1]) + (((2n₋+1)/(2k₋-1)) * B[k-1, n-1]) + B[k, n-2]
        end
    end

    for n in 4: (β_s + 1)
        for k in 3:n
            k₋ = k-1
            n₋ = n-1
            B[k, n] = B[n, k] * (-1)^(n₋ + k₋) * ((2k₋+1)/(2n₋+1))
        end
    end

    return B[1:α_s + β_s + 2, 1:β_s+1]
end


"""
Calculates the γ_{left} coefficients for functions f ⋆ g (both supported on 
[-1,1] on the interval [0,2] expanded in the Legendre basis. This is seen in 
(4.2). Set α_s and β_s to the highest degree used in the expansion for f and
g respectively + 1.
"""
function gammaleft_matrix(f, g; α_s=100, β_s=100)
    α = legendrecoeff(f, N=α_s + β_s + 4)
    β = legendrecoeff(g, N=β_s + 1)
    return bleft_matrix(α, α_s=α_s, β_s=β_s) * β
end


"""
Calculates B^{right}, whose entries satisfy an almost identical recurence 
relation to those in  B^{left}
"""
function bright_matrix(α; α_s=100, β_s=100)
    if β_s == 0
        return (α[1] + (α[2]/3))
    end

    B = zeros(Float64, α_s + β_s + 3, β_s+1)
    B[1, 1] = α[1] + α[2]/3

    for k in 2:(α_s + β_s + 3)
        k₋ = k-1
        B[k, 1] = -(α[k₋]/(2k₋-1)) + (α[k+1]/(2k₋ +3))
    end
    
    B[1,2] = -B[2,1]/3

    for k in 2:(α_s + β_s + 2)
        k₋ = k - 1
        B[k, 2] = (B[k₋, 1]/(2k₋-1)) - (B[k+1, 1]/(2k₋+3)) + B[k, 1]
    end

    for n in 3: (β_s +1)
        B[1, n] = B[n, 1] * (-1)^(n-1) * (1/(2*(n-1)+1))
        B[2, n] = B[n, 2] * (-1)^(n-1+1) * (3/(2*(n-1)+1))
    end

    for n in 3: (β_s + 1)
        for k in n:(α_s + β_s + 2)
            k₋ = k-1
            n₋ = n-2
            B[k, n] = (-((2n₋+1)/(2k₋+3)) * B[k+1, n-1]) + (((2n₋+1)/(2k₋-1)) * B[k-1, n-1]) + B[k, n-2]
        end
    end

    for n in 4: (β_s + 1)
        for k in 3:n
            k₋ = k-1
            n₋ = n-1
            B[k, n] = B[n, k] * (-1)^(n₋ + k₋) * ((2k₋+1)/(2n₋+1))
        end
    end

    return B[1:α_s + β_s + 2, 1:β_s+1]
end


"""
Calculates the γ_{right} coefficients for functions f ⋆ g (both supported on 
[-1,1] on the interval [2,0] expanded in the Legendre basis. This is  similar 
to (4.2). Set α_s and β_s to the highest degree used in the expansion for f and
g respectively + 1.
"""
function gammaright_matrix(f, g; α_s=100, β_s=100)
    α = legendrecoeff(f, N=α_s + β_s + 4)
    β = legendrecoeff(g, N=β_s + 1)
    return bright_matrix(α, α_s=α_s, β_s=β_s) * β
end


"""
Returns the Legendre series of the left side of f ⋆ g (i.e on the interval 
[-2,0])
"""
function left_conv_unit(f, g; α_s=100, β_s=100)
    γₗ = gammaleft_matrix(f, g, α_s=α_s, β_s=β_s)
    # length is α_s + β_s + 2
    T = legendre(-2..0)
    return T[:, 1:α_s + β_s + 2] * γₗ
end


"""
Returns the Legendre series of the right side of f ⋆ g (i.e on the interval 
[0,2])
"""
function right_conv_unit(f, g; α_s=100, β_s=100)
    γᵣ = gammaright_matrix(f, g, α_s=α_s, β_s=β_s)
    return legendre(0..2)[:, 1:α_s + β_s + 2] * γᵣ
end

"""
Returns an anonymous function that combines the left and right sides of the 
convolutions.
"""
function conv_unit(f, g, x::Real; α_s=100, β_s=100)
    #println(x)
    if x in -2..0
        return left_conv_unit(f, g, α_s=α_s, β_s=β_s)[x]
    elseif x in 0..2
        return right_conv_unit(f, g, α_s=α_s, β_s=β_s)[x]
    else
        return 0
    end
end

function inv_phi_func(dom_f)
    return x  -> (((dom_f[2] - dom_f[1]) / 2) * (x + 1)) + dom_f[1]
end

"""
Returns f ⋆ g where f and g are supported on arbitrary intervals of the same
length. (4.9)
"""
function legendre_same_length(f, g, dom_f, dom_g; α_s=100, β_s=100)
    @assert dom_f[2] - dom_f[1] == dom_g[2] - dom_g[1]
    L = dom_f[2] - dom_f[1]
    #Inverse of each of the phi functions in the paper
    # ϕ_f_inv = x -> (((dom_f[2] - dom_f[1]) / 2) * (x + 1)) + dom_f[1]
    # ϕ_g_inv = x -> (((dom_g[2] - dom_g[1]) / 2) * (x + 1)) + dom_g[1]
    ϕ_f_inv = inv_phi_func(dom_f)
    ϕ_g_inv = inv_phi_func(dom_g)
    fᵣ = x -> f(ϕ_f_inv(x))
    gᵣ = x -> g(ϕ_g_inv(x))
    dom_full = [dom_f[1] + dom_g[1], dom_f[2] + dom_g[2]]
    ϕ_map = x -> (2*(x-dom_full[1])/(dom_full[2] - dom_full[1])) - 1
    # return (L/2) * left_conv_unit(fᵣ, gᵣ, N=N), (L/2) * right_conv_unit(fᵣ, gᵣ, N=N)
    return x -> (L/2) * conv_unit(fᵣ, gᵣ, 2 * ϕ_map(x), α_s=α_s, β_s=β_s)
end


"""
Calculates the function h as in (5 ii).
"""
function h_12(f, g, dom_f, dom_g,  x; α_s=100, β_s=100)
    g1 = x -> g(x)*(x in dom_g[1]..(dom_g[1] + dom_f[2] - dom_f[1]))
    g2 = x -> g(x)*(x in (dom_f[1] + dom_g[1])..(dom_g[2] - dom_f[2] + 2*dom_f[1]))
    g3 = x -> g(x)*(x in (dom_g[2] - dom_f[2] + dom_f[1])..dom_g[2])
    g4 = x -> g(x)*(x in (dom_g[2] - dom_f[2] + dom_f[1])..dom_g[2])
    f1 = x -> f(x)*(x in (2*dom_f[2] + dom_g[1] - dom_g[2] - dom_f[1])..dom_f[2])

    if x in (dom_f[1] + dom_g[1])..(dom_f[2] + dom_g[1])
        new_dom = [dom_f[1] + dom_g[1], dom_f[2] + dom_g[1]]
        return legendre_same_length(f, g1, new_dom, new_dom, α_s=α_s, β_s=β_s)(x)

    elseif x in (dom_f[2] + dom_g[1])..(dom_f[1] + dom_g[2])
        new_dom = [dom_f[2] + dom_g[1], dom_f[1] + dom_g[2]]
        return legendre_same_length(f1, g2, new_dom, new_dom, α_s=α_s, β_s=β_s)(x) + 
               legendre_same_length(f, g3, new_dom, new_dom, α_s=α_s, β_s=β_s)(x)

    elseif x in (dom_f[1] + dom_g[2])..(dom_f[2] + dom_g[2])
        new_dom = [dom_f[1] + dom_g[2], dom_f[2] + dom_g[2]]
        return legendre_same_length(f, g4, new_dom, new_dom,  α_s=α_s, β_s=β_s)(x)

    else
        return 0.0

    end
end


"""
Calculates the convolution of two functions defined on arbitrary intervals as
in section 5.
"""
function legendre_conv(f, g, dom_f, dom_g; α_s=α_s, β_s=β_s)
    rat = (dom_g[2] - dom_g[1])/(dom_f[2] - dom_f[1])
    #@assert rat >= 1
    if rat < 1
        return legendre_conv(g, f, dom_g, dom_f, α_s=α_s, β_s=β_s)
    end
    r = Int(modf(rat)[2])
    if dom_f[2] - dom_f[1] == dom_g[2] - dom_g[1]
        return legendre_same_length(f, g, dom_f, dom_g, α_s=α_s, β_s=β_s)
    else
        if modf(rat)[1] == 0.0
        # d-c / b-a > 1 and integer
        # partition g into (d-c)/(b-a) subdomains and add
            funcs = Array{Function}(undef, r)
            for j in 1:r
                first = dom_g[1] + (j-1)*(dom_f[2] - dom_f[1])
                last = dom_g[1] + j*(dom_f[2]-dom_f[1])
                dⱼ = first..last
                gⱼ = x -> g(x)*(x in dⱼ)
                funcs[j] = gⱼ
            end
            h = x -> sum(legendre_same_length(f, funcs[j], dom_f, dom_g,
                                              α_s=α_s, β_s=β_s)(x) for j in 1:r)
            return h
        elseif 1 < rat < 2
        # d-c/b-a > 2
        # h is piecewise on 3 intervals
            return x -> h_12(f, g, dom_f, dom_g, x, α_s=α_s, β_s=β_s)
        else
        # split into sum satisfying conditions 1 and 2 
            i₁ = dom_g[1]..(dom_g[1] + (r-1)*(dom_f[2]-dom_f[1]))
            i₂ = (dom_g[1] + (r-1)*(dom_f[2]-dom_f[1]))..dom_g[2]
            g1 = x -> g(x)*(x in i₁)
            g2 = x -> g(x)*(x in i₂)
            h1 = x -> legendre_same_length(f, g1, dom_f, dom_g, α_s=α_s,
                                           β_s=β_s)(x)
            h2 = x -> h_12(f, g2, dom_f, dom_g, x, α_s=α_s, β_s=β_s)
            h = x -> h1(x) + h2(x)
            return h
        end  
    end
end


"""
Calculates f ⋆ g supported on dom_f and dom_g respectively, expressed as a
legendre series.
"""
function legendre_conv_series(f, g, dom_f, dom_g; α_s=α_s, β_s=β_s)
    h = legendre_conv(f, g, dom_f, dom_g, α_s=α_s, β_s=β_s)
    start = dom_f[1] + dom_g[1]
    stop = dom_f[2] + dom_g[2]
    legendreseries(h, interval=start..stop, N=α_s+β_s+2)
end

#left_conv_unit(sin, cos; α_s=100, β_s=100)


# n = 1000
# h= left_conv_unit(x->sin(x), x->cos(x), α_s = 3, β_s = 2)
# #θ = range(-1, 0, length=2n+1)[1:end-1]
# h[-1]

#≈ π*sin.(θ)
# T = legendre(-2..0)
# h = T[:,1:3]*[1,1,1]
# h[-1]

h =legendre_conv(sin, cos, [-1, 1], [-1, 1], α_s = 10, β_s = 11)

h.(range(-2, 2, length=100))[50]