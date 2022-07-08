using ClassicalOrthogonalPolynomials, Plots, SpecialFunctions, Polynomials, 
      Parameter

# compute fourier transform of m degree Legendre poly
function legendreft(m) 
    f = x -> 2 * (-1im)^m * sphericalbesselj(m, x)
    return f
end

# compute legendre series of a poly
function legendrecoeff(f)
    T = Legendre()
    x = axes(T, 1)
    c = T \ f.(x)
    return c
end

function legendreseries(f)
    T = Legendre()
    c = legendrecoeff(f)
    g = T * c
    return g
end

# compute fourier transform of poly using leg series
#legendreseries(x -> x^2)

struct Poly
    domain::Vector
    p::Polynomial
end
degree!(pol::Poly) = degree(pol.p)
l_coeffs(pol::Poly) = legendrecoeff(pol.p)
l_series(pol::Poly) =  legendreseries(pol.p)

function eval!(pol::Poly, x::Real)
    if x>=pol.domain[1] && x<=pol.domain[2]
        return pol.p(x)
    else 
        println("Outside domain")
        return 0
    end
end 

function xshift(f::Poly, c::Real)
    #given f(x) finds the polynomial f(x+c)
    N = degree!(f) + 1
    xs = range(start=f.domain[1], stop=f.domain[2], length=N)
    ys = (f.p).(xs)
    xs = xs .- c
    return Poly(f.domain .- c, fit(xs, ys))
end

function xreflect(f::Poly)
    #given f(x) finds the polynomial f(x+c)
    N = degree!(f) + 1
    xs = range(start=f.domain[1], stop=f.domain[2], length=N)
    ys = (f.p).(xs)
    xs = xs * -1
    return Poly(f.domain * -1, fit(xs, ys))
end


function poly_conv(f::Poly, g::Poly, c::Real)
    poly_conv_inner(f, g, c, f.domain, g.domain)
end


function poly_conv_inner(f::Poly, g::Poly, c::Real, d₁::Vector, d₂::Vector)
    @assert (c >= d₁[1] + d₂[1]) && (c <= d₁[2] + d₂[2])
    a = max(d₁[1],c-d₂[2])
    b = min(d₁[2], c-d₂[1])
    @assert a <= b
    g = xreflect(xshift(g, c))
    return integrate(f.p*g.p, a, b)
end

function piecewise_conv(f::Poly, g::Poly)
   skip = false
   arr = Poly.([[0, 1]]*3, zeros(3))
   #pol = [Poly([0,0], 0), Poly([0,0], 0), Poly([0,0], 0)]
   a, b, c, d = f.domain[1], f.domain[2], g.domain[1], g.domain[2]
   domains = [[a+c, b+c], [b+c, a+d], [a+d, b+d]]
   for i in 1:3
    dom = domains[i]
    if dom[1] != dom[2]
        deg = degree!(f) + degree!(g) + mod(i,2)
        xs = range(start=dom[1], stop=dom[2], length=deg+1)
        lam = x -> poly_conv(f, g, x)
        ys = lam.(xs)
        println(fit(xs, ys))
        arr[i] = Poly(dom, fit(xs, ys))
        # append!(arr, Poly(dom, fit(xs, ys)))
    else
        skip = true
    end
    #pol[i] = Poly(dom, fit(xs, ys))
   end
   if skip
    arr = [arr[1], arr[3]]
   end
   return arr
end

# FOR LEGENDRE SERIES ON [1,1], then extend to on [a,b] [c,d] where d-c=b-a 

# find k,n-th entry of B left
function bleft(k::Int, n::Int, f::Poly, g::Poly)
    α = legendrecoeff(f)
    β = legendrecoeff(g)
    if k > n
        if n == 0
            if k == 0
                return α[1] - α[2]/3
            else
                return α[k-1]/(2k-1) - α[k+1]/(2k+3)
            end
        if n == 1
            if k == 0
                return -bleft(1, 0, f, g)/3
            else
                return bleft(k-1, 0, f, g)/(2k-1) - bleft(k,0,f,g) - bleft(k+1,0,f,g)/(2k+3)
            end
        else
            return -(2n-1)/(2k+3) * bleft(k+1,n-1,f,g) + (2n-1)/(2k-1) * bleft(k-1,n-1,f,g) + bleft(k, n-2,f,g)
        end
    else
        return bleft(n,k,f,g) * (-1)^(n+k) * (2k+1)/(2n+1)
    end
end
end


# find k,n-th entry of B right
function bright(k::Int, n::Int, f::Poly, g::Poly)
    α = legendrecoeff(f)
    β = legendrecoeff(g)
    if k > n
        if n == 0
            if k == 0
                return α[1] + α[2]/3
            else
                return α[k-1]/(2k-1) + α[k+1]/(2k+3)
            end
        if n == 1
            if k == 0
                return -bright(1, 0, f, g)/3
            else
                return bright(k-1, 0, f, g)/(2k-1) + bright(k,0,f,g) - bright(k+1,0,f,g)/(2k+3)
            end
        else
            return -(2n-1)/(2k+3) * bright(k+1,n-1,f,g) + (2n-1)/(2k-1) * bright(k-1,n-1,f,g) + bright(k, n-2,f,g)
        end
    else
        return bright(n,k,f,g) * (-1)^(n+k) * (2k+1)/(2n+1)
    end
end
end

# find γₖ left  
function gammaleft(k::Int, f::Poly, g::Poly)
    # find degree of f 
    # take sum and use bleft 
    N = degree!(f)
    β = legendrecoeff(g)
    ret = 0 
    for i in 0:N
        ret += β[i+1] * bleft(k, i, f, g)
    end
    return ret
end


# find γₖ right
function gammaright(k::Int, f::Poly, g::Poly)
    N = degree!(f)
    β = legendrecoeff(g)
    ret = 0 
    for i in 0:N
        ret += β[i+1] * bright(k, i, f, g)
    end
    return ret
end


function legendreconv_1_minus_1(f::Poly, g::Poly)
    T = Legendre()
    lam_left = k -> gammaleft(k, f, g)
    γ_left = lam_left.(0:degree!(f) + degree!(g) + 1)
    #Calculate h but then we need to map x -> x+1
    pre_h_left = T * γ_left
    h_left = x ->  pre_h_left(x+1)
    h_left = T / T \ h_left

    lam_right = k -> gammaright(k, f, g)
    γ_right = lam_right.(0:degree!(f) + degree!(g) + 1)
    #Calculate h but then we need to map x -> x+1
    pre_h_right = T * γ_right
    h_right = x ->  pre_h_right(x-1)
    h_right = T / T \ h_right

end

