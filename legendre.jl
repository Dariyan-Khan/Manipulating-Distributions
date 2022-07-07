using ClassicalOrthogonalPolynomials, Plots, SpecialFunctions, Polynomials, 
      Parameters

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

function shift(f::Poly, c::Real)
    #given f(x) finds the polynomial f(x+c)
    N = degree(f) + 1
    xs = range(start=f.domain[1], stop=f.domain[2], length=N)
    ys = (f.p).(xs)
    xs = xs .- c
    return fit(xs, ys)
end




function poly_conv(f::Poly, g::Poly, x::Real, d₁::Vector, d₂::Vector)
    @assert length(d₁) == 2 && d₁[1] <= d₁[2]
    @assert length(d₂) == 2 && d₂[1] <= d₂[2]
    @assert x >= d₁[1] + d₂[2] && x <= d₁[2] + d₂[2]
    a = max(d₁[1],x-d₂[2])
    b = min(d₁[2], x-d₂[1])
    
    return 
end

function conv(f::Poly, g::Poly)
    
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

end

# find γₖ left  
function gammaleft(k::Int, f::Poly, g::Poly)
    # find degree of f 
    # take sum and use bleft 
    N = degree(f)
    β = legendrecoeff(g)
    ret = 0 
    for i in 0:N
        ret += β[i+1] * bleft(k, i, f, g)
    end
    return ret
end

# find γₖ right
function gammaright(k::Int, f::Poly, g::Poly)

end