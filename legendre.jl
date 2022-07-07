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




a = Poly([1, 2], Polynomial([1, 2]))

















