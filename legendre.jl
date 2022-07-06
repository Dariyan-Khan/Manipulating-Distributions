using ClassicalOrthogonalPolynomials, Plots, SpecialFunctions

# compute fourier transform of m degree Legendre poly
function legendreft(m) 
    f = x -> 2 * (-1im)^m * sphericalbesselj(m, x)
    return f
end

# compute legendre series of a poly
function legendreseries(f)
    T = Legendre()
    x = axes(T, 1)
    c = T \ f.(x)
    return c
end
# compute fourier transform of poly using leg series