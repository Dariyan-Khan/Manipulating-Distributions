using ClassicalOrthogonalPolynomials, Plots, SpecialFunctions

# compute fourier transform of m degree Legendre poly
function LegFT(m) 
    f = x -> 2 * (-1im)^m * sphericalbesselj(m, x)
    return f
end

# compute legendre series of a poly

# compute fourier transform of poly using leg series