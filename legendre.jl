using ClassicalOrthogonalPolynomials, Plots, SpecialFunctions

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
    c = legendrecoeff(f)
    g = T * c
    return g
end

# compute fourier transform of poly using leg series
function bleft(k, n, f, g)
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

function gammaleft(k, f, g)
    # find degree of f 
    # take sum and use bleft 
end