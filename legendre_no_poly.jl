using ClassicalOrthogonalPolynomials, Plots


function legendreft(m) 
    f = x -> 2 * (-1im)^m * sphericalbesselj(m, x)
    return f
end

# compute legendre series of a poly
function legendrecoeff(f; N=nothing)
    T = Legendre()
    x = axes(T, 1)
    N == nothing ? c = T \ f.(x) : T[:,1:N] \ f.(x)
end

function legendreseries(f; N=nothing)
    T = Legendre()
    c = legendrecoeff(f, N=N)
    N == nothing ? g = T * c : T[:, 1:N] * c
    return g
end

function bleft(k::Int, n::Int, f, g; N=100)
    α = legendrecoeff(f, N=N)
    β = legendrecoeff(g, N=N)

    if k > n
        if n == 0
            if k == 0
                return α[1] - α[2]/3
            else
                return α[k-1]/(2k-1) - α[k+1]/(2k+3)
            end
        elseif n == 1
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
