using ClassicalOrthogonalPolynomials, FastTransforms, Test, QuasiArrays

P = Legendre()
T = ChebyshevT()

x = axes(P,1)
P \ exp.(x)


# finite number of coefficients
n = 10
Pₙ = P[:,1:n]
Tₙ = T[:,1:n]
f = Pₙ * (Pₙ \ exp.(x))

x = ChebyshevGrid{1}(100)
f[x] # uses Clenshaw to evaluate

Pₙ \ f # Legendre coefficients
@test leg2cheb(Pₙ \ f) ≈ Tₙ \ f # Chebyshev T coefficients

@test ichebyshevtransform([leg2cheb(Pₙ \ f); zeros(90)]) ≈ f[x]

# ichebyshevtransform + leg2cheb is roughly O(n*log(n)) for evaluating at n points
# Clenshaw is roughly O(n^2)

# different interval

P₀₁ = legendre(0..1)
T₀₁ = chebyshevt(0..1)
x = axes(P₀₁,1)

f = P₀₁ * (P₀₁ \ exp.(x))

f[0.1] ≈ exp(0.1)

c = 5
fₛ = legendre(c..(1+c)) * (P₀₁ \ f)
@test fₛ[0.1+c] ≈ f[0.1]



# indefinite integration

m = 5; n = 3




function Pconv(x)
    T = chebyshevt(-1 .. (x+1))
    t = axes(T,1)
    Pₘ = P[:,m+1]
    Pₙ = P[:,n+1]
    f = t -> (Pₘ[t] * Pₙ[x-t])::Float64
    g = T / T \ f.(t) # represents Pₘ(t)*Pₙ(x-t) for t in -1 .. (x_1)
    sum(g)::Float64 # Integral
end


P₋₂ = legendre(-2..0)
x₋₂ = axes(P₋₂,1)

P₋₂ \ Pconv.(x₋₂)


## Piecewise ClassicalOrthogonalPolynomials

using PiecewiseOrthogonalPolynomials
# Discontinuous
C₀ = ContinuousPolynomial{0}([0, 1, 2]) # Legendre on 0..1 and 1..2
# Continuous
C₁ = ContinuousPolynomial{1}([0, 1, 2]) # (1-x^2) * Jacobi(1,1) on 0..1 and 1..2, plus piecewise affine

using Plots

x = range(0, 2; length=100)
plot(x, C₀[x,1:6])

plot(x, C₁[x,1:7])