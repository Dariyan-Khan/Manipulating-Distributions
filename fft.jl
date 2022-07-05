using LinearAlgebra
using Plots

f = x -> sin(3x) + sin(30x)

function dftmatrix(n)
    z = range(0, 2π, length=n+1)[1:end-1]
    [exp(-im*(k-1)*z[j]) for k=1:n, j=1:n]/sqrt(n)
end

function fftmatrix(n)
    σ = [1:2:2n-1;2:2:2n]
    P_σ = I(2n)[σ,:]
    Qₙ = dftmatrix(n)
    Dₙ = Diagonal([exp(im*k*π/n) for k=0:n-1])
    Q₂ₙ_adj = (1 / sqrt(2)) * (P_σ'*[Qₙ' Qₙ'; Qₙ'*Dₙ -Qₙ'*Dₙ])
    Q₂ₙ_adj'

end

n=31
z = range(0, 2π, length=2n+1)[1:end-1]
samples = f.(z)
A = fftmatrix(n)
f̂ = (1 / sqrt(2n)) * A * samples

# print(abs.(f̂))
# plot(0:length(f̂)-1, abs.(f̂))

function fourierconv(f, g, n)
    z = range(0, 2π, length=2n+1)[1:end-1]
    A = fftmatrix(n)
    f̂ = A*(f.(z))
    ĝ = A*(g.(z))
    ĉ = [f̂[k]*ĝ[k] for k in 1:2n]
    return ĉ
end

function conv(f, g, n)
    A = fftmatrix(n)
    ĉ = fourierconv(f, g, n)
    return A \ ĉ
end

plot(0:length(f̂)-1, abs.(f̂))




