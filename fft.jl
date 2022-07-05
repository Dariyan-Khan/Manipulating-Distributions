using LinearAlgebra

f = x -> sin(x) + sin(3*x)

function dftmatrix(n)
    z = range(0, 2π, length=n+1)[1:end-1]
    [exp(-im*(k-1)*z[j]) for k=1:n, j=1:n]/sqrt(n)
end

function fftmatrix(n)
    σ = [1:2:2n-1;2:2:2n]
    P_σ = I(2n)[σ,:]
    Qₙ = dftmatrix(n)
    Dₙ = Diagonal([exp(im*k*π/n) for k=0:n-1])
    (P_σ'*[Qₙ' Qₙ'; Qₙ'*Dₙ -Qₙ'*Dₙ])
end

z = range(0, 2π, length=n+1)[1:end-1]
f.(z)

A = fftmatrix(n)




