module ManipulatingDistributions
using LinearAlgebra, ClassicalOrthogonalPolynomials, Distributions
export fft_conv, legendre_conv


include("fft.jl")
include("legendre.jl")

end # module
