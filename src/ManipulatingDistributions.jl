module ManipulatingDistributions
using LinearAlgebra, ClassicalOrthogonalPolynomials
export fft_conv, legendre_conv


include("fft.jl")
include("legendre.jl")

end # module
