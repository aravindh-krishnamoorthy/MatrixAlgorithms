# This exception is thrown by functions that fail when the input matrix
# contains repeated eigenvalues.
struct RepeatedEigenvalueException <: Exception
end

# Schur-Parett algorithms
include("schur_parlett.jl")
