module MatrixAlgorithms

using LinearAlgebra

########################################
# potri:
# Positive definite matrix inverse
########################################
include("../potri2/potri2.jl")
########################################
# trsr:
# Matrix square root
########################################
include("../trsr/trsr.jl")
########################################
# matrixfunctions:
# Functions of matrices
################################################################################
# PACKAGE-SPECIFIC EXCEPTIONS
################################################################################
# The following exception is thrown by functions that fail when the input matrix
# contains repeated eigenvalues.
################################################################################
struct RepeatedEigenvalueException <: Exception
end
################################################################################
# ALGORITHMS
################################################################################
## Schur-Parlett algorithms
################################################################################
include("../matrixfunctions/schur_parlett.jl")

end
