################################################################################
# This file is a part of the Julia package: MatrixFunctions.jl
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
#
################################################################################

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
include("schur_parlett.jl")
