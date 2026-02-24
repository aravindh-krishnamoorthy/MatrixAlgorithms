################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
# Contains contributions by ChatGPT 5.2 Thinking
################################################################################

################################################################################
# Invert a positive definite matrix using the algorithm in
# "Matrix Inversion Using Cholesky Decomposition", by Aravindh Krishnamoorthy
#   and Deepak Menon, arXiv:1111.4144.
################################################################################
# Blocked version
################################################################################
const BLAS = LinearAlgebra.BLAS
function potri2_blocked!(uplo::Char, X::StridedMatrix{T}; bs::Int=64) where {T<:LinearAlgebra.BlasFloat}
    n = size(X, 1)
    @assert size(X, 2) == n

    @inbounds for i = 1:n
        X[i, :] = X[i, i]*X[i, :]
    end

    b   = zeros(T, n)
    rhs = similar(b)

    @views for j = n:-1:1
        x = view(rhs, 1:j)
        copyto!(x, view(b, 1:j))
        x[1:j] = -x[1:j]
        x[j] += one(T)

        A = view(X, 1:j, 1:j)
        BLAS.trsv!('U', 'N', 'N', A, x)

        copyto!(view(X, j, 1:j), conj(x))

        if j > 1
            A2 = view(X, 1:j-1, j:n)
            v  = view(X, j:n, j-1)
            y  = view(b, 1:j-1)
            V  = reshape(v, :, 1)
            Y  = reshape(y, :, 1)
            BLAS.gemm!('N', 'N', one(T), A2, V, zero(T), Y)
        end
    end

    @inbounds for i = 1:n
        for j = i+1:n
            X[i, j] = conj(X[j, i])
        end
    end

    return X
end
