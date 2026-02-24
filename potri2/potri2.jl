################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
# Contains contributions by ChatGPT 5.2 Thinking
################################################################################

using LinearAlgebra
using Libdl

include("potri2_blocked.jl")

################################################################################
# Invert a positive definite matrix using the algorithm in
# "Matrix Inversion Using Cholesky Decomposition", by Aravindh Krishnamoorthy
#   and Deepak Menon, arXiv:1111.4144.
################################################################################
# Reference version
################################################################################
function potri2!(uplo::Char, X::AbstractMatrix{T}) where {T}
    n = size(X,1)
    d = zeros(T, n)
    z = zero(T)
    @inbounds if uplo == 'U'
        for j = 1:n
            Xjj = inv(X[j,j])
            X[j,j] = Xjj
            d[j] = Xjj
            fill!(@view(X[j+1:n, j]), z)
        end
        for j = n:-1:1
            for k = n:-1:j+1
                xkj = X[k,j]
                @simd for i = 1:j
                    X[j,i] = X[j,i] - X[i,k]*xkj
                end
            end
            for k = j:-1:1
                tmp = X[j,k]*d[k]
                X[j,k] = conj(tmp)
                @simd for i = 1:k-1
                    X[j,i] = X[j,i] - X[i,k]*tmp
                end
            end
        end
        for i = 1:n, j = i+1:n
            X[i,j] = X[j,i]'
        end
        return X
    else
        for j = 1:n
            Xjj = inv(X[j,j])
            X[j,j] = Xjj
            d[j] = Xjj
            fill!(@view(X[1:j-1, j]), z)
        end
        for j = n:-1:1
            for k = n:-1:j+1
                xjk = X[j,k]
                @simd for i = 1:j
                    X[i,j] = X[i,j] - X[k,i]*xjk
                end
            end
            for k = j:-1:1
                tmp = X[k,j]*d[k]
                X[k,j] = conj(tmp)
                @simd for i = 1:k-1
                    X[i,j] = X[i,j] - X[k,i]*tmp
                end
            end
        end
        for i = 1:n, j = 1:i-1
            X[i,j] = X[j,i]'
        end
        return X
    end
end
