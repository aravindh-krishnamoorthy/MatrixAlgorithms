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
function _invtri_unblocked_lower!(A::AbstractMatrix{T}) where {T}
    n = size(A,1)
    dinv = similar(A, n)
    @inbounds for i=1:n
        dinv[i] = inv(A[i,i])
        A[i,i] = dinv[i]
    end
    @inbounds for j=1:n
        for i=j+1:n
            s = zero(T)
            @simd for k=j:i-1
                s += A[i,k]*A[k,j]
            end
            A[i,j] = -dinv[i]*s
        end
    end
    A
end

function _invtri_unblocked_upper!(A::AbstractMatrix{T}) where {T}
    n = size(A,1)
    dinv = similar(A, n)
    @inbounds for i=1:n
        dinv[i] = inv(A[i,i])
        A[i,i] = dinv[i]
    end
    @inbounds for j=n:-1:1
        for i=j-1:-1:1
            s = zero(T)
            @simd for k=i+1:j
                s += A[i,k]*A[k,j]
            end
            A[i,j] = -dinv[i]*s
        end
    end
    A
end

function _invtri_rec_lower!(A::AbstractMatrix{T}, bs::Int) where {T}
    n = size(A,1)
    n <= bs && return _invtri_unblocked_lower!(A)
    n1 = n >>> 1
    n2 = n - n1
    A11 = @view A[1:n1, 1:n1]
    A21 = @view A[n1+1:n, 1:n1]
    A22 = @view A[n1+1:n, n1+1:n]
    _invtri_rec_lower!(A22, bs)
    _invtri_rec_lower!(A11, bs)
    T1 = similar(A21)
    mul!(T1, A21, A11)
    mul!(A21, A22, T1)
    @inbounds @simd for idx in eachindex(A21)
        A21[idx] = -A21[idx]
    end
    A
end

function _invtri_rec_upper!(A::AbstractMatrix{T}, bs::Int) where {T}
    n = size(A,1)
    n <= bs && return _invtri_unblocked_upper!(A)
    n1 = n >>> 1
    n2 = n - n1
    A11 = @view A[1:n1, 1:n1]
    A12 = @view A[1:n1, n1+1:n]
    A22 = @view A[n1+1:n, n1+1:n]
    _invtri_rec_upper!(A11, bs)
    _invtri_rec_upper!(A22, bs)
    T1 = similar(A12)
    mul!(T1, A11, A12)
    mul!(A12, T1, A22)
    @inbounds @simd for idx in eachindex(A12)
        A12[idx] = -A12[idx]
    end
    A
end

function potri2_blocked!(uplo::Char, X::AbstractMatrix{T}; bs::Int=64) where {T}
    n = size(X,1)
    if uplo == 'L'
        _invtri_rec_lower!(X, bs)
        Y = similar(X)
        mul!(Y, adjoint(X), X)
        copyto!(X, Y)
    else
        _invtri_rec_upper!(X, bs)
        Y = similar(X)
        mul!(Y, X, adjoint(X))
        copyto!(X, Y)
    end
    X
end
