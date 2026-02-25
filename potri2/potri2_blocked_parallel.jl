################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
# Contains reviewed and tested contributions from ChatGPT 5.2 Thinking
################################################################################

################################################################################
# Invert a positive definite matrix using the algorithm in
# "Matrix Inversion Using Cholesky Decomposition", by Aravindh Krishnamoorthy
#   and Deepak Menon, arXiv:1111.4144.
################################################################################
# Blocked version
################################################################################
const BLAS = LinearAlgebra.BLAS
function potri2_blocked_parallel!(uplo::Char, X::StridedMatrix{T}; bs::Int=64) where {T<:LinearAlgebra.BlasFloat}
    n = size(X, 1)
    @assert size(X, 2) == n
    @assert bs > 0

    if uplo == 'U'
        @inbounds for i = 1:n
            di = X[i, i]
            for k = i:n
                X[i, k] = di * X[i, k]
            end
        end
    else
        @inbounds for i = 1:n
            di = X[i, i]
            for k = i:n
                X[i, k] = di * conj(X[k, i])
            end
        end
    end

    Bblk    = Matrix{T}(undef, n, bs)
    bvec    = Vector{T}(undef, bs)
    rhs     = Vector{T}(undef, bs)

    @views for j = n:-bs:1
        jb = min(j, bs)
        jr1 = j -jb + 1
        jr  = jr1:j
        j0  = j - jb

        if j < n
            BLAS.gemm!('N', 'N', one(T), view(X, jr, j+1:n), view(X, j+1:n, jr), zero(T), view(Bblk, jr, :))
        else
            fill!(Bblk, zero(T))
        end
        _potri2_inner!(view(X, jr, jr), view(Bblk, jr, :), view(bvec, 1:jb), view(rhs, 1:jb))
        if j0 > 0
            B = view(Bblk, 1:j0, 1:jb)
            BLAS.gemm!('N', 'N', one(T), view(X, 1:j0, j+1:n), view(X, j+1:n, jr), zero(T), B)
            BLAS.gemm!('N', 'N', one(T), view(X, 1:j0, jr), view(X, jr, jr), one(T), B)
            BLAS.trsm!('L', 'U', 'N', 'N', -one(T), view(X, 1:j0, 1:j0), B)
            @views X[jr1:jr1+jb-1, 1:j0] .= B[1:j0, 1:jb]'
        end
    end
    X .= Hermitian(X, :L)
    return X
end

@inline function potri2_bd!(X::AbstractMatrix{T}, j::Integer, jb::Integer) where {T}
end

@inline function potri2_bo!(uplo::Char, X::AbstractMatrix{T}, i::Integer, ib::Integer, j::Integer, jb::Integer) where {T}
end

@inline function _potri2_inner!(R::StridedMatrix{T}, Tinit::StridedMatrix{T},
                                         b::StridedVector{T}, rhs::StridedVector{T}) where {T<:LinearAlgebra.BlasFloat}
    n = size(R, 1)
    @inbounds for i = 1:n
        b[i] = Tinit[i, n]
    end

    @views for j = n:-1:1
        x = view(rhs, 1:j)
        copyto!(x, view(b, 1:j))
        @inbounds for i = 1:j
            x[i] = -x[i]
        end
        @inbounds x[j] += one(T)

        A = view(R, 1:j, 1:j)
        BLAS.trsv!('U', 'N', 'N', A, x)

        @inbounds for i = 1:j
            R[j, i] = conj(x[i])
        end

        if j > 1
            A2 = view(R, 1:j-1, j:n)
            v  = view(R, j:n, j-1)
            y  = view(b, 1:j-1)
            BLAS.gemv!('N', one(T), A2, v, zero(T), y)
            @inbounds for i = 1:j-1
                y[i] += Tinit[i, j-1]
            end
        end
    end

    @inbounds for i = 1:n
        for j = i+1:n
            R[i, j] = conj(R[j, i])
        end
    end

    return R
end
