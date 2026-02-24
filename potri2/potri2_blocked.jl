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
function potri2_blocked!(uplo::Char, X::StridedMatrix{T}; bs::Int=64) where {T<:LinearAlgebra.BlasFloat}
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
                X[i, k] = di * X[k, i]
            end
        end
    end

    Tb_full = Matrix{T}(undef, bs, bs)
    Rtmp    = Matrix{T}(undef, bs, bs)
    Bblk    = Matrix{T}(undef, n, bs)
    bvec    = Vector{T}(undef, bs)
    rhs     = Vector{T}(undef, bs)

    @views for j = n:-bs:1
        nb = min(j, bs)
        jr1 = j - nb + 1
        jr  = jr1:j
        j0  = j - nb

        Tb = view(Tb_full, 1:nb, 1:nb)

        if j < n
            A = view(X, jr, j+1:n)
            C = view(X, j+1:n, jr)
            BLAS.gemm!('N', 'N', one(T), A, C, zero(T), Tb)
        else
            fill!(Tb, zero(T))
        end

        copyto!(view(Rtmp, 1:nb, 1:nb), view(X, jr, jr))
        _potri2_inner!(view(Rtmp, 1:nb, 1:nb), Tb, view(bvec, 1:nb), view(rhs, 1:nb))
        copyto!(Tb, view(Rtmp, 1:nb, 1:nb))

        @inbounds for ii = 1:nb
            r = jr1 + ii - 1
            for jj = 1:ii
                c = jr1 + jj - 1
                X[r, c] = Tb[ii, jj]
            end
        end

        if j0 > 0
            B = view(Bblk, 1:j0, 1:nb)

            A1 = view(X, 1:j0, jr)
            BLAS.gemm!('N', 'N', one(T), A1, Tb, zero(T), B)

            if j < n
                A2 = view(X, 1:j0, j+1:n)
                C2 = view(X, j+1:n, jr)
                BLAS.gemm!('N', 'N', one(T), A2, C2, one(T), B)
            end

            @inbounds for col = 1:nb, row = 1:j0
                B[row, col] = -B[row, col]
            end

            U = view(X, 1:j0, 1:j0)
            BLAS.trsm!('L', 'U', 'N', 'N', one(T), U, B)

            @inbounds for col = 1:nb
                rrow = jr1 + col - 1
                for row = 1:j0
                    X[rrow, row] = conj(B[row, col])
                end
            end
        end
    end

    @inbounds for i = 1:n
        for j = i+1:n
            X[i, j] = conj(X[j, i])
        end
    end

    return X
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
