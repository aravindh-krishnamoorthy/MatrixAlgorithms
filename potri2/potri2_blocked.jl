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
@inline function _potri2_inner!(R::StridedMatrix{T}, Tinit::StridedMatrix{T},
                                         b::StridedVector{T}, rhs::StridedVector{T}) where {T<:Number}
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

        ldiv!(UpperTriangular(view(R, 1:j, 1:j)), x)

        @inbounds for i = 1:j
            R[j, i] = conj(x[i])
        end

        if j > 1
            y  = view(b, 1:j-1)
            mul!(y, view(R, 1:j-1, j:n), view(R, j:n, j-1))
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

function potri2_blocked!(uplo::Char, X::StridedMatrix{T}; bs::Int=64) where {T<:Number}
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

    if n <= bs
        Tinit = zeros(T, n, n)
        bvec  = Vector{T}(undef, n)
        rhs   = Vector{T}(undef, n)
        return _potri2_inner!(X, Tinit, bvec, rhs)
    end

    Bblk    = Matrix{T}(undef, n, bs)
    bvec    = Vector{T}(undef, bs)
    rhs     = Vector{T}(undef, bs)

    @views for j = n:-bs:1
        nb = min(j, bs)
        jr1 = j - nb + 1
        jr  = jr1:j
        j0  = j - nb

        if j < n
            mul!(view(Bblk, 1:j, 1:nb), view(X, 1:j, j+1:n), view(X, j+1:n, jr), one(T), zero(T))
        else
            fill!(Bblk, zero(T))
        end
        _potri2_inner!(view(X, jr, jr), view(Bblk, jr, :), view(bvec, 1:nb), view(rhs, 1:nb))
        if j0 > 0
            B = view(Bblk, 1:j0, 1:nb)
            mul!(B, view(X, 1:j0, jr), view(X, jr, jr), one(T), one(T))
            ldiv!(UpperTriangular(view(X, 1:j0, 1:j0)), rmul!(B, -one(T)))
            @views X[jr1:jr1+nb-1, 1:j0] .= B[1:j0, 1:nb]'
        end
    end
    X .= Hermitian(X, :L)
    return X
end
