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
# Blocked parallel version refactored from potri2_blocked_parallel.jl.
# The underlying blocked algorithm is unchanged. The diagonal and off-diagonal
# computations are factored into block-local kernels, and the block traversal is
# scheduled by anti-diagonal wavefronts.
################################################################################
const BLAS = LinearAlgebra.BLAS

@inline _bstart(b::Integer, bs::Integer) = (b - 1) * bs + 1
@inline _bstop(b::Integer, bs::Integer, n::Integer) = min(b * bs, n)
@inline _brange(b::Integer, bs::Integer, n::Integer) = _bstart(b, bs):_bstop(b, bs, n)
@inline _bsize(b::Integer, bs::Integer, n::Integer) = length(_brange(b, bs, n))

@inline function potri2_bd!(
    X::StridedMatrix{T},
    br::Integer,
    bs::Integer,
    Bblk::StridedMatrix{T},
    bvec::StridedVector{T},
    rhs::StridedVector{T},
) where {T<:LinearAlgebra.BlasFloat}
    n  = size(X, 1)
    rr = _brange(br, bs, n)
    jb = length(rr)
    j  = last(rr)

    Tinit = view(Bblk, 1:jb, 1:jb)
    if j < n
        BLAS.gemm!('N', 'N', one(T), view(X, rr, j+1:n), view(X, j+1:n, rr), zero(T), Tinit)
    else
        fill!(Tinit, zero(T))
    end

    _potri2_inner!(view(X, rr, rr), Tinit, view(bvec, 1:jb), view(rhs, 1:jb))
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

@inline function potri2_bo!(
    X::StridedMatrix{T},
    brow::Integer,
    bcol::Integer,
    bs::Integer,
    W::StridedMatrix{T},
) where {T<:LinearAlgebra.BlasFloat}
    n = size(X, 1)
    @assert brow > bcol

    rr = _brange(brow, bs, n)
    cr = _brange(bcol, bs, n)
    jb = length(rr)
    ib = length(cr)
    j  = last(rr)

    B = view(W, 1:ib, 1:jb)

    # Trailing update from blocks to the right of the current block-row.
    if j < n
        BLAS.gemm!('N', 'N', one(T), view(X, cr, j+1:n), view(X, j+1:n, rr), zero(T), B)
    else
        fill!(B, zero(T))
    end

    # Local contribution from the current diagonal block.
    BLAS.gemm!('N', 'N', one(T), view(X, cr, rr), view(X, rr, rr), one(T), B)

    # Block backward-substitution contribution from already computed blocks
    # on the same block-row: sum_{p=bcol+1}^{brow-1} X[cr, pr] * X[rr, pr]'.
    for bp = bcol+1:brow-1
        pr = _brange(bp, bs, n)
        BLAS.gemm!('N', 'C', one(T), view(X, cr, pr), view(X, rr, pr), one(T), B)
    end

    # Solve with the local leading diagonal block from the still-unprocessed
    # triangular factor: B <- -X[cr, cr]^{-1} * B.
    BLAS.trsm!('L', 'U', 'N', 'N', -one(T), view(X, cr, cr), B)

    # Store only the computed output block in the lower triangle.
    @views X[rr, cr] .= B'
    return X
end

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

    nblk = cld(n, bs)
    nt   = Threads.maxthreadid()
    Bwrk = [Matrix{T}(undef, bs, bs) for _ = 1:nt]
    bwrk = [Vector{T}(undef, bs) for _ = 1:nt]
    xwrk = [Vector{T}(undef, bs) for _ = 1:nt]

    # Anti-diagonal wavefront on the lower-triangular block matrix.
    # Wave index is brow + bcol, traversed from bottom-right to top-left.
    for wave = 2nblk:-1:2
        clo = max(1, wave - nblk)
        chi = min(nblk, fld(wave, 2))

        Threads.@threads for bcol = clo:chi
            brow = wave - bcol
            if brow < bcol
                continue
            end

            tid = Threads.threadid()
            if brow == bcol
                potri2_bd!(X, brow, bs, Bwrk[tid], bwrk[tid], xwrk[tid])
            else
                potri2_bo!(X, brow, bcol, bs, Bwrk[tid])
            end
        end
    end

    X .= Hermitian(X, :L)
    return X
end
