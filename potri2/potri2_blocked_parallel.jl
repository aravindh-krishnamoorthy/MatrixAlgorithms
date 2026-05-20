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
# Blocked parallel version
#
# Notes:
# - This refactor keeps the blocked algorithm from potri2_blocked_parallel.jl.
# - Anti-diagonal wavefront scheduling is applied on OUTER square blocks of size
#   bs_col.
# - bs_row is used only as an INNER row-panel size inside potri2_bo!, so the
#   off-diagonal block solve can be performed in smaller row chunks.
# - This implementation is for uplo = 'L'.
################################################################################

const BLAS = LinearAlgebra.BLAS

@inline _bstart(b::Integer, bs::Integer) = (b - 1) * bs + 1
@inline _bstop(b::Integer, bs::Integer, n::Integer) = min(b * bs, n)
@inline _brange(b::Integer, bs::Integer, n::Integer) = _bstart(b, bs):_bstop(b, bs, n)

@inline function potri2_bd!(
    X::StridedMatrix{T},
    bdiag::Integer,
    bs_col::Integer,
    W::StridedMatrix{T},
    bvec::StridedVector{T},
    rhs::StridedVector{T},
) where {T<:LinearAlgebra.BlasFloat}
    n  = size(X, 1)
    jr = _brange(bdiag, bs_col, n)
    jb = length(jr)
    j  = last(jr)

    Tinit = view(W, 1:jb, 1:jb)
    if j < n
        BLAS.gemm!('N', 'N', one(T), view(X, jr, j+1:n), view(X, j+1:n, jr), zero(T), Tinit)
    else
        fill!(Tinit, zero(T))
    end

    _potri2_inner!(view(X, jr, jr), Tinit, view(bvec, 1:jb), view(rhs, 1:jb))
    return X
end

@inline function _potri2_inner!(
    R::StridedMatrix{T},
    Tinit::StridedMatrix{T},
    b::StridedVector{T},
    rhs::StridedVector{T},
) where {T<:LinearAlgebra.BlasFloat}
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
    bs_col::Integer,
    bs_row::Integer,
    W::StridedMatrix{T},
) where {T<:LinearAlgebra.BlasFloat}
    n = size(X, 1)
    @assert brow > bcol

    # Current diagonal block columns (the RHS block columns)
    jr = _brange(brow, bs_col, n)
    jb = length(jr)
    j  = last(jr)

    # Earlier block rows / columns for this off-diagonal block
    cr = _brange(bcol, bs_col, n)
    cb = length(cr)

    # Workspace B corresponds to the original whole-column solve restricted to
    # rows cr and columns jr:
    #
    #   U(cr,cr) * B = -( X(cr,j+1:n)X(j+1:n,jr) + X(cr,jr)X(jr,jr)
    #                    + sum_{p=bcol+1}^{brow-1} X(cr,pr) * B_pr )
    #
    # The final output block is X(jr,cr) = B'.
    B = view(W, 1:cb, 1:jb)

    # Trailing update contribution: X(cr, j+1:n) * X(j+1:n, jr)
    if j < n
        BLAS.gemm!('N', 'N', one(T), view(X, cr, j+1:n), view(X, j+1:n, jr), zero(T), B)
    else
        fill!(B, zero(T))
    end

    # Local diagonal-block contribution: X(cr, jr) * X(jr, jr)
    BLAS.gemm!('N', 'N', one(T), view(X, cr, jr), view(X, jr, jr), one(T), B)

    # Contributions from already-computed blocks to the right in the same output row:
    # B_pr is stored as X(jr, pr)', so use gemm with transposed second factor.
    for bp = bcol+1:brow-1
        pr = _brange(bp, bs_col, n)
        BLAS.gemm!('N', 'C', one(T), view(X, cr, pr), view(X, jr, pr), one(T), B)
    end

    # Now solve U(cr,cr) * B = -B using INNER row panels of height bs_row.
    # This is just blocked back-substitution inside the current outer block.
    for pend = cb:-bs_row:1
        pstart = max(1, pend - bs_row + 1)
        ploc   = pstart:pend
        pr     = cr[ploc]
        P      = view(B, ploc, :)

        # Coupling from panels to the right inside the same outer block.
        if pend < cb
            qloc = (pend + 1):cb
            qr   = cr[qloc]
            BLAS.gemm!('N', 'N', one(T), view(X, pr, qr), view(B, qloc, :), one(T), P)
        end

        # Diagonal panel solve.
        BLAS.trsm!('L', 'U', 'N', 'N', -one(T), view(X, pr, pr), P)
    end

    @views X[jr, cr] .= B'
    return X
end

function potri2_blocked_parallel!(
    uplo::Char,
    X::StridedMatrix{T};
    bs_col::Int = 64,
    bs_row::Int = bs_col,
) where {T<:LinearAlgebra.BlasFloat}
    n = size(X, 1)
    @assert size(X, 2) == n
    @assert bs_col > 0
    @assert bs_row > 0

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

    nblk = cld(n, bs_col)
    nt   = Threads.maxthreadid()
    Wwrk = [Matrix{T}(undef, bs_col, bs_col) for _ = 1:nt]
    bwrk = [Vector{T}(undef, bs_col) for _ = 1:nt]
    xwrk = [Vector{T}(undef, bs_col) for _ = 1:nt]

    # Anti-diagonal wavefront on OUTER square blocks of size bs_col.
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
                potri2_bd!(X, brow, bs_col, Wwrk[tid], bwrk[tid], xwrk[tid])
            else
                potri2_bo!(X, brow, bcol, bs_col, bs_row, Wwrk[tid])
            end
        end
    end

    X .= Hermitian(X, :L)
    return X
end
