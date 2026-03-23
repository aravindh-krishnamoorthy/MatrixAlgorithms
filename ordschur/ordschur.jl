
################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
################################################################################

using LinearAlgebra

################################################################################
#
# Reordering of eigenvalues of Schur decomposition. Algorithm based on:
#   Daniel Kressner, "Block algorithms for reordering standard and generalized Schur forms,"
#   ACM Transactions on Mathematical Software Volume 32 Issue 4 pp 521–532 https://doi.org/10.1145/1186785.1186787
#
################################################################################
function ordschur(S::LinearAlgebra.Schur, t::AbstractVector{<:Integer})
    T = copy(S.T)
    Z = copy(S.Z)

    n = size(T, 1)
    n == size(T, 2) || throw(ArgumentError("S.T must be square"))
    size(Z, 1) == n && size(Z, 2) == n || throw(ArgumentError("S.Z has incompatible size"))
    length(t) == n || throw(ArgumentError("t must have length equal to size(S.T,1)"))
    sort(collect(t)) == collect(1:n) || throw(ArgumentError("t must be a permutation of 1:n"))

    RT = typeof(real(zero(eltype(T))))
    tol = sqrt(eps(RT)) * max(opnorm(T, 1), one(RT))

    # Detect Schur blocks in the initial quasi-triangular form.
    # block_ranges[k] is the index range of block k.
    # block_id[i] tells which block contains index i.
    block_ranges = UnitRange{Int}[]
    block_id = zeros(Int, n)

    i = 1
    bid = 0
    while i <= n
        bid += 1
        if eltype(T) <: Real && i < n && abs(T[i+1, i]) > tol
            push!(block_ranges, i:i+1)
            block_id[i] = bid
            block_id[i+1] = bid
            i += 2
        else
            push!(block_ranges, i:i)
            block_id[i] = bid
            i += 1
        end
    end

    # Validate that every 2x2 real Schur block stays adjacent in t.
    pos_in_t = zeros(Int, n)
    for k in 1:n
        pos_in_t[t[k]] = k
    end

    if eltype(T) <: Real
        for r in block_ranges
            if length(r) == 2
                i, j = first(r), last(r)
                if abs(pos_in_t[i] - pos_in_t[j]) != 1
                    throw(ArgumentError("indices $i and $j belong to the same 2x2 real Schur block and must appear adjacently in t"))
                end
            end
        end
    end

    # Convert entry permutation t into a block permutation tb.
    tb = Int[]
    seen = falses(length(block_ranges))
    for idx in t
        b = block_id[idx]
        if !seen[b]
            push!(tb, b)
            seen[b] = true
        end
    end

    # current[pos] = original block index currently sitting at block position pos
    current = collect(1:length(block_ranges))

    # Realize tb by adjacent block swaps.
    for pos in 1:length(tb)
        want = tb[pos]
        curpos = findfirst(==(want), current)
        curpos === nothing && error("internal permutation tracking failure")

        while curpos > pos
            Aind = block_ranges[curpos - 1]
            Bind = block_ranges[curpos]

            p = length(Aind)
            q = length(Bind)

            a0 = first(Aind)
            b1 = last(Bind)
            rind = a0:b1

            Ablk = Matrix(T[Aind, Aind])
            Bblk = Matrix(T[Bind, Bind])
            Cblk = Matrix(T[Aind, Bind])

            # Solve Sylvester equation
            # Ablk * X - X * Bblk = -Cblk
            K = kron(Matrix{eltype(T)}(I, q, q), Ablk) -
                kron(transpose(Bblk), Matrix{eltype(T)}(I, p, p))
            X = reshape(-(K \ vec(Cblk)), p, q)

            # Build an orthonormal basis whose first q columns span [X; I].
            # Need a full square orthogonal/unitary matrix on the active window.
            Y = vcat(X, Matrix{eltype(T)}(I, q, q))
            U, _, _ = svd(Y; full=true)
            Q = U

            # Similarity update on the active window
            T[rind, a0:end] = Q' * T[rind, a0:end]
            T[:, rind] = T[:, rind] * Q
            Z[:, rind] = Z[:, rind] * Q

            m = p + q

            # Clean numerical fill below first subdiagonal in the active window
            for ii in 1:m
                for jj in 1:(ii - 2)
                    T[rind[ii], rind[jj]] = zero(eltype(T))
                end
            end

            # Zero the off-diagonal block below the diagonal blocks after the swap
            new_top = a0:(a0 + q - 1)
            new_bot = (a0 + q):b1
            T[new_bot, new_top] .= zero(eltype(T))

            # For complex Schur form, restore strict upper triangularity
            if !(eltype(T) <: Real)
                for ii in 2:m
                    T[rind[ii], rind[ii - 1]] = zero(eltype(T))
                end
            end

            # Standardize real 2x2 blocks after the swap
            if eltype(T) <: Real
                # If a 2x2 block exists at the top, eliminate any spurious entries below its first subdiagonal
                if q == 2
                    T[a0+1, a0] = T[a0+1, a0]
                end
                # If a 2x2 block exists at the bottom, same
                if p == 2
                    T[b1, b1-1] = T[b1, b1-1]
                end
            end

            # Update the block partition after adjacent swap
            block_ranges[curpos - 1] = a0:(a0 + q - 1)
            block_ranges[curpos] = (a0 + q):b1
            current[curpos - 1], current[curpos] = current[curpos], current[curpos - 1]
            curpos -= 1
        end
    end

    vals = ComplexF64[]
    i = 1
    while i <= n
        if eltype(T) <: Real && i < n && abs(T[i+1, i]) > tol
            λ = eigvals(Matrix(T[i:i+1, i:i+1]))
            push!(vals, ComplexF64(λ[1]), ComplexF64(λ[2]))
            i += 2
        else
            push!(vals, ComplexF64(T[i, i]))
            i += 1
        end
    end

    return LinearAlgebra.Schur(T, Z, vals)
end
