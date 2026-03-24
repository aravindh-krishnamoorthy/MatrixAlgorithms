################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
################################################################################

using LinearAlgebra
using LinearAlgebra: checksquare

################################################################################
#
# Reordering of eigenvalues of Schur decomposition. Algorithm based on:
#   Daniel Kressner, "Block algorithms for reordering standard and generalized
#   Schur forms," ACM Transactions on Mathematical Software Volume 32 Issue 4
#   pp 521-532 https://doi.org/10.1145/1186785.1186787
# NOTE: Includes inputs from ChatGPT 5.4 Thinking
#
################################################################################
@views @inbounds function ordschur!(S::LinearAlgebra.Schur, p::AbstractVector{<:Integer})
    T = S.T
    Z = S.Z
    vals = S.values

    n = checksquare(T)
    size(Z,1) == n && size(Z,2) == n ||
        throw(ArgumentError("ordschur!: S.Z has incompatible size."))
    length(vals) == n ||
        throw(ArgumentError("ordschur!: S.values has incompatible length."))
    length(p) == n ||
        throw(ArgumentError("ordschur!: permutation p must have length equal to size(S.T,1)."))
    sort!(collect(p)) == collect(1:n) ||
        throw(ArgumentError("ordschur!: permutation p must contain each index 1:n exactly once."))

    # Identify 1x1 and 2x2 blocks
    sizes = ones(Int, n)
    i = 1
    while i < n
        if eltype(T) <: Real && !iszero(T[i+1,i])
            sizes[i] = 2
            sizes[i+1] = 0
            i += 2
        else
            sizes[i] = 1
            i += 1
        end
    end
    if i == n
        sizes[n] = 1
    end

    # A 2x2 real Schur block must move as a unit
    if eltype(T) <: Real
        pinv = similar(p)
        for i in 1:n
            pinv[p[i]] = i
        end
        i = 1
        while i < n
            if sizes[i] == 2
                abs(pinv[i] - pinv[i+1]) == 1 ||
                    throw(ArgumentError("ordschur!: indices $i and $(i+1) belong to the same 2x2 real Schur block and must remain adjacent in p."))
                i += 2
            else
                i += 1
            end
        end
    end

    # Initial block starts and block ids per entry
    nb = count(!iszero, sizes)
    blocks = Vector{Int}(undef, nb)
    block_id = zeros(Int, n)

    b = 0
    for i in 1:n
        if !iszero(sizes[i])
            b += 1
            blocks[b] = i
            block_id[i:i+sizes[i]-1] .= b
        end
    end

    # Desired block permutation induced by p
    pb = Vector{Int}(undef, nb)
    k = 0
    lastb = 0
    for idx in p
        b = block_id[idx]
        if b != lastb
            k += 1
            pb[k] = b
            lastb = b
        end
    end
    k == nb || throw(ArgumentError("ordschur!: invalid block permutation induced by p."))

    # current[pos] = original block id currently sitting at block position pos
    # pos_of[b] = current position of original block id b
    current = collect(1:nb)
    pos_of = collect(1:nb)

    # Workspace
    Δ = Matrix{eltype(T)}(I, 2, 2)
    M_K0  = zeros(eltype(T), 4, 4)
    M_K1  = zeros(eltype(T), 4, 4)
    M_rhs = zeros(eltype(T), 4)
    M_X   = zeros(eltype(T), 2, 2)
    M_Q   = zeros(eltype(T), 4, 4)

    # Realize desired block order by adjacent swaps
    for pos in 1:nb
        want = pb[pos]
        curpos = pos_of[want]

        while curpos > pos
            i = blocks[curpos-1]
            j = blocks[curpos]
            s1, s2 = sizes[i], sizes[j]
            i1, i2 = i, i+s1-1
            j1, j2 = j, j+s2-1
            rind = i1:j2
            m = s1 + s2

            K0  = M_K0[1:s1*s2, 1:s1*s2]
            K1  = M_K1[1:s1*s2, 1:s1*s2]
            rhs = M_rhs[1:s1*s2]
            X   = M_X[1:s1, 1:s2]
            Q   = M_Q[1:m, 1:m]

            # Solve A*X - X*B = -C
            kron!(K0, Δ[1:s2,1:s2], T[i1:i2,i1:i2])
            K0 .-= kron!(K1, transpose(T[j1:j2,j1:j2]), Δ[1:s1,1:s1])
            rhs .= -vec(T[i1:i2,j1:j2])
            ldiv!(lu!(K0), rhs)
            X .= reshape(rhs, s1, s2)

            # Build orthogonal/unitary similarity
            fill!(Q, zero(eltype(T)))
            Q[1:s1, 1:s2] .= X
            Q[s1+1:m, 1:s2] .= Δ[1:s2,1:s2]
            Q[1:s1, s2+1:m] .= Δ[1:s1,1:s1]
            F = qr!(Q)

            lmul!(adjoint(F.Q), T[rind,i1:n])
            rmul!(T[1:n,rind], F.Q)
            rmul!(Z[:,rind], F.Q)

            # Restore Schur structure
            for k in 1:m-2
                T[rind[k+2:end], rind[k]] .= zero(eltype(T))
            end
            T[i1+s2:j2, i1:i1+s2-1] .= zero(eltype(T))
            if !(eltype(T) <: Real)
                T[rind[2:end], rind[1:end-1]] .= zero(eltype(T))
            end

            # Update bookkeeping after swapping adjacent blocks of sizes s1 and s2
            sizes[i1:j2] .= 0
            sizes[i1] = s2
            sizes[i1+s2] = s1
            blocks[curpos-1] = i1
            blocks[curpos] = i1 + s2

            b1 = current[curpos-1]
            b2 = current[curpos]
            current[curpos-1] = b2
            current[curpos] = b1
            pos_of[b1] = curpos
            pos_of[b2] = curpos - 1

            curpos -= 1
        end
    end

    permute!(vals, p)
    return S
end
