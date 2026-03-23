################################################################################
#
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
#
################################################################################

import Base: sqrt
using LinearAlgebra
using JordanForm

################################################################################
#
# Square root of a matrix
# Reference [1]: Smith, M. I. (2003). A Schur Algorithm for Computing Matrix pth Roots.
#   SIAM Journal on Matrix Analysis and Applications (Vol. 24, Issue 4, pp. 971–989).
#   https://doi.org/10.1137/s0895479801392697
#
# NOTE: The matrix is lifted to the complex field if the diagonal values of the
#         triangular Schur matrix does not have a square root in type T
#
################################################################################
function sqrtm(A::AbstractMatrix{T}; atol::Real = 0,
               rtol::Real = min(size(A)...)*eps(real(float(one(T))))) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("sqrt: Matrix A must be square."))
    symmetric = ishermitian(A)
    if symmetric
        e, V = eigen(A)
        negative = any(e .< 0)
        if negative
            Q = Diagonal(sqrt.(complex.(e)))
            if isreal(V)
                return complex.(V*real(Q)*V', V*imag(Q)*V')
            else
                return V*Q*V'
            end
        else
            Q = Diagonal(sqrt.(e))
            return V*Q*V'
        end
    elseif isreal(A)
        S = schur(real(A))
        scale = isempty(S.values) ? zero(real(float(one(T)))) : maximum(abs, S.values)
        tol = max(atol, rtol*scale)
        negative = false
        for i = 1:n
            if isreal(S.values[i]) && real(S.values[i]) < -tol
                negative = true
                break
            end
        end
        if negative
            S = Schur{Complex}(S)
            scale = isempty(S.values) ? zero(real(float(one(eltype(S.values))))) : maximum(abs, S.values)
            tol = max(atol, rtol*scale)
        end
        rankA = count(λ -> abs(λ) > tol, S.values)
        f = rankA < n - 1 ? X -> gantmacher!(X; atol=atol, rtol=rtol) : trsr!
        return S.Z*f(S.T)*S.Z'
    else # complex A
        S = schur(A)
        scale = isempty(S.values) ? zero(real(float(one(T)))) : maximum(abs, S.values)
        tol = max(atol, rtol*scale)
        rankA = count(λ -> abs(λ) > tol, S.values)
        f = rankA < n - 1 ? X -> gantmacher!(X; atol=atol, rtol=rtol) : trsr!
        return S.Z*f(S.T)*S.Z'
    end
end

################################################################################
#
# Square root of a quasi upper triangular matrix (output of Schur decomposition)
# Reference [1]: Smith, M. I. (2003). A Schur Algorithm for Computing Matrix pth Roots.
#   SIAM Journal on Matrix Analysis and Applications (Vol. 24, Issue 4, pp. 971–989).
#   https://doi.org/10.1137/s0895479801392697
#
# VERSION: Pure Julia version for both real- and complex-valued inputs
# NOTE: It is assumed that the diagonal elements of A have a square root in type T
#
################################################################################
@views @inbounds function trsr!(A::AbstractMatrix{T}) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("trsr!: Matrix A must be square."))
    # Choose complex or real dot product based on T
    dot = T <: Complex ? BLAS.dotu : BLAS.dot
    # Square roots of 1x1 and 2x2 diagonal blocks
    i = 1
    sizes = ones(Int,n)
    while i < n
        if !iszero(A[i+1,i])
            μ = sqrt(-real(A[i,i+1]*A[i+1,i]))
            r = sqrt(hypot(A[i,i], μ))
            θ = atan(μ, real(A[i,i]))
            s, c = sincos(θ/2)
            α, β′ = r*c, r*s/µ
            A[i,i] = α
            A[i+1,i+1] = α
            A[i,i+1] = β′*A[i,i+1]
            A[i+1,i] = β′*A[i+1,i]
            sizes[i] = 2
            sizes[i+1] = 0
            i += 2
        else
            A[i,i] = sqrt(A[i,i])
            sizes[i] = 1
            i += 1
        end
    end
    if i == n
        A[n,n] = sqrt(A[n,n])
        sizes[i] = 1
    end
    # Algorithm 4.3 in Reference [1]
    Δ = I(4)
    M_L₀ = zeros(T,4,4)
    M_L₁ = zeros(T,4,4)
    for k = 1:n-1
        for i = 1:n-k
            if sizes[i] == 0 || sizes[i+k] == 0 continue end
            i₁, i₂, j₁, j₂, s₁, s₂ = i, i+sizes[i]-1, i+k, i+k+sizes[i+k]-1, sizes[i], sizes[i+k]
            k₁, k₂ = i+s₁, i+k-1
            L₀ = M_L₀[1:s₁*s₂,1:s₁*s₂]
            L₁ = M_L₁[1:s₁*s₂,1:s₁*s₂]
            if s₁ == 1 && s₂ == 1
                Bᵢⱼ⁽⁰⁾ = dot(A[i₁,k₁:k₂], A[k₁:k₂,j₁])
                A[i₁,j₁] = (A[i₁,j₁] - Bᵢⱼ⁽⁰⁾)/(A[i₁,i₁] + A[j₁,j₁])
            else
                # Compute Bᵢⱼ⁽⁰⁾ and update A[i₁:i₂,j₁:j₂]
                mul!(A[i₁:i₂,j₁:j₂], A[i₁:i₂,k₁:k₂], A[k₁:k₂,j₁:j₂], T(-1.0), T(+1.0))
                # Solve Uᵢ,ᵢ₊ₖ using Reference [1, (4.10)]
                kron!(L₀, Δ[1:s₂,1:s₂], A[i₁:i₂,i₁:i₂])
                L₀ .+= kron!(L₁, transpose(A[j₁:j₂,j₁:j₂]), Δ[1:s₁,1:s₁])
                ldiv!(lu!(L₀), A[i₁:i₂,j₁:j₂][:])
            end
        end
    end
    return A
end

################################################################################
#
# Square root of a matrix using Gantmacher's algorithm
# Reference [2]: Gantmacher, F. R. The Theory of Matrices. Vol. 2, Chapter VIII, Sec. 7.
# NOTE: Jordan form version (numerically unstable) for singular/rank-deficient inputs
#       and will need to be improved.
# NOTE: This first version of the code was generated by ChatGPT 5.2 Thinking
#
################################################################################
function gantmacher!(A::AbstractMatrix{T}; atol::Real = 0,
                     rtol::Real = min(size(A)...)*eps(real(float(one(T))))) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("gantmacher!: Matrix A must be square."))

    S, J = jordan_form(A)
    J = Matrix(J)
    S = Matrix(S)

    RT = real(float(one(T)))
    CT = Complex{RT}

    function blockdiag(blocks::Vector{<:AbstractMatrix})
        N = sum(size(B, 1) for B in blocks)
        X = zeros(CT, N, N)
        p = 1
        for B in blocks
            q = p + size(B, 1) - 1
            X[p:q, p:q] .= CT.(B)
            p = q + 1
        end
        return X
    end

    function jordanblock(λ, k::Int)
        B = zeros(CT, k, k)
        for i = 1:k
            B[i, i] = λ
        end
        for i = 1:k-1
            B[i, i+1] = one(CT)
        end
        return B
    end

    function parseblocks(JM::AbstractMatrix)
        vals = diag(JM)
        scale = isempty(vals) ? zero(RT) : maximum(abs, vals)
        tol = max(atol, rtol*scale)

        blocks = Tuple{Int,Int,eltype(JM)}[]
        i = 1
        while i <= size(JM, 1)
            λ = JM[i, i]
            j = i
            while j < size(JM, 1) &&
                  abs(JM[j+1, j+1] - λ) <= tol &&
                  abs(JM[j, j+1] - one(eltype(JM))) <= tol
                j += 1
            end
            push!(blocks, (i, j - i + 1, λ))
            i = j + 1
        end
        return blocks, tol
    end

    function nonzero_block_sqrt(JB::AbstractMatrix)
        k = size(JB, 1)
        λ = CT(JB[1, 1])
        μ = sqrt(λ)
        N = CT.(JB) - λ*Matrix{CT}(I, k, k)
        X = zeros(CT, k, k)
        term = Matrix{CT}(I, k, k)
        coeff = one(CT)
        X .+= μ*term
        for j = 1:k-1
            term = term*(N/λ)
            coeff *= (CT(1)/CT(2) - CT(j - 1))/CT(j)
            X .+= μ*coeff*term
        end
        return X
    end

    function zeropairs(sizes::Vector{Int})
        s = sort(copy(sizes), rev=true)
        pairs = Tuple{Int,Int}[]
        i = 1
        while i <= length(s)
            if i == length(s)
                s[i] == 1 || throw(ArgumentError("gantmacher!: zero Jordan blocks do not admit a square root."))
                push!(pairs, (1, 0))
                i += 1
            else
                a, b = s[i], s[i+1]
                (a - b <= 1) || throw(ArgumentError("gantmacher!: zero Jordan blocks do not admit a square root."))
                push!(pairs, (a, b))
                i += 2
            end
        end
        return pairs
    end

    function zero_pair_root(s::Int, t::Int)
        if t == 0
            return zeros(CT, 1, 1)
        end
        r = s + t
        Jnil = jordanblock(zero(CT), r)
        perm = vcat(collect(1:2:r), collect(2:2:r))
        P = zeros(CT, r, r)
        for j = 1:r
            P[perm[j], j] = one(CT)
        end
        return P' * Jnil * P
    end

    blocks, tol = parseblocks(J)

    zero_ids = Int[]
    nonzero_ids = Int[]
    for (k, (_, sz, λ)) in enumerate(blocks)
        if abs(λ) <= tol
            push!(zero_ids, k)
        else
            push!(nonzero_ids, k)
        end
    end

    zero_ids = sort(zero_ids; by = i -> blocks[i][2], rev=true)
    order = vcat(nonzero_ids, zero_ids)

    perm = Int[]
    for idx in order
        i, sz, _ = blocks[idx]
        append!(perm, i:i+sz-1)
    end
    P = zeros(CT, n, n)
    for j = 1:n
        P[perm[j], j] = one(CT)
    end

    Jp = P' * CT.(J) * P
    pblocks, _ = parseblocks(Jp)

    pieces = Matrix{CT}[]

    nz = length(nonzero_ids)
    for k = 1:nz
        i, sz, _ = pblocks[k]
        JB = Jp[i:i+sz-1, i:i+sz-1]
        push!(pieces, nonzero_block_sqrt(JB))
    end

    zero_sizes = [blocks[i][2] for i in zero_ids]
    pairs = zeropairs(zero_sizes)

    for (s, t) in pairs
        push!(pieces, zero_pair_root(s, t))
    end

    Xp = blockdiag(pieces)
    XJ = P * Xp * P'

    return CT.(S) * XJ / CT.(S)
end
