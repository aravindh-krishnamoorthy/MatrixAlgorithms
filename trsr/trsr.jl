################################################################################
#
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
#
################################################################################

import Base: sqrt
using LinearAlgebra

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
        rankA = count(λ -> abs(λ) > tol, S.values)
        use_gantmacher = rankA < n - 1
        if negative || use_gantmacher
            S = Schur{Complex}(S)
        end
        if use_gantmacher
            return S.Z*gantmacher!(S.T; atol=atol, rtol=rtol)*S.Z'
        else
            return S.Z*trsr!(S.T)*S.Z'
        end
    else # complex A
        S = schur(A)
        scale = isempty(S.values) ? zero(real(float(one(T)))) : maximum(abs, S.values)
        tol = max(atol, rtol*scale)
        rankA = count(λ -> abs(λ) > tol, S.values)
        if rankA < n - 1
            return S.Z*gantmacher!(S.T; atol=atol, rtol=rtol)*S.Z'
        else
            return S.Z*trsr!(S.T)*S.Z'
        end
    end
end

################################################################################
#
# Square root of a matrix using a hybrid Schur-Gantmacher algorithm
#
# VERSION: Floating-point hybrid using Schur reordering for the nonzero spectrum,
#          numerical recovery of the nilpotent partition from nullities of powers,
#          and a Gantmacher-style construction on the zero-primary block.
# NOTE: This branch is tolerance-dependent and intended for singular / near-singular
#       floating-point inputs. It may fail when the zero-eigenvalue structure is
#       numerically ambiguous.
#
################################################################################
function gantmacher!(A::AbstractMatrix{T}; atol::Real = 0,
                     rtol::Real = min(size(A)...)*eps(real(float(one(T))))) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("gantmacher!: Matrix A must be square."))

    RT = typeof(real(float(one(T))))
    Tc = Complex{RT}
    Ac = Matrix{Tc}(A)

    function _basis_complement(B, S; need::Union{Nothing,Int}=nothing)
        Bc = Matrix{Tc}(B)
        Sc = Matrix{Tc}(S)
        if size(Bc, 2) == 0
            return zeros(Tc, size(Bc, 1), 0)
        end
        Y = copy(Bc)
        if size(Sc, 2) > 0
            Y .-= Sc * (Sc' * Y)
        end
        F = qr(Y)
        Q = Matrix(F.Q)
        R = Matrix(F.R)
        r = 0
        t = max(atol, rtol * (isempty(R) ? zero(RT) : maximum(abs, diag(R))))
        for j = 1:min(size(R,1), size(R,2))
            if abs(R[j,j]) > t
                r += 1
            end
        end
        if need !== nothing
            r < need && throw(ArgumentError("gantmacher!: numerically ambiguous nilpotent chain structure."))
            r = need
        end
        return Q[:, 1:r]
    end

    function _nullity_basis(M)
        N = nullspace(M; atol=atol, rtol=rtol)
        return Matrix{Tc}(N), size(N, 2)
    end

    function _blockdiag(blocks)
        N = sum(size(B, 1) for B in blocks)
        X = zeros(Tc, N, N)
        p = 1
        for B in blocks
            q = p + size(B, 1) - 1
            X[p:q, p:q] .= Matrix{Tc}(B)
            p = q + 1
        end
        return X
    end

    function _jordanblock(λ, k::Int)
        J = zeros(Tc, k, k)
        for i = 1:k
            J[i, i] = λ
        end
        for i = 1:k-1
            J[i, i+1] = one(Tc)
        end
        return J
    end

    S = schur(Ac)
    scale = isempty(S.values) ? zero(RT) : maximum(abs, S.values)
    tol0 = max(atol, rtol*scale)
    select = map(λ -> abs(λ) > tol0, S.values)
    S = ordschur(S, select)

    r = count(select)
    Tord = Matrix(S.T)
    Qord = Matrix(S.Z)

    if r == n
        return Qord * trsr!(copy(Tord)) * Qord'
    elseif r == 0
        T1 = zeros(Tc, 0, 0)
        B = zeros(Tc, 0, n)
        T0 = Tord
    else
        T1 = Tord[1:r, 1:r]
        B = Tord[1:r, r+1:n]
        T0 = Tord[r+1:n, r+1:n]
    end

    U1 = r == 0 ? zeros(Tc, 0, 0) : trsr!(copy(T1))

    m0 = size(T0, 1)
    if m0 == 0
        Uord = U1
        return Qord * Uord * Qord'
    end

    K = Vector{Matrix{Tc}}(undef, m0 + 1)
    ν = zeros(Int, m0 + 1)
    K[1] = zeros(Tc, m0, 0)
    ν[1] = 0

    P = Matrix{Tc}(I, m0, m0)
    stabilized = false
    for k = 1:m0
        P = P * T0
        K[k+1], ν[k+1] = _nullity_basis(P)
        if ν[k+1] == m0
            for j = k+1:m0
                K[j+1] = K[k+1]
                ν[j+1] = m0
            end
            stabilized = true
            break
        end
    end
    stabilized || throw(ArgumentError("gantmacher!: failed to identify a numerically nilpotent zero-primary block."))

    d = zeros(Int, m0)
    for k = 1:m0
        d[k] = ν[k+1] - ν[k]
    end

    b = zeros(Int, m0)
    for k = 1:m0-1
        b[k] = d[k] - d[k+1]
    end
    b[m0] = d[m0]

    powers = Vector{Matrix{Tc}}(undef, m0 + 1)
    powers[1] = Matrix{Tc}(I, m0, m0)
    for k = 1:m0
        powers[k+1] = powers[k] * T0
    end

    heads = [zeros(Tc, m0, 0) for _ = 1:m0]
    for k = m0:-1:1
        need = b[k]
        if need == 0
            continue
        end
        G = zeros(Tc, m0, 0)
        for j = k+1:m0
            if size(heads[j], 2) > 0
                G = hcat(G, powers[j-k+1] * heads[j])
            end
        end
        Sspan = hcat(K[k], G)
        heads[k] = _basis_complement(K[k+1], Sspan; need=need)
    end

    chains = Vector{Vector{Tc}}()
    parts = Int[]
    for k = m0:-1:1
        H = heads[k]
        for j = 1:size(H, 2)
            h = H[:, j]
            for p = k-1:-1:0
                push!(chains, powers[p+1] * h)
            end
            push!(parts, k)
        end
    end
    V = hcat(chains...)
    size(V, 2) == m0 || throw(ArgumentError("gantmacher!: failed to construct a full nilpotent chain basis."))

    if rank(V; atol=atol, rtol=rtol) < m0
        throw(ArgumentError("gantmacher!: numerically singular chain basis for the zero-primary block."))
    end

    sort!(parts, rev=true)
    pieces = Matrix{Tc}[]
    i = 1
    while i <= length(parts)
        if i == length(parts)
            parts[i] == 1 || throw(ArgumentError("gantmacher!: inferred nilpotent block partition does not admit a square root."))
            push!(pieces, zeros(Tc, 1, 1))
            i += 1
        else
            s, t = parts[i], parts[i+1]
            (s - t <= 1) || throw(ArgumentError("gantmacher!: inferred nilpotent block partition does not admit a square root."))
            rpair = s + t
            Jnil = _jordanblock(zero(Tc), rpair)
            perm = vcat(collect(1:2:rpair), collect(2:2:rpair))
            Pperm = zeros(Tc, rpair, rpair)
            for j = 1:rpair
                Pperm[perm[j], j] = one(Tc)
            end
            push!(pieces, Pperm' * Jnil * Pperm)
            i += 2
        end
    end

    YJ = _blockdiag(pieces)
    U0 = V * YJ / V

    if r == 0
        Uord = U0
    else
        Ksys = kron(Matrix{Tc}(I, m0, m0), U1) + kron(transpose(U0), Matrix{Tc}(I, r, r))
        y = Ksys \ vec(B)
        Y = reshape(y, r, m0)
        Uord = [U1 Y; zeros(Tc, m0, r) U0]
    end

    X = Qord * Uord * Qord'

    A_scale = max(opnorm(Ac, 2), one(RT))
    res = opnorm(X * X - Ac, 2) / A_scale
    res <= max(sqrt(eps(RT)), 10 * max(atol, rtol)) ||
        throw(ArgumentError("gantmacher!: hybrid Schur-Gantmacher construction failed residual check."))

    return X
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
            α, β′ = r*c, r*s/μ
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
            if sizes[i] == 0 || sizes[i+k] == 0
                continue
            end
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
