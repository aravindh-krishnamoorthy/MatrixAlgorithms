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
function sqrtm(A::AbstractMatrix{T}, f::Function = trsr!) where {T}
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
        negative = false
        for i = 1:n
            if isreal(S.values[i]) && real(S.values[i]) < 0
                negative = true
                break
            end
        end
        if negative
            S = Schur{Complex}(S)
        end
        return S.Z*f(S.T)*S.Z'
    else # complex A
        S = schur(A)
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
# NOTE: The function was modified to use Schur output instead of Jordon form using
#       ChatGPT 5.4 Thinking
#
################################################################################
function gantmacher!(A::AbstractMatrix{T}; atol::Real = 0,
                     rtol::Real = min(size(A)...)*eps(real(float(one(T))))) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("gantmacher!: Matrix A must be square."))

    RT = typeof(real(float(one(T))))
    Tc = Complex{RT}

    Ac = Matrix{Tc}(undef, n, n)
    @inbounds for j = 1:n
        for i = 1:n
            Ac[i, j] = A[i, j]
        end
    end

    S = schur(Ac)

    scale_vals = zero(RT)
    @inbounds for i = 1:length(S.values)
        ai = abs(S.values[i])
        if ai > scale_vals
            scale_vals = ai
        end
    end
    scale_mat = opnorm(Ac, 1)
    scale = max(scale_vals, scale_mat)
    tol0 = max(atol, rtol * scale)

    select = Vector{Bool}(undef, n)
    r = 0
    @inbounds for i = 1:n
        keep = abs(S.values[i]) > tol0
        select[i] = keep
        if keep
            r += 1
        end
    end

    S = ordschur(S, select)

    Tord = Matrix{Tc}(S.T)
    Qord = Matrix{Tc}(S.Z)

    if r == n
        return Qord * trsr!(copy(Tord)) * Qord'
    end

    m0 = n - r

    T1 = Matrix{Tc}(undef, r, r)
    T0 = Matrix{Tc}(undef, m0, m0)
    B  = Matrix{Tc}(undef, r, m0)

    if r > 0
        @views copyto!(T1, Tord[1:r, 1:r])
        @views copyto!(B,  Tord[1:r, r+1:n])
    end
    if m0 > 0
        @views copyto!(T0, Tord[r+1:n, r+1:n])
    end

    U1 = r == 0 ? zeros(Tc, 0, 0) : trsr!(copy(T1))

    if m0 == 0
        return Qord * U1 * Qord'
    end

    K = Vector{Matrix{Tc}}(undef, m0 + 1)
    ν = Vector{Int}(undef, m0 + 1)
    K[1] = zeros(Tc, m0, 0)
    ν[1] = 0

    Pk = Matrix{Tc}(I, m0, m0)
    stabilized = false
    for k = 1:m0
        Pk = Pk * T0
        N = nullspace(Pk; atol=atol, rtol=rtol)
        Nk = Matrix{Tc}(N)
        K[k + 1] = Nk
        νk = size(Nk, 2)
        ν[k + 1] = νk
        if νk == m0
            for j = k+1:m0
                K[j + 1] = Nk
                ν[j + 1] = m0
            end
            stabilized = true
            break
        end
    end
    stabilized || throw(ArgumentError("gantmacher!: failed to identify a numerically nilpotent zero-primary block."))

    d = Vector{Int}(undef, m0)
    @inbounds for k = 1:m0
        d[k] = ν[k + 1] - ν[k]
    end

    b = Vector{Int}(undef, m0)
    @inbounds for k = 1:m0-1
        b[k] = d[k] - d[k + 1]
    end
    b[m0] = d[m0]

    powers = Vector{Matrix{Tc}}(undef, m0 + 1)
    powers[1] = Matrix{Tc}(I, m0, m0)
    for k = 1:m0
        powers[k + 1] = powers[k] * T0
    end

    heads = Vector{Matrix{Tc}}(undef, m0)
    for k = 1:m0
        heads[k] = zeros(Tc, m0, 0)
    end

    for k = m0:-1:1
        need = b[k]
        need == 0 && continue

        gcols = 0
        for j = k+1:m0
            gcols += size(heads[j], 2)
        end

        G = zeros(Tc, m0, gcols)
        col = 1
        for j = k+1:m0
            Hj = heads[j]
            hjcols = size(Hj, 2)
            if hjcols > 0
                block = powers[j - k + 1] * Hj
                @views G[:, col:col + hjcols - 1] .= block
                col += hjcols
            end
        end

        Kk  = K[k]
        Kk1 = K[k + 1]
        nk  = size(Kk, 2)
        nk1 = size(Kk1, 2)

        Sspan = zeros(Tc, m0, nk + gcols)
        if nk > 0
            @views Sspan[:, 1:nk] .= Kk
        end
        if gcols > 0
            @views Sspan[:, nk+1:nk+gcols] .= G
        end

        nk1 == 0 && throw(ArgumentError("gantmacher!: numerically ambiguous nilpotent chain structure."))

        local Sorth
        if size(Sspan, 2) > 0
            F = svd(Sspan; full=false)
            σ = F.S
            σmax = isempty(σ) ? zero(RT) : maximum(abs, σ)
            tolS = max(atol, rtol * σmax)
            rs = 0
            @inbounds for i = 1:length(σ)
                if abs(σ[i]) > tolS
                    rs += 1
                end
            end
            Sorth = rs == 0 ? zeros(Tc, m0, 0) : Matrix{Tc}(F.U[:, 1:rs])
        else
            Sorth = zeros(Tc, m0, 0)
        end

        Y = copy(Kk1)
        if size(Sorth, 2) > 0
            Y .-= Sorth * (Sorth' * Y)
        end

        FY = svd(Y; full=false)
        σY = FY.S
        σYmax = isempty(σY) ? zero(RT) : maximum(abs, σY)
        tolY = max(atol, rtol * σYmax)

        rloc = 0
        @inbounds for i = 1:length(σY)
            if abs(σY[i]) > tolY
                rloc += 1
            end
        end

        rloc < need && throw(ArgumentError("gantmacher!: numerically ambiguous nilpotent chain structure."))
        heads[k] = need == 0 ? zeros(Tc, m0, 0) : Matrix{Tc}(FY.U[:, 1:need])
    end

    nchains = 0
    for k = 1:m0
        nchains += size(heads[k], 2)
    end
    nchains == 0 && throw(ArgumentError("gantmacher!: failed to construct any nilpotent chains."))

    chains = Vector{Vector{Tc}}(undef, m0)
    parts = Vector{Int}(undef, nchains)

    chain_idx = 1
    part_idx = 1
    for k = m0:-1:1
        H = heads[k]
        hcols = size(H, 2)
        for j = 1:hcols
            h = @view H[:, j]
            for p = k-1:-1:0
                chains[chain_idx] = powers[p + 1] * h
                chain_idx += 1
            end
            parts[part_idx] = k
            part_idx += 1
        end
    end

    V = Matrix{Tc}(undef, m0, m0)
    for j = 1:m0
        V[:, j] = chains[j]
    end

    size(V, 2) == m0 || throw(ArgumentError("gantmacher!: failed to construct a full nilpotent chain basis."))
    rank(V; atol=atol, rtol=rtol) < m0 &&
        throw(ArgumentError("gantmacher!: numerically singular chain basis for the zero-primary block."))

    sort!(parts, rev=true)

    npieces = 0
    i = 1
    while i <= length(parts)
        if i == length(parts)
            parts[i] == 1 || throw(ArgumentError("gantmacher!: inferred nilpotent block partition does not admit a square root."))
            npieces += 1
            i += 1
        else
            s = parts[i]
            t = parts[i + 1]
            (s - t <= 1) || throw(ArgumentError("gantmacher!: inferred nilpotent block partition does not admit a square root."))
            npieces += 1
            i += 2
        end
    end

    pieces = Vector{Matrix{Tc}}(undef, npieces)
    piece_idx = 1
    i = 1
    while i <= length(parts)
        if i == length(parts)
            pieces[piece_idx] = zeros(Tc, 1, 1)
            piece_idx += 1
            i += 1
        else
            s = parts[i]
            t = parts[i + 1]
            rpair = s + t

            Jnil = zeros(Tc, rpair, rpair)
            for ii = 1:rpair-1
                Jnil[ii, ii + 1] = one(Tc)
            end

            perm = Vector{Int}(undef, rpair)
            idx = 1
            for ii = 1:2:rpair
                perm[idx] = ii
                idx += 1
            end
            for ii = 2:2:rpair
                perm[idx] = ii
                idx += 1
            end

            Pperm = zeros(Tc, rpair, rpair)
            for j = 1:rpair
                Pperm[perm[j], j] = one(Tc)
            end

            pieces[piece_idx] = Pperm' * Jnil * Pperm
            piece_idx += 1
            i += 2
        end
    end

    YJ = zeros(Tc, m0, m0)
    p = 1
    for idx = 1:npieces
        Bi = pieces[idx]
        q = p + size(Bi, 1) - 1
        @views YJ[p:q, p:q] .= Bi
        p = q + 1
    end

    size(YJ, 1) == m0 || throw(ArgumentError("gantmacher!: internal size mismatch in nilpotent model block."))

    U0 = V * YJ / V

    local Y
    if r == 0
        Y = zeros(Tc, 0, m0)
    elseif m0 == 0
        Y = zeros(Tc, r, 0)
    else
        I0 = Matrix{Tc}(I, m0, m0)
        I1 = Matrix{Tc}(I, r, r)
        Ksys = kron(I0, U1) + kron(transpose(U0), I1)
        y = Ksys \ vec(B)
        Y = reshape(y, r, m0)
    end

    Uord = zeros(Tc, n, n)
    if r > 0
        @views Uord[1:r, 1:r] .= U1
        @views Uord[1:r, r+1:n] .= Y
    end
    @views Uord[r+1:n, r+1:n] .= U0

    X = Qord * Uord * Qord'

    A_scale = max(opnorm(Ac, 2), one(RT))
    res = opnorm(X * X - Ac, 2) / A_scale
    res <= max(sqrt(eps(RT)), 10 * max(atol, rtol)) ||
        throw(ArgumentError("gantmacher!: hybrid Schur-Gantmacher construction failed residual check."))

    return X
end
