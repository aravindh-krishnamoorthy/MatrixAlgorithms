################################################################################
# This file is a part of the Julia package: MatrixFunctions.jl
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
#
# Functions in this file:
#   1. fm_schur_parlett_recurrence(f::Function, X::AbstractMatrix)
#   2. fm_schur_parlett_block(f::Function, X::AbstractMatrix)
#
# In all functions:
# - Intermediate results are computed in the input element type unless
#     conversion to complex type is needed.
################################################################################

using LinearAlgebra
using TaylorSeries

################################################################################
# Schur-Parlett algorithm based on Algorithm 4.13 In
#   Matrix Functions, Higham, 2008, SIAM, ISBN 978-0-89871-646-7.
#
# NOTE: This function implements a basic algorithm that is used internally as 
#       a reference and for testing; this function must not be used in
#       production code.
################################################################################
function fm_schur_parlett_recurrence(f::Function, X::AbstractMatrix)
    m,n = size(X)
    @assert m == n
    @assert m > 0

    # Schur decomposition
    S = schur(X)
    if !istriu(S.T)
        # T is only triangular in the complex field
        S = Schur{Complex}(S)
    end

    # Iterate to obtain the matrices
    T, Z, λ = S

    # For diagonalisable matrices, use the straightforward formula
    if isdiag(T)
        # Move to complex field in case the current domain is insufficient
        # to compute f.(λ)
        λ′ = λ
        try
            λ′ = f.(λ)
        catch e
            if isa(e, DomainError)
                λ = complex.(λ)
                λ′ = f.(λ)
            else
                throw(e)
            end
        end
        D = Diagonal(λ′)
        F = Z*D*Z'
        return F
    end
    
    # Move to complex field in case the current domain is insufficient
    # to compute f.(λ)
    λ′ = λ
    try
        λ′ = f.(λ)
    catch e
        if isa(e, DomainError)
            λ = complex.(λ)
            λ′ = f.(λ)
        else
            throw(e)
        end
    end
    F = zeros(eltype(λ′), n, n)
    F[diagind(F)] = λ′

    # Compute the upper triangular elements via Algorithm 4.13
    for i = 1:n
        F[i,i] = f(λ[i])
    end
    for j = 2:n
        for i = j-1:-1:1
            if λ[i] == λ[j]
                throw(RepeatedEigenvalueException())
            else
                δ = λ[i] - λ[j]
                # TODO: Improve numerical stability
                F[i,j] = T[i,j]*(λ′[i] - λ′[j])/δ + (sum(F[i,i+1:j-1].*T[i+1:j-1,j]) - sum(T[i,i+1:j-1].*F[i+1:j-1,j]))/δ
            end
        end
    end
    return Z*F*Z'
end

################################################################################
#
# Schur-Parlett block algorithm based on Algorithm 9.6 in
#   Matrix Functions, Higham, 2008, SIAM, ISBN 978-0-89871-646-7.
#
# NOTE: This algorithm does not perform reordering and blocking in Steps 3--5
#       (Algorithm 3.6 + confluent permutation search + reordering). Instead
#       this algorithm assumes 'schur' to return an appropriately blocked T.
#
################################################################################
function fm_schur_parlett_block(f::Function, X::AbstractMatrix; MAXITER::Int = 100)
    m,n = size(X)
    @assert m == n
    @assert m > 0

    # Schur decomposition
    S = schur(X)
    T, Z, λ = S.T, S.Z, S.values

    # For diagonalisable matrices, use the straightforward formula
    if isdiag(T)
        # Move to complex field in case the current domain is insufficient
        # to compute f.(λ)
        λ′ = λ
        try
            λ′ = f.(λ)
        catch e
            if isa(e, DomainError)
                λ = complex.(λ)
                λ′ = f.(λ)
            else
                throw(e)
            end
        end
        D = Diagonal(λ′)
        F = Z*D*Z'
        return F
    end

    if ~istriu(T)
        # If T is quasi upper triangular, convert to complex-valued
        # upper triangular as required by Algorithm 9.6
        S = Schur{Complex}(S)
        T, Z, λ = S.T, S.Z, S.values
    end

    # Algorithm 9.5
    # Assign diagonal elements to groups with blocking parameter δ = 0.1
    δ = 0.1
    q = zeros(Int,n)
    ql = zeros(Int,n)
    D = [abs(λ[i] - λ[j]) for i = 1:n, j = 1:n]
    for i = 1:n
        if q[i] == 0
            q[i] = maximum(q) + 1
            ql[q[i]] += 1
        end
        for j = i:n
            if q[j] == 0
                if D[i,j] < δ
                    q[j] = q[i]
                    ql[q[j]] += 1
                end
            end
        end
    end

    #TODO: Reordering. Here, we assume that schur decomposition already
    #TODO: reorders as per requirement.

    display(ql)
    
    ix = 1
    S = fill!(similar(T), 0)
    for i=1:n
        # Break if done
        if ql[i] == 0 break end
        
        # Current m×m diagonal matrix
        ixe = ix+ql[i]-1
        Tii = T[ix:ixe,ix:ixe]
        
        # Algorithm 9.4 for diagonal blocks
        E = Matrix{eltype(Tii)}(I, ql[i], ql[i])
        σ = sum(diag(Tii))/ql[i]
        Mii = Tii - σ*E
        tol = sqrt(eps())
        F0 = f(σ)*E
        F1 = F0
        Pii = Mii
        ϵ = σ + Taylor1(MAXITER)
        ft = f(ϵ)
        for s=1:MAXITER
            F1 = F0  + ft[s]*Pii
            Pii = Pii*Mii/(s+1)
            #TODO: Improve exit condition based on Algorithm 9.4
            if norm(F1 - F0) <= tol*norm(F1)
                break
            end
            F0 = F1
        end
        S[ix:ixe,ix:ixe] = F1

        # Off-diagonal blocks
        jx = j-ql[i]
        for j = i-1:-1:1
            jxe = jx-ql[j]+1
            
        end

        # Next block
        ix = ixe + 1
    end

    F = Z*S*Z'
    return F
end
