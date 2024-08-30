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
function fm_schur_parlett_block(f::Function, X::AbstractMatrix)
    m,n = size(X)
    @assert m == n
    @assert m > 0

    # Schur decomposition
    T, Z, λ = schur(X)

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

    # Assign diagonal elements to groups with blocking parameter δ = 0.1
    δ = 0.1
    S = ones(n,1)
    C = zeros(n,1)
    C[1] = 1
    for i = 2:n
        if norm(λ[i] - λ[i-1]) ≤ δ
            S[i] = S[i-1]
            C[S[i]] = C[S[i]] + 1
            continue
        else
            S[i] = S[i-1] + 1
            C[S[i]] = 1
        end
    end
    C = cumsum(C)
    mb = S[n]
end
