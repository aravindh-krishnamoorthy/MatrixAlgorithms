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
# Square root of a quasi upper triangular matrix (output of Schur decomposition)
# Reference [1]: Smith, M. I. (2003). A Schur Algorithm for Computing Matrix pth Roots.
#   SIAM Journal on Matrix Analysis and Applications (Vol. 24, Issue 4, pp. 971–989).
#   https://doi.org/10.1137/s0895479801392697
#
# VERSION: Eventual FORTRAN version for both real- and complex-valued inputs
# NOTE: It is assumed that the diagonal elements of A have a square root in type T
#
################################################################################
@views @inbounds function trsrf!(A::AbstractMatrix{T}) where {T}
    m, n = size(A)
    (m == n) || throw(ArgumentError("trsrf!: Matrix A must be square."))
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
                BLAS.gemm!('N', 'N', T(-1.0), A[i₁:i₂,k₁:k₂], A[k₁:k₂,j₁:j₂], T(+1.0), A[i₁:i₂,j₁:j₂])
                # Solve Uᵢ,ᵢ₊ₖ
                _, scale = LAPACK.trsyl!('N', 'N', A[i₁:i₂,i₁:i₂], A[j₁:j₂,j₁:j₂], A[i₁:i₂,j₁:j₂])
                rmul!(A[i₁:i₂,j₁:j₂], inv(scale))
            end
        end
    end
    return A
end
