################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
# Contains contributions by ChatGPT 5.2 Thinking
################################################################################

################################################################################
# Invert a positive definite matrix using the algorithm in
# "Matrix Inversion Using Cholesky Decomposition", by Aravindh Krishnamoorthy
#   and Deepak Menon, arXiv:1111.4144.
################################################################################
# Parallel version
################################################################################
function potri2_bd!(uplo::Char, X::AbstractMatrix{T}, j::Integer, jb::Integer) where {T}
    n = size(X,1)
    if uplo == 'U'
    else # uplo == 'L'
        for k = j:-1:j-jb+1
            for i = j-jb+1:k
                X[i,k] = X[i,k] - dot(X[k+1:n,i], X[k,k+1:n])
            end
            for i = k-1:-1:j-jb+1
                X[i,k] = X[i,k] - dot(X[i+1:k,i], X[i+1:k,k])
            end
        end
    end
end

function potri2_bo!(uplo::Char, X::AbstractMatrix{T}, i::Integer, ib::Integer, j::Integer, jb::Integer) where {T}
    n = size(X,1)
    if uplo == 'U'
    else # uplo == 'L'
        for k = j-jb+1:j
            if k < n
                X[i-ib+1:i,k] = X[i-ib+1:i,k] - X[k+1:n,i-ib+1:i]'*X[k:k,k+1:n]'
            end
        end
        for l = i:-1:i-ib+1
            for k = j-jb+1:j
                X[l,k] = X[l,k] - dot(X[l+1:k,l], X[l+1:k,k])
            end
        end
    end
end

function potri2_parallel!(uplo::Char, X::AbstractMatrix{T}) where {T}
    nb = 32
    n = size(X,1)
    if uplo == 'U'
        return potri2!(uplo, X)
    else # uplo == 'L'
        for j = 1:n
            X[j,j] = 1/X[j,j]
            X[j+1:n,j] = X[j+1:n,j]*X[j,j]
            X[j,j+1:n] .= 0
            X[j,j] = X[j,j]*X[j,j]
        end
        # Anti-diagonal blocks can be run parallelly.
        ib = 0
        for j = n:-nb:1
            Threads.@threads for k = n:-nb:j
                ii = k
                ib = min(ii,nb)
                jj = j+(n-k)
                jb = min(jj,nb)
                if ii == jj
                    potri2_bd!(uplo, X, jj, jb)
                elseif ii > jj
                    potri2_bo!(uplo, X, jj, jb, ii, ib)
                end
            end
        end
        beg = ib
        for i = n-nb:-nb:1
            Threads.@threads for k = beg:nb:i
                ii = i-(k-beg)
                ib = min(ii,nb)
                jj = k
                jb = min(jj,nb)
                if ii == jj
                    potri2_bd!(uplo, X, jj, jb)
                elseif ii > jj
                    potri2_bo!(uplo, X, jj, jb, ii, ib)
                end
            end
        end
        return X
    end
end
