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
# Blocked version
################################################################################
function potri2_blocked!(uplo::Char, X::AbstractMatrix{T}; bs::Int=64) where {T}
    n = size(X,1)
    z = zero(T)
    d = zeros(T, n)

    Tbuf = Matrix{T}(undef, bs, bs)
    Pbuf = Matrix{T}(undef, n, bs)

    @inbounds if uplo == 'U'
        for j = 1:n
            Xjj = inv(X[j,j])
            X[j,j] = Xjj
            d[j] = Xjj
            for i = j+1:n
                X[i,j] = z
            end
        end

        @views for je = n:-bs:1
            jb = max(1, je-bs+1)
            J  = jb:je
            mJ = length(J)

            for ke = n:-bs:(je+1)
                kb = max(je+1, ke-bs+1)
                K  = kb:ke

                Irect = 1:jb-1
                if !isempty(Irect)
                    C  = view(X, J, Irect)
                    BK = view(X, K, J)
                    AI = view(X, Irect, K)
                    mul!(C, transpose(BK), transpose(AI), -one(T), one(T))
                end

                if !isempty(K)
                    TJ = view(Tbuf, 1:mJ, 1:mJ)
                    fill!(TJ, z)

                    A = view(X, J, K)
                    B = view(X, K, J)
                    mul!(TJ, A, B, one(T), zero(T))

                    for jj_loc = 1:mJ
                        jj = jb + jj_loc - 1
                        for ii_loc = 1:jj_loc
                            ii = jb + ii_loc - 1
                            X[jj, ii] -= TJ[ii_loc, jj_loc]
                        end
                    end
                end
            end

            for j = je:-1:jb
                if j < je
                    Kb = (j+1):je
                    y  = view(Pbuf, 1:j, 1)
                    mul!(y, view(X, 1:j, Kb), view(X, Kb, j), one(T), zero(T))
                    for i = 1:j
                        X[j,i] -= y[i]
                    end
                end

                for k = j:-1:1
                    tmp = X[j,k] * d[k]
                    X[j,k] = conj(tmp)
                    for i = 1:k-1
                        X[j,i] = X[j,i] - X[i,k] * tmp
                    end
                end
            end
        end

        @inbounds for i = 1:n, j = i+1:n
            X[i,j] = X[j,i]'
        end
        return X

    else
        for j = 1:n
            Xjj = inv(X[j,j])
            X[j,j] = Xjj
            d[j] = Xjj
            for i = 1:j-1
                X[i,j] = z
            end
        end

        @views for je = n:-bs:1
            jb = max(1, je-bs+1)
            J  = jb:je
            mJ = length(J)

            for ke = n:-bs:(je+1)
                kb = max(je+1, ke-bs+1)
                K  = kb:ke

                Irect = 1:jb-1
                if !isempty(Irect)
                    C = view(X, Irect, J)
                    A = view(X, K, Irect)
                    B = view(X, J, K)
                    mul!(C, transpose(A), transpose(B), -one(T), one(T))
                end

                if !isempty(K)
                    TJ = view(Tbuf, 1:mJ, 1:mJ)
                    fill!(TJ, z)

                    A = view(X, K, J)
                    B = view(X, J, K)
                    mul!(TJ, transpose(A), transpose(B), one(T), zero(T))

                    for jj_loc = 1:mJ
                        jj = jb + jj_loc - 1
                        for ii_loc = 1:jj_loc
                            ii = jb + ii_loc - 1
                            X[ii, jj] -= TJ[ii_loc, jj_loc]
                        end
                    end
                end
            end

            for j = je:-1:jb
                if j < je
                    Kb = (j+1):je
                    m  = length(Kb)

                    v = view(Pbuf, 1:m, 1)
                    for t = 1:m
                        v[t] = X[j, first(Kb) + t - 1]
                    end

                    y = view(Pbuf, 1:j, 2)
                    mul!(y, transpose(view(X, Kb, 1:j)), v, one(T), zero(T))

                    for i = 1:j
                        X[i,j] -= y[i]
                    end
                end

                for k = j:-1:1
                    tmp = X[k,j] * d[k]
                    X[k,j] = conj(tmp)
                    for i = 1:k-1
                        X[i,j] = X[i,j] - X[k,i] * tmp
                    end
                end
            end
        end

        @inbounds for i = 1:n, j = 1:i-1
            X[i,j] = X[j,i]'
        end
        return X
    end
end
