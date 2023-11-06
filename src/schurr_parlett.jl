# Schur-Parlett algorithm based on Algorithm 9.6 in Matrix Functions, Higham.
function fm_schur_parlett(f::Function, X::AbstractMatrix)
    T, Q = schur(X)
    # For diagonalisable matrices, use the straightforward formula
    if isdiag(T)
        F = Q*f.(T)*Q'
        return F
    end
end
