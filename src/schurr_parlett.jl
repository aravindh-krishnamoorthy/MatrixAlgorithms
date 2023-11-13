# Schur-Parlett algorithm based on Algorithm 9.6 in Matrix Functions, Higham.
function fm_schur_parlett(f::Function, X::AbstractMatrix)
    S = schur(X)
    # For diagonalisable matrices, use the straightforward formula
    if isdiag(S.T)
        F = S.Z*f.(S.T)*S.Z'
        return F
    elseif !istriu(S.T)
        # T is only triangular in the complex field
        S = Schur{Complex}(S)
    end
end
