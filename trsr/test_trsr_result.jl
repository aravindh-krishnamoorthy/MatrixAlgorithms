using LinearAlgebra
using BenchmarkTools

for iter = 1:1000
    println("$iter/1000...")
    for N in [8, 64, 1000]
        # Real trsr!
        X = randn(N,N)
        X1 = sqrtm(X)
        @assert X1*X1 ≈ X
        X1 = sqrtm(X+X')
        @assert X1*X1 ≈ X+X'
        X1 = sqrtm(X+N*I(N))
        @assert X1*X1 ≈ X+N*I(N)
        X1 = sqrtm(X+X'+N*I(N))
        @assert X1*X1 ≈ X+X'+N*I(N)
        # Complex trsr!
        X = complex.(randn(N,N),randn(N,N))
        X1 = sqrtm(X)
        @assert X1*X1 ≈ X
        X1 = sqrtm(X+X')
        @assert X1*X1 ≈ X+X'
        X1 = sqrtm(X+N*I(N))
        @assert X1*X1 ≈ X+N*I(N)
        X1 = sqrtm(X+X'+N*I(N))
        @assert X1*X1 ≈ X+X'+N*I(N)
    end
end
