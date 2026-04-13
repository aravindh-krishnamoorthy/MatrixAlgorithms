using LinearAlgebra
using BenchmarkTools

include("trsr.jl")

for N in [8, 64, 1000]
    # Real
    println("$N: Real/No structure")
    X = randn(N,N)
    @btime sqrt($X)
    @btime sqrtm($X)
    println("$N: Real/Symmetric")
    @btime sqrt($X+$X')
    @btime sqrtm($X+$X')
    println("$N: Real/No structure/PD")
    @btime sqrt($X+$N*I($N))
    @btime sqrtm($X+$N*I($N))
    println("$N: Real/Symmetric/PD")
    @btime sqrt($X+$X'+$N*I($N))
    @btime sqrtm($X+$X'+$N*I($N))
    # Complex
    println("$N: Complex/No structure")
    X = complex.(randn(N,N),randn(N,N))
    @btime sqrt($X)
    @btime sqrtm($X)
    println("$N: Complex/Hermitian")
    @btime sqrt($X+$X')
    @btime sqrtm($X+$X')
    println("$N: Complex/No structure/PD")
    @btime sqrt($X+$N*I($N))
    @btime sqrtm($X+$N*I($N))
    println("$N: Complex/Hermitian/PD")
    @btime sqrt($X+$X'+$N*I($N))
    @btime sqrtm($X+$X'+$N*I($N))
end
