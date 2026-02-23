using Random
using LinearAlgebra
include("../src/MatrixAlgorithms.jl")

rng = MersenneTwister(555);
N = 1024
A = randn(rng,N,N)
A = A*A'

T = Matrix(cholesky(A).L) ;

T1 = LAPACK.potri!('L', copy(T))
T1 = tril(T1)
display(T1)
T2 = MatrixAlgorithms.potri2_parallel!('L',copy(T))
T2 = triu(T2)'
display(T2)
display(T1-T2)
norm(T1-T2)
