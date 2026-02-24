################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
################################################################################

using LinearAlgebra
using BenchmarkTools
using Printf
# using Plots
# unicodeplots()

include("../src/MatrixAlgorithms.jl")
# BLAS.set_num_threads(1)

if length(ARGS) > 0
    MS = let expr = Meta.parse(ARGS[1])
        @assert expr.head == :vect
        Int.(expr.args)
    end
else
    MS = [2, 4, 8, 16, 32, 64, 100, 128, 256, 500, 512, 1024]
end

# Real arithmetic
RU = zeros(length(MS))
RL = zeros(length(MS))
PU = zeros(length(MS))
PL = zeros(length(MS))
FU = zeros(length(MS))
FL = zeros(length(MS))
MKLU = zeros(length(MS))
MKLL = zeros(length(MS))
for i in 1:length(MS)
    N = MS[i] 
    println("N=$N...")

    # Real
    X = rand(N, N)
    X = X*X'
    U = Matrix(cholesky(X).U)
    L = Matrix(cholesky(X).L)
    X0 = LAPACK.potri!('U', copy(U))
    X1 = MatrixAlgorithms.potri2!('U', copy(U))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2!('L', copy(L))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2_parallel!('U', copy(U))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2_parallel!('L', copy(L))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2_blocked!('U', copy(U))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2_blocked!('L', copy(L))
    display(norm(triu(X1) - X0))
    # Timing
    b = @benchmark LAPACK.potri!('U', copy($U)) ;
    MKLU[i] = mean(b).time
    b = @benchmark LAPACK.potri!('L', copy($L)) ;
    MKLL[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2!('U', copy($U)) ;
    RU[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2!('L', copy($L)) ;
    RL[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2_parallel!('U', copy($U)) ;
    PU[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2_parallel!('L', copy($L)) ;
    PL[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2_blocked!('U', copy($U)) ;
    FU[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2_blocked!('L', copy($L)) ;
    FL[i] = mean(b).time
end

println("| N | REF/U (ns) | JREF/U (ns) | PAR/U (ns) | BLO/U (ns) | REF/L (ns) | JREF/L (ns) | PAR/L (ns) | BLO/L (ns)")
println("| :--- | :--- | :--- | :--- | :--- |")
for i = 1:length(MS)
    Printf.@printf("| %d | %.2g | %.2g (%.2g) | %.2g (%.2g) | %.2g (%.2g) | %.2g | %.2g (%.2g) | %.2g (%.2g) | %.2g (%.2g) |\n", MS[i], MKLU[i], RU[i], RU[i]/MKLU[i], PU[i], PU[i]/MKLU[i], FU[i], FU[i]/MKLU[i], MKLL[i], RL[i], RL[i]/MKLL[i], PL[i], PL[i]/MKLL[i], FL[i], FL[i]/MKLL[i])
end

# Complex arithmetic
RU = zeros(length(MS))
RL = zeros(length(MS))
PU = zeros(length(MS))
PL = zeros(length(MS))
FU = zeros(length(MS))
FL = zeros(length(MS))
MKLU = zeros(length(MS))
MKLL = zeros(length(MS))
for i in 1:length(MS)
    N = MS[i] 
    println("N=$N...")

    # Real
    X = complex.(rand(N, N), rand(N, N))
    X = X*X'
    U = Matrix(cholesky(X).U)
    L = Matrix(cholesky(X).L)
    X0 = LAPACK.potri!('U', copy(U))
    X1 = MatrixAlgorithms.potri2!('U', copy(U))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2!('L', copy(L))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2_parallel!('U', copy(U))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2_parallel!('L', copy(L))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2_blocked!('U', copy(U))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2_blocked!('L', copy(L))
    display(norm(triu(X1) - X0))
    # Timing
    b = @benchmark LAPACK.potri!('U', copy($U)) ;
    MKLU[i] = mean(b).time
    b = @benchmark LAPACK.potri!('L', copy($L)) ;
    MKLL[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2!('U', copy($U)) ;
    RU[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2!('L', copy($L)) ;
    RL[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2_parallel!('U', copy($U)) ;
    PU[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2_parallel!('L', copy($L)) ;
    PL[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2_blocked!('U', copy($U)) ;
    FU[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.potri2_blocked!('L', copy($L)) ;
    FL[i] = mean(b).time
end

println("| N | REF/U (ns) | JREF/U (ns) | PAR/U (ns) | BLO/U (ns) | REF/L (ns) | JREF/L (ns) | PAR/L (ns) | BLO/L (ns)")
println("| :--- | :--- | :--- | :--- | :--- |")
for i = 1:length(MS)
    Printf.@printf("| %d | %.2g | %.2g (%.2g) | %.2g (%.2g) | %.2g (%.2g) | %.2g | %.2g (%.2g) | %.2g (%.2g) | %.2g (%.2g) |\n", MS[i], MKLU[i], RU[i], RU[i]/MKLU[i], PU[i], PU[i]/MKLU[i], FU[i], FU[i]/MKLU[i], MKLL[i], RL[i], RL[i]/MKLL[i], PL[i], PL[i]/MKLL[i], FL[i], FL[i]/MKLL[i])
end
