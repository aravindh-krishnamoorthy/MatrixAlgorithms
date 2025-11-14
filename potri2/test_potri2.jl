################################################################################
# This file is a part of the package: MatrixAlgorithms
# Released under the MIT license, see LICENSE file for details.
# Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
################################################################################

using MKL
using LinearAlgebra
using BenchmarkTools
using Printf
# using Plots
# unicodeplots()

include("../src/MatrixAlgorithms.jl")
MKL.set_num_threads(1)

MS = [2, 4, 8, 16, 32, 64, 100, 128, 256, 500, 512, 1024]
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
    X1 = MatrixAlgorithms.potri2_freference!('U', copy(U))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.potri2_freference!('L', copy(L))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.dpotri2!('U', copy(U))
    display(norm(triu(X1) - X0))
    X1 = MatrixAlgorithms.dpotri2!('L', copy(L))
    display(norm(triu(X1) - X0))
    # Timing
    b = @benchmark LAPACK.potri!('U', copy($U)) ;
    MKLU[i] = mean(b).time
    b = @benchmark LAPACK.potri!('L', copy($L)) ;
    MKLL[i] = mean(b).time
    # b = @benchmark MatrixAlgorithms.potri2!('U', copy($U)) ;
    # JU[i] = mean(b).time
    # b = @benchmark MatrixAlgorithms.potri2!('L', copy($L)) ;
    # JL[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.dpotri2!('U', copy($U); rl=false) ;
    FU[i] = mean(b).time
    b = @benchmark MatrixAlgorithms.dpotri2!('L', copy($L); rl=false) ;
    FL[i] = mean(b).time
end

# plt = scatter(MS, MKLU, m=:square, label="MKLU")
# plt = scatter!(plt, MS, MKLL, m=:square, label="MKLL")
# plt = scatter!(plt, MS, FU, m=:cross, label="FU")
# plt = scatter!(plt, MS, FL, m=:cross, label="FL")
# plt = plot!(plt, title="MKL potri vs potri2 in Julia", xlabel="n", ylabel="Time (ns)", xscale=:log2, yscale=:log10, grid=true)

println("| N | MKL/U (ns) | MKL/L (ns) | Fortran/U (ns) | Fortran/L (ns) |")
println("| :--- | :--- | :--- | :--- | :--- |")
for i = 1:length(MS)
    Printf.@printf("| %d | %.2g | %.2g | %.2g (%.2g) | %.2g (%.2g) |\n", MS[i], MKLU[i], MKLL[i], FU[i], FU[i]/MKLU[i], FL[i], FL[i]/MKLL[i])
end
