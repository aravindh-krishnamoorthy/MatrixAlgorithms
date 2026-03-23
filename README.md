# MatrixAlgorithms
_Staging area for matrix algebra algorithms in **M**ATLAB, **J**ulia, **F**ortran, **R**ust, and **C**/C++_

<div align="center">

  | Directory | Language | Description | Target | Development Stage |
  |---|---|---|---|---|
  | [potri2](https://github.com/aravindh-krishnamoorthy/MatrixAlgorithms/tree/main/potri2) | J | Efficient implementation of LAPACK's `potri` function | LAPACK | Partial |
  | [trsr](https://github.com/aravindh-krishnamoorthy/MatrixAlgorithms/tree/main/trsr) | J-F | Efficient square root of matrices |  | Partial |
  | [matrixfunctions](https://github.com/aravindh-krishnamoorthy/MatrixAlgorithms/tree/main/matrixfunctions) | J | Matrix functions and their Fréchet derivatives |  | Partial |
  | [ordschur](https://github.com/aravindh-krishnamoorthy/MatrixAlgorithms/tree/main/ordschur) | J | Efficient reordering of Schur triangular matrices (inspired by https://github.com/JuliaLang/julia/issues/53239) |  | Not yet started |
  
</div>

---

## Julia Usage
This package can be loaded under Julia as follows:
```julia
julia> include("/<path>/MatrixAlgorithms/src/MatrixAlgorithms.jl")
```
