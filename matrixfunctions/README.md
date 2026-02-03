<h1 align="center">MatrixFunctions.jl</h1>
<p align="center">
A Julia package for computing scalar and matrix-valued functions of matrix variables and their Fréchet derivatives.<br>
Released under the MIT license, see LICENSE file for details.<br>
Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.<br>
</p>

_**Note:** This package is still under development. The features with check marks below are available. Those without are to be developed._

## Table of Contents
- [Preparation](#preparation)
- [Implementation](#implementation)
- [Usage](#usage)
- [References](#references)

## Preparation
### The Book [1]: Reading progress
- [X] Chapter 1: Theory of Matrix Functions
- [X] Chapter 4: Techniques for General Functions
- [X] Chapter 9: Schur-Parlett Algorithm
- [ ] Chapter 13: Function of Matrix Times Vector: $f(\boldsymbol{A})\boldsymbol{b}$
- [ ] Chapter 2: Applications
- [ ] Chapter 3: Conditioning
- [ ] Others

## Implementation
### Matrix-valued functions of matrix variables
- Matrix function `fm(f, X)`
    - [X]  `fm_schur_parlett_recurrence(f, X)` File: [schur_parlett.jl](https://github.com/aravindh-krishnamoorthy/MatrixFunctions.jl/blob/main/src/schur_parlett.jl)
    - [ ] ⚡ `fm_schur_parlett_block(f, X)` File: [schur_parlett.jl](https://github.com/aravindh-krishnamoorthy/MatrixFunctions.jl/blob/main/src/schur_parlett.jl)
    - [ ] ...
    - [ ] Tests
    - [ ] Code optimisation
   
### Fréchet derivatives of matrix functions
- Fréchet derivative of a matrix function `fm_frechet(f, X, E)`
    - [ ] `fm_frechet_[...](f, X, E)` 
    - [ ] ...
    - [ ] Tests
    - [ ] Code optimisation

### Function of matrix times vector
- Matrix times vector $f(\boldsymbol{X})\boldsymbol{y}$ `fmv(f, X, y)`
    - [ ] `fmv_[..](f, X, y)` 
    - [ ] ...
    - [ ] Tests
    - [ ] Code optimisation

### Matrix polynomials
- Matrix polynomial `pm(p, X)`
    - [ ] `pm_[...](f, X, E)`  
    - [ ] ...
    - [ ] Tests
    - [ ] Code optimisation
- Fréchet derivative of a matrix polynomial `pm_frechet(f, X, E)`
    - [ ] `pm_frechet_[...](f, X, E)`  
    - [ ] ...
    - [ ] Tests
    - [ ] Code optimisation

## References
| | |
| --- | --- |
| [1] | Higham, Nicholas J. Functions of matrices: theory and computation. Society for Industrial and Applied Mathematics, 2008. |

<!---
## Notes
- From https://github.com/JuliaLang/julia/discussions/43982#discussioncomment-6678802
  -  The block method given above (Eq 3.16)
  -  Daleckii-Krein using the eigendecomposition for normal matrices (Corollary 3.12)
  -  Using the Schur decomposition, if you know the Frechet derivative of f(T) for triangular T (Problem 3.2)
--->
