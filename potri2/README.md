# potri2
An efficient version of LAPACK's `potri` using the algorithm in "Matrix Inversion Using Cholesky Decomposition," by Aravindh Krishnamoorthy and Deepak Menon, [arXiv:1111.4144](https://arxiv.org/abs/1111.4144).

**Note:** This function is only partially implemented.

- For $N < 32,$ the scalar version in `dpotri2s.f90` is used, which is partially optimized.
- For $N \geq 32,$ the block version in `dpotri2b.f90` is used, which, while being nearly optimized. For parallelization, see PR https://github.com/aravindh-krishnamoorthy/MatrixAlgorithms/pull/1
- Only the IEEE double precision version is currently being implemented. Once complete, the single precision and complex-valued variants will be implemented.

## IEEE double precision

Linux (x86_64-linux-gnu) on 8×11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz with MKL single thread

| N | MKL/U (ns) | MKL/L (ns) | Fortran/U (ns) | Fortran/L (ns) |
| :--- | :--- | :--- | :--- | :--- |
| 2 | 2.5e+02 | 2.5e+02 | 36 (0.14) | 38 (0.15) |
| 4 | 4.1e+02 | 4.2e+02 | 65 (0.16) | 75 (0.18) |
| 8 | 7.2e+02 | 7.8e+02 | 2.1e+02 (0.29) | 2.3e+02 (0.3) |
| 16 | 2e+03 | 2.1e+03 | 1.1e+03 (0.55) | 1.2e+03 (0.56) |
| 32 | 5.2e+03 | 4.6e+03 | 6.1e+03 (1.2) | 5.5e+03 (1.2) |
| 64 | 1.7e+04 | 1.7e+04 | 3.3e+04 (2) | 3.1e+04 (1.8) |
| 100 | 3.9e+04 | 4.1e+04 | 1e+05 (2.6) | 9.9e+04 (2.4) |
| 128 | 6.4e+04 | 6.6e+04 | 1.8e+05 (2.9) | 1.7e+05 (2.6) |
| 256 | 3.4e+05 | 3.6e+05 | 1.1e+06 (3.3) | 1e+06 (2.8) |
| 500 | 2.1e+06 | 2.2e+06 | 5.2e+06 (2.5) | 5.4e+06 (2.5) |
| 512 | 2.6e+06 | 2.6e+06 | 7.7e+06 (2.9) | 6e+06 (2.3) |
| 1024 | 2.2e+07 | 1.8e+07 | 3.7e+07 (1.7) | 3.6e+07 (2) |


