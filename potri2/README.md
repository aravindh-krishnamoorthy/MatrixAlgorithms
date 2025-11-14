# potri2
An efficient version of LAPACK's `potri` using the algorithm in "Matrix Inversion Using Cholesky Decomposition," by Aravindh Krishnamoorthy and Deepak Menon, [arXiv:1111.4144](https://arxiv.org/abs/1111.4144).

**Note:** This function is only partially implemented.

- For $N < 32,$ the scalar version in `dpotri2s.f90` is used, which is partially optimized.
- For $N \geq 32,$ the block version in `dpotri2b.f90` is used, which, while being nearly optimized, is not yet parallelized.
- Only the IEEE double precision version is currently being implemented. Once complete, the single precision and complex-valued variants will be implemented.

## IEEE double precision

Linux (x86_64-linux-gnu) on 8×11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz with MKL single thread

| N | MKL/U (ns) | MKL/L (ns) | Fortran/U (ns) | Fortran/L (ns) |
| :--- | :--- | :--- | :--- | :--- |
| 2 | 2.6e+02 | 2.5e+02 | 36 (0.14) | 39 (0.15) |
| 4 | 4.3e+02 | 4.2e+02 | 66 (0.15) | 74 (0.18) |
| 8 | 8.8e+02 | 7.9e+02 | 2.2e+02 (0.25) | 2.5e+02 (0.32) |
| 16 | 2e+03 | 2e+03 | 1.2e+03 (0.59) | 1.3e+03 (0.63) |
| 32 | 5.1e+03 | 5e+03 | 6.3e+03 (1.2) | 6e+03 (1.2) |
| 64 | 1.9e+04 | 1.7e+04 | 3.5e+04 (1.8) | 3.3e+04 (1.9) |
| 100 | 4.1e+04 | 4.3e+04 | 1e+05 (2.5) | 9.6e+04 (2.2) |
| 128 | 6.5e+04 | 6.8e+04 | 1.8e+05 (2.8) | 1.7e+05 (2.5) |
| 256 | 3.5e+05 | 3.7e+05 | 1.2e+06 (3.3) | 9.9e+05 (2.7) |
| 500 | 2.1e+06 | 2.2e+06 | 5.3e+06 (2.5) | 5.3e+06 (2.4) |
| 512 | 2.3e+06 | 2.4e+06 | 6.8e+06 (3) | 6e+06 (2.5) |
| 1024 | 1.8e+07 | 1.9e+07 | 4.9e+07 (2.7) | 4.6e+07 (2.5) |


