# potri2
An efficient version of LAPACK's `potri` using the algorithm in "Matrix Inversion Using Cholesky Decomposition," by Aravindh Krishnamoorthy and Deepak Menon, [arXiv:1111.4144](https://arxiv.org/abs/1111.4144).

**Note:** This function is only partially implemented.

- For $N < 32,$ the scalar version in `dpotri2s.f90` is used, which is partially optimized.
- For $N \geq 32,$ the block version in `dpotri2b.f90` is used, which is not yet optimized or parallelized.
- Only the IEEE double precision version is currently being implemented. Once complete, the single precision and complex-valued variants will be implemented.

## IEEE double precision

Linux (x86_64-linux-gnu) on 8×11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz with MKL single thread

| N | MKL/U (ns) | MKL/L (ns) | Fortran/U (ns) | Fortran/L (ns) |
| :--- | :--- | :--- | :--- | :--- |
| 2 | 2.1e+02 | 2e+02 | 35 (0.17) | 38 (0.19) |
| 4 | 3.5e+02 | 4e+02 | 66 (0.19) | 74 (0.18) |
| 8 | 7.4e+02 | 7.6e+02 | 2e+02 (0.27) | 2.4e+02 (0.32) |
| 16 | 2e+03 | 1.9e+03 | 1.1e+03 (0.56) | 1.2e+03 (0.62) |
| 32 | 4.7e+03 | 4.6e+03 | 6.5e+03 (1.4) | 5.4e+03 (1.2) |
| 64 | 1.7e+04 | 1.8e+04 | 3.5e+04 (2.1) | 3.2e+04 (1.8) |
| 128 | 6e+04 | 6.6e+04 | 1.8e+05 (2.9) | 1.6e+05 (2.5) |
| 256 | 3.3e+05 | 3.5e+05 | 1.1e+06 (3.3) | 9.7e+05 (2.7) |
| 512 | 2.2e+06 | 2.3e+06 | 7e+06 (3.2) | 6.6e+06 (2.9) |
| 1024 | 1.8e+07 | 1.7e+07 | 5.4e+07 (3) | 5.3e+07 (3) |

