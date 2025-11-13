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
| 2 | 2.5e+02 | 2.4e+02 | 37 (0.15) | 38 (0.15) |
| 4 | 4.3e+02 | 4.1e+02 | 67 (0.15) | 79 (0.19) |
| 8 | 8e+02 | 8.4e+02 | 2.3e+02 (0.28) | 2.5e+02 (0.29) |
| 16 | 2e+03 | 2e+03 | 1.2e+03 (0.6) | 1.3e+03 (0.62) |
| 32 | 5.4e+03 | 5.3e+03 | 6e+03 (1.1) | 5.6e+03 (1.1) |
| 64 | 1.8e+04 | 1.6e+04 | 4.5e+04 (2.5) | 4.7e+04 (2.8) |
| 128 | 6.3e+04 | 6.9e+04 | 2.3e+05 (3.6) | 2.2e+05 (3.3) |
| 256 | 3.3e+05 | 3.6e+05 | 1.4e+06 (4.3) | 1.2e+06 (3.4) |
| 512 | 2.2e+06 | 2.3e+06 | 8.7e+06 (3.9) | 8e+06 (3.4) |
| 1024 | 1.8e+07 | 1.8e+07 | 6.6e+07 (3.7) | 6.5e+07 (3.6) |
