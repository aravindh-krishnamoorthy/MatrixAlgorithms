# potri2
An efficient version of LAPACK's `potri` using the algorithm in "Matrix Inversion Using Cholesky Decomposition," by Aravindh Krishnamoorthy and Deepak Menon, [arXiv:1111.4144](https://arxiv.org/abs/1111.4144).

**Note:** This function is only partially implemented.

- For $N < 32,$ the scalar version in `dpotri2s.f90` is used, which is partially optimized.
- For $N \geq 32,$ the block version in `dpotri2b.f90` is used, which is not yet optimized or parallelized.
- Only the IEEE double precision version is currently being implemented. Once complete, the single precision and complex-valued variants will be implemented.

## IEEE double precision

Linux (x86_64-linux-gnu) on 8×11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz for commit [2d93904](https://github.com/aravindh-krishnamoorthy/MatrixAlgorithms/commit/2d93904edbdbd0ce4657899aff0e3d7eb7df8e62)

| N | MKL/U (ns) | MKL/L (ns) | Fortran/U (ns) | Fortran/L (ns) |
| :--- | :--- | :--- | :--- | :--- |
| 2 | 2.4e+02 | 2.2e+02 | 37 (0.15) | 37 (0.17) |
| 4 | 4e+02 | 3.8e+02 | 70 (0.17) | 85 (0.22) |
| 8 | 8.2e+02 | 7.5e+02 | 2.1e+02 (0.25) | 3.4e+02 (0.45) |
| 16 | 2.1e+03 | 1.9e+03 | 1.3e+03 (0.6) | 1.2e+03 (0.67) |
| 32 | 5.1e+03 | 4.9e+03 | 5.9e+03 (1.2) | 5.6e+03 (1.1) |
| 64 | 1.7e+04 | 2.9e+04 | 4.3e+04 (2.5) | 4.5e+04 (1.5) |
| 128 | 1.1e+05 | 1.1e+05 | 2.3e+05 (2.1) | 2.2e+05 (2) |
| 256 | 5.3e+05 | 5.6e+05 | 1.4e+06 (2.6) | 1.2e+06 (2.2) |
| 512 | 1.8e+06 | 1.9e+06 | 8.5e+06 (4.7) | 7.9e+06 (4.2) |
| 1024 | 1.1e+07 | 1.1e+07 | 1e+08 (9.5) | 8.7e+07 (7.6) |
