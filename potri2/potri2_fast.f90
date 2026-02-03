! potri2_fast.f90
!
! Invert an SPD/HPD matrix given its Cholesky factor (triangular) without BLAS/LAPACK.
! Algorithm (faster + cache-friendlier than the Julia reference):
!   1) Invert the triangular factor in-place (forward substitution, lower-triangular).
!   2) Form A^{-1} = inv(L)^H * inv(L) (compute into the *upper* triangle to avoid overwrites),
!      then Hermitian-symmetrize.
!
! For best performance: provide the *lower* Cholesky factor (uplo='L').
! If uplo='U', this code copies U^H into a temporary lower-triangular matrix and proceeds.
!
! By ChatGPT 5.2 Thinking on 03.02.2026 13:30
! Not yet reviewed

module potri2_fast
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: potri2

  interface potri2
     module procedure potri2_real64
     module procedure potri2_cplx64
  end interface

contains

  subroutine potri2_real64(uplo, A, info)
    character(len=1), intent(in)    :: uplo
    real(real64),     intent(inout) :: A(:,:)
    integer,          intent(out)   :: info
    integer :: n
    real(real64), allocatable :: W(:,:)

    info = 0
    n = size(A,1)
    if (size(A,2) /= n) then
      info = -2; return
    end if

    select case (uplo)
    case ('L','l')
      call inv_lower_real64(A, info)
      if (info /= 0) return
      call gram_upper_real64(A)
      call symm_from_upper_real64(A)

    case ('U','u')
      allocate(W(n,n))
      call copy_Ut_to_L_real64(A, W)          ! W := transpose(U) in lower storage
      call inv_lower_real64(W, info)
      if (info /= 0) return
      call gram_upper_real64(W)
      call symm_from_upper_real64(W)
      A = W
      deallocate(W)

    case default
      info = -1
    end select
  end subroutine potri2_real64


  subroutine potri2_cplx64(uplo, A, info)
    character(len=1), intent(in)      :: uplo
    complex(real64),  intent(inout)   :: A(:,:)
    integer,          intent(out)     :: info
    integer :: n
    complex(real64), allocatable :: W(:,:)

    info = 0
    n = size(A,1)
    if (size(A,2) /= n) then
      info = -2; return
    end if

    select case (uplo)
    case ('L','l')
      call inv_lower_cplx64(A, info)
      if (info /= 0) return
      call gram_upper_cplx64(A)
      call herm_from_upper_cplx64(A)

    case ('U','u')
      allocate(W(n,n))
      call copy_Uh_to_L_cplx64(A, W)          ! W := U^H in lower storage (conjugate-transpose)
      call inv_lower_cplx64(W, info)
      if (info /= 0) return
      call gram_upper_cplx64(W)
      call herm_from_upper_cplx64(W)
      A = W
      deallocate(W)

    case default
      info = -1
    end select
  end subroutine potri2_cplx64


  !-------------------------
  ! Lower-triangular inversion in-place (non-unit diagonal).
  ! A on entry: lower-triangular Cholesky factor L (upper part ignored).
  ! A on exit : inv(L) in lower triangle. Upper part is zeroed.
  !
  ! Column-by-column forward substitution for inv(L):
  !   inv(L)(i,j) = -( sum_{k=j}^{i-1} L(i,k)*inv(L)(k,j) ) / L(i,i),  i>j
  !   inv(L)(j,j) = 1 / L(j,j)
  !-------------------------

  subroutine inv_lower_real64(A, info)
    real(real64), intent(inout) :: A(:,:)
    integer,      intent(out)   :: info
    integer :: n, i, j, k
    real(real64) :: diag, s

    info = 0
    n = size(A,1)

    ! Invert L in-place (lower triangle)
    do j = 1, n
      diag = A(j,j)
      if (diag == 0.0_real64) then
        info = j; return
      end if
      A(j,j) = 1.0_real64 / diag

      do i = j+1, n
        s = 0.0_real64
        !$omp simd reduction(+:s)
        do k = j, i-1
          s = s + A(i,k) * A(k,j)
        end do
        A(i,j) = -s / A(i,i)
      end do
    end do

    ! Zero upper triangle
    do j = 2, n
      do i = 1, j-1
        A(i,j) = 0.0_real64
      end do
    end do
  end subroutine inv_lower_real64


  subroutine inv_lower_cplx64(A, info)
    complex(real64), intent(inout) :: A(:,:)
    integer,         intent(out)   :: info
    integer :: n, i, j, k
    complex(real64) :: diag, s

    info = 0
    n = size(A,1)

    do j = 1, n
      diag = A(j,j)
      if (diag == (0.0_real64,0.0_real64)) then
        info = j; return
      end if
      A(j,j) = (1.0_real64,0.0_real64) / diag

      do i = j+1, n
        s = (0.0_real64,0.0_real64)
        !$omp simd reduction(+:s)
        do k = j, i-1
          s = s + A(i,k) * A(k,j)
        end do
        A(i,j) = -s / A(i,i)
      end do
    end do

    do j = 2, n
      do i = 1, j-1
        A(i,j) = (0.0_real64,0.0_real64)
      end do
    end do
  end subroutine inv_lower_cplx64


  !-------------------------
  ! Form upper triangle of inv(A) without overwriting inv(L) (stored in lower).
  !
  ! For real:    C = inv(L)^T * inv(L)
  ! For complex: C = inv(L)^H * inv(L)
  !
  ! Compute only i<=j:
  !   C(i,j) = sum_{k=j..n} conj(inv(L)(k,i)) * inv(L)(k,j)
  !
  ! Writes to A(i,j) (upper) while reading A(k,*) (lower), so no overwrite hazard.
  !-------------------------

  subroutine gram_upper_real64(A)
    real(real64), intent(inout) :: A(:,:)
    integer :: n, i, j, k
    real(real64) :: s

    n = size(A,1)

    !$omp parallel do schedule(static) private(i,k,s)
    do j = 1, n
      do i = 1, j
        s = 0.0_real64
        !$omp simd reduction(+:s)
        do k = j, n
          s = s + A(k,i) * A(k,j)
        end do
        A(i,j) = s
      end do
    end do
    !$omp end parallel do
  end subroutine gram_upper_real64


  subroutine gram_upper_cplx64(A)
    complex(real64), intent(inout) :: A(:,:)
    integer :: n, i, j, k
    complex(real64) :: s

    n = size(A,1)

    !$omp parallel do schedule(static) private(i,k,s)
    do j = 1, n
      do i = 1, j
        s = (0.0_real64,0.0_real64)
        !$omp simd reduction(+:s)
        do k = j, n
          s = s + conjg(A(k,i)) * A(k,j)
        end do
        A(i,j) = s
      end do
    end do
    !$omp end parallel do
  end subroutine gram_upper_cplx64


  !-------------------------
  ! Symmetrize/Hermitianize from computed upper triangle.
  !-------------------------

  subroutine symm_from_upper_real64(A)
    real(real64), intent(inout) :: A(:,:)
    integer :: n, i, j
    n = size(A,1)
    do j = 1, n
      do i = j+1, n
        A(i,j) = A(j,i)
      end do
    end do
  end subroutine symm_from_upper_real64


  subroutine herm_from_upper_cplx64(A)
    complex(real64), intent(inout) :: A(:,:)
    integer :: n, i, j
    n = size(A,1)
    do j = 1, n
      do i = j+1, n
        A(i,j) = conjg(A(j,i))
      end do
      A(j,j) = cmplx(real(A(j,j),kind=real64), 0.0_real64, kind=real64)
    end do
  end subroutine herm_from_upper_cplx64


  !-------------------------
  ! Copy helpers for uplo='U' fast-path conversion to a lower factor.
  !-------------------------

  subroutine copy_Ut_to_L_real64(U, L)
    real(real64), intent(in)  :: U(:,:)
    real(real64), intent(out) :: L(:,:)
    integer :: n, i, j
    n = size(U,1)
    L = 0.0_real64
    do j = 1, n
      do i = 1, j
        L(j,i) = U(i,j)     ! transpose(U) into lower storage
      end do
    end do
  end subroutine copy_Ut_to_L_real64


  subroutine copy_Uh_to_L_cplx64(U, L)
    complex(real64), intent(in)  :: U(:,:)
    complex(real64), intent(out) :: L(:,:)
    integer :: n, i, j
    n = size(U,1)
    L = (0.0_real64,0.0_real64)
    do j = 1, n
      do i = 1, j
        L(j,i) = conjg(U(i,j))   ! U^H into lower storage
      end do
    end do
  end subroutine copy_Uh_to_L_cplx64

end module potri2_fast
