!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2B(UPLO, N, A, LDA, INFO)
   IMPLICIT NONE

   CHARACTER UPLO
   INTEGER INFO, LDA, N
   DOUBLE PRECISION TMP
   DOUBLE PRECISION A(LDA, *)

   INTEGER NB, IB, JB, BEG
   INTEGER I, J, K, II, JJ
   INTEGER I1, I2, J1, J2

   IF (N.LT.256) THEN
      NB = 32
   ELSE
      NB = 128
   END IF

   IF (UPLO .EQ. 'U') THEN
      DO CONCURRENT(J=1:N)
         DO CONCURRENT(I=J + 1:N)
            TMP = A(I, J)
            A(I, J) = A(J, I)
            A(J, I) = TMP
         END DO
      END DO
   END IF

   DO J = 1, N
      TMP = 1.0D0/A(J, J)
      A(J + 1:N, J) = A(J + 1:N, J)*TMP
      A(J, J) = TMP*TMP
      A(J, J + 1:N) = 0.0D0
   END DO

   ! Anti-diagonal blocks can be run in parallel
   IB = 0
   DO J = N, 1, -NB
      ! Each K index can be run in parallel
      !$omp parallel do default(shared) private(K,II,JJ,IB,JB,INFO) schedule(static)
      DO K = N, J, -NB
         II = K
         IB = MIN(II, NB)
         JJ = J + (N - K)
         JB = MIN(JJ, NB)
         IF (II .EQ. JJ) THEN
            CALL DPOTRI2BD(UPLO, N, A, LDA, INFO, II, IB)
         ELSE IF (II .GT. JJ) THEN
            CALL DPOTRI2BO(UPLO, N, A, LDA, INFO, II, IB, JJ, JB)
         END IF
      END DO
      !$omp end parallel do
   END DO
   BEG = N - NB*((N-1)/NB)
   IF (BEG .EQ. 0) BEG = NB
   DO I = N - NB, 1, -NB
      ! Each K index can be run in parallel
      !$omp parallel do default(shared) private(K,II,JJ,IB,JB,INFO) schedule(static)
      DO K = BEG, I, NB
         II = I - (K - BEG)
         IB = MIN(II, NB)
         JJ = K
         JB = MIN(JJ, NB)
         IF (II .EQ. JJ) THEN
            CALL DPOTRI2BD(UPLO, N, A, LDA, INFO, II, IB)
         ELSE IF (II .GT. JJ) THEN
            CALL DPOTRI2BO(UPLO, N, A, LDA, INFO, II, IB, JJ, JB)
         END IF
      END DO
      !$omp end parallel do
   END DO

   DO CONCURRENT(I=1:N)
      A(I, 1:I - 1) = A(1:I - 1, I)
   END DO

   INFO = 0
   RETURN
END SUBROUTINE DPOTRI2B

SUBROUTINE DPOTRI2BD(UPLO, N, A, LDA, INFO, J, JB)
   IMPLICIT NONE

   CHARACTER UPLO
   INTEGER INFO, LDA, N, J, JB
   DOUBLE PRECISION A(LDA, *)

   INTEGER I, K, LEN, NCOL
   DOUBLE PRECISION DDOT
   EXTERNAL DDOT, DGEMV, DTRSV

   DO K = J, J - JB + 1, -1
      IF (K < N) THEN
         LEN = N - K
         NCOL = K - (J - JB + 1) + 1
         CALL DGEMV('T', LEN, NCOL, -1.0D0, A(K + 1, J - JB + 1), LDA, A(K, K + 1), LDA, 1.0D0, A(J - JB + 1, K), 1)
      END IF
      LEN = K - (J - JB + 1) + 1
      CALL DTRSV('L', 'T', 'U', LEN, A(J - JB + 1, J - JB + 1), LDA, A(J - JB + 1, K), 1)
   END DO

   INFO = 0
   RETURN
END SUBROUTINE DPOTRI2BD

SUBROUTINE DPOTRI2BO(UPLO, N, A, LDA, INFO, J, JB, I, IB)
   IMPLICIT NONE

   CHARACTER UPLO
   INTEGER INFO, LDA, N, I, IB, J, JB
   DOUBLE PRECISION A(LDA, *)

   INTEGER K, L, LEN
   DOUBLE PRECISION DDOT
   EXTERNAL DDOT, DGEMV

   DO K = J - JB + 1, J
      IF (K < N) THEN
         CALL DGEMV('T', N - K, IB, -1.0D0, A(K + 1, I - IB + 1), LDA, A(K, K + 1), LDA, 1.0D0, A(I - IB + 1, K), 1)
      END IF
   END DO

   DO K = J - JB + 1, J
      IF (K - I > 0) THEN
         CALL DGEMV('T', K - I, IB, -1.0D0, A(I + 1, I - IB + 1), LDA, A(I + 1, K), 1, 1.0D0, A(I - IB + 1, K), 1)
      END IF
      LEN = MIN(K, I) - (I - IB + 1) + 1
      CALL DTRSV('L', 'T', 'U', LEN, A(I - IB + 1, I - IB + 1), LDA, A(I - IB + 1, K), 1)
   END DO

   INFO = 0
   RETURN
END SUBROUTINE DPOTRI2BO
