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

   INTEGER NB, IB, JB
   INTEGER I, J
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

   DO J = N, 1, -NB
      JB = MIN(J, NB)
      CALL DPOTRI2BD(UPLO, N, A, LDA, INFO, J, JB)
      DO I = J - JB, 1, -NB
         IB = MIN(I, NB)
         CALL DPOTRI2BO(UPLO, N, A, LDA, INFO, J, JB, I, IB)
      END DO
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
   DOUBLE PRECISION DDOT, WORK(JB, JB)
   EXTERNAL DDOT, DGEMV, DTRSV, DGEMM

   IF (J < N) THEN
      CALL DGEMM('T','T', JB, JB, N-J, -1.0D0, A(J+1, J-JB+1), LDA, A(J-JB+1, J+1), LDA, 0.0D0, WORK, JB)

      DO K = J, J-JB+1, -1
         NCOL = K - (J-JB+1) + 1
         DO I = 1, NCOL
            A(J-JB+1 + I - 1, K) = A(J-JB+1 + I - 1, K) + WORK(I, K - (J-JB+1) + 1)
         END DO
      END DO
   END IF

   DO K = J, J-JB+1, -1
      IF (K < J) THEN
         LEN  = J - K
         NCOL = K - (J-JB+1) + 1
         CALL DGEMV('T', LEN, NCOL, -1.0D0, A(K+1, J-JB+1), LDA, A(K, K+1), LDA, 1.0D0, A(J-JB+1, K), 1)
      END IF

      LEN = K - (J-JB+1) + 1
      CALL DTRSV('L','T','U', LEN, A(J-JB+1, J-JB+1), LDA, A(J-JB+1, K), 1)
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
   DOUBLE PRECISION DDOT, WORK(IB,JB)
   EXTERNAL DDOT, DGEMV, DGEMM

   IF (J < N) THEN
      CALL DGEMM('T','T', IB, JB, N-J, -1.0D0, A(J+1, I-IB+1), LDA, A(J-JB+1, J+1), LDA, 0.0D0, WORK, IB)
      DO K = J - JB + 1, J
         DO L = 1, IB
            A(I - IB + 1 + L - 1, K) = A(I - IB + 1 + L - 1, K) + WORK(L, K - (J - JB + 1) + 1)
         END DO
      END DO
   END IF
   DO K = J - JB + 1, J
      IF (K < J) THEN
         LEN = J - K
         CALL DGEMV('T', LEN, IB, -1.0D0, A(K + 1, I - IB + 1), LDA, A(K, K + 1), LDA, 1.0D0, A(I - IB + 1, K), 1)
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
