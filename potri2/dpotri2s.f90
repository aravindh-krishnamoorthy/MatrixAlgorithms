!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2S(UPLO, N, A, LDA, INFO)
   IMPLICIT NONE

   CHARACTER UPLO
   INTEGER INFO, LDA, N
   DOUBLE PRECISION A(LDA, *)

   INTEGER I, J, K
   DOUBLE PRECISION ONE, ZERO
   PARAMETER(ONE=1.0, ZERO=0.0)

   IF (UPLO .EQ. 'U') THEN
      DO I = 1, N
         A(I, I) = 1/A(I, I)
         DO CONCURRENT(J=I + 1:N)
            A(I, J) = A(I, J)*A(I, I)
            A(J, I) = 0
         END DO
         A(I, I) = A(I, I)*A(I, I)
      END DO
      DO J = N, 1, -1
         DO K = N, J + 1, -1
            DO CONCURRENT(I=1:J)
               A(J, I) = A(J, I) - A(I, K)*A(K, J)
            END DO
         END DO
         DO K = J, 1, -1
            DO CONCURRENT(I=1:K - 1)
               A(J, I) = A(J, I) - A(I, K)*A(J, K)
            END DO
         END DO
      END DO
      DO CONCURRENT(I=1:N)
         DO CONCURRENT(J=I + 1:N)
            A(I, J) = A(J, I)
         END DO
      END DO
   ELSE ! UPLO.EQ.'L'
      DO J = 1, N
         A(J, J) = 1/A(J, J)
         DO CONCURRENT(I=J + 1:N)
            A(I, J) = A(I, J)*A(J, J)
            A(J, I) = 0
         END DO
         A(J, J) = A(J, J)*A(J, J)
      END DO
      DO J = N, 1, -1
         DO K = N, J + 1, -1
            DO CONCURRENT(I=1:J)
               A(I, J) = A(I, J) - A(K, I)*A(J, K)
            END DO
         END DO
         DO K = J, 1, -1
            DO CONCURRENT(I=1:K - 1)
               A(I, J) = A(I, J) - A(K, I)*A(K, J)
            END DO
         END DO
      END DO
      DO CONCURRENT(I=1:N)
         DO CONCURRENT(J=1:I - 1)
            A(I, J) = A(J, I)
         END DO
      END DO
   END IF
   INFO = 0
   RETURN
END
