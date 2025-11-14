!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2(UPLO, N, A, LDA, INFO)
   IMPLICIT NONE

   CHARACTER UPLO
   INTEGER INFO, LDA, N
   DOUBLE PRECISION A(LDA, *)

   EXTERNAL DPOTRI2S, DPOTRI2B

   INTEGER NB
   PARAMETER(NB=32)

   IF (NB .LE. 1 .OR. NB .GE. N) THEN
      ! Scalar version
      CALL DPOTRI2S(UPLO, N, A, LDA, INFO)
   ELSE
      ! Block version
      CALL DPOTRI2B(UPLO, N, A, LDA, INFO)
   END IF
   INFO = 0
   RETURN
END
