!###############################################################################
! This file is a part of the package: MatrixAlgorithms
! Released under the MIT license, see LICENSE file for details.
! Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
!###############################################################################

SUBROUTINE DPOTRI2B(UPLO, N, A, LDA, INFO)
    IMPLICIT           NONE

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    DOUBLE PRECISION   TMP
    DOUBLE PRECISION   A( LDA, * )

    INTEGER            NB, IB, JB
    INTEGER            I, J
    PARAMETER          ( NB = 32 )

    IF (UPLO.EQ.'U') THEN
        DO J = 1, N
            DO I = J+1, N
                TMP    = A(I,J)
                A(I,J) = A(J,I)
                A(J,I) = TMP
            END DO
        END DO
    END IF

    DO CONCURRENT (J = 1:N)
        A(J,J) = 1/A(J,J)
        A(J+1:N,J) = A(J+1:N,J)*A(J,J)
        A(J,J+1:N) = 0
        A(J,J) = A(J,J)*A(J,J)
    END DO
    DO J = N, 1-NB, -NB
        JB = MIN(J,NB)
        CALL DPOTRI2BD(UPLO, N, A, LDA, INFO, J, JB)
        DO I = J-JB, 1, -NB
            IB = MIN(I,NB)
            CALL DPOTRI2BO(UPLO, N, A, LDA, INFO, J, JB, I, IB)
        END DO
    END DO
    DO CONCURRENT (I = 1:N)
        A(I,1:I-1) = A(1:I-1,I)
    END DO
    INFO = 0
    RETURN
END

SUBROUTINE DPOTRI2BD(UPLO, N, A, LDA, INFO, J, JB)
    IMPLICIT           NONE

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N, J, JB
    DOUBLE PRECISION   A( LDA, * )
    INTEGER            I, K

    DO K = J, J-JB+1, -1
        DO CONCURRENT (I = J-JB+1:K)
            A(I,K) = A(I,K) - DOT_PRODUCT(A(K+1:N,I), A(K,K+1:N))
        END DO
        DO I = K-1, J-JB+1, -1
            A(I,K) = A(I,K) - DOT_PRODUCT(A(I+1:K,I), A(I+1:K,K))
        END DO
    END DO
    INFO = 0
    RETURN
END SUBROUTINE

SUBROUTINE DPOTRI2BO(UPLO, N, A, LDA, INFO, J, JB, I, IB)
    IMPLICIT           NONE

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N, I, IB, J, JB
    DOUBLE PRECISION   A( LDA, * )
    INTEGER            K, L

    DO CONCURRENT (K = J-JB+1:J)
        A(I-IB+1:I,K) = A(I-IB+1:I,K) - MATMUL(A(K,K+1:N), A(K+1:N,I-IB+1:I))
    END DO
    DO L = I, I-IB+1, -1
        DO CONCURRENT (K = J-JB+1:J)
            A(L,K) = A(L,K) - DOT_PRODUCT(A(L+1:K,L), A(L+1:K,K))
        END DO
    END DO
    INFO = 0
    RETURN
END SUBROUTINE
