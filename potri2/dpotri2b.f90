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
    PARAMETER          ( NB = 128 )

    IF (UPLO.EQ.'U') THEN
        DO CONCURRENT (J = 1:N)
            DO CONCURRENT (I=J+1:N)
                TMP = A(I,J)
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
END SUBROUTINE DPOTRI2B

SUBROUTINE DPOTRI2BD(UPLO, N, A, LDA, INFO, J, JB)
    IMPLICIT           NONE

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N, J, JB
    DOUBLE PRECISION   A( LDA, * )

    INTEGER            I, K, LEN
    DOUBLE PRECISION   DDOT
    EXTERNAL           DDOT

    DO K = J, J-JB+1, -1

        DO I = J-JB+1, K
            LEN = N - K
            IF (LEN > 0) THEN
                A(I,K) = A(I,K) - DDOT(LEN, A(K+1,I), 1, A(K,K+1), LDA)
            END IF
        END DO

        DO I = K-1, J-JB+1, -1
            LEN = K - I
            IF (LEN > 0) THEN
                A(I,K) = A(I,K) - DDOT(LEN, A(I+1,I), 1, A(I+1,K), 1)
            END IF
        END DO

    END DO

    INFO = 0
    RETURN
END SUBROUTINE DPOTRI2BD

SUBROUTINE DPOTRI2BO(UPLO, N, A, LDA, INFO, J, JB, I, IB)
    IMPLICIT           NONE

    CHARACTER          UPLO
    INTEGER            INFO, LDA, N, I, IB, J, JB
    DOUBLE PRECISION   A( LDA, * )

    INTEGER            K, L, LEN
    DOUBLE PRECISION   DDOT
    EXTERNAL           DDOT, DGEMV

    DO K = J-JB+1, J
        IF (K < N) THEN
            CALL DGEMV('T', N-K, IB, -1.0D0, A(K+1, I-IB+1), LDA, &
                       A(K, K+1), LDA, 1.0D0, A(I-IB+1, K), 1)
        END IF
    END DO

    DO L = I, I-IB+1, -1
        DO K = J-JB+1, J
            LEN = K - L
            IF (LEN > 0) THEN
                A(L,K) = A(L,K) - DDOT(LEN, A(L+1,L), 1, A(L+1,K), 1)
            END IF
        END DO
    END DO

    INFO = 0
    RETURN
END SUBROUTINE DPOTRI2BO
