        MODULE DOT
        CONTAINS
        SUBROUTINE MAT_VEC_DOT(A,B,i,a_min,a_max,b_min,b_max,out)
                REAL*8, DIMENSION(:,:), INTENT(IN) :: A
                REAL*8, DIMENSION(:), INTENT(IN) :: B
                INTEGER, INTENT(IN) :: i, a_min, a_max, b_min, b_max
                REAL*8, INTENT(OUT) :: out
                
                INTEGER :: j

                out = 0
                
                DO j=1,a_max-a_min+1
                        out = out + A(i,a_min-1+j) * B(b_min-1+j)
                END DO

        END SUBROUTINE
        END MODULE DOT

        SUBROUTINE BANDED_MULT_IN_PLACE(A, B, N_DIM, M_DIM, BANDWIDTH)
                USE DOT
                INTEGER, INTENT(IN) :: N_DIM, M_DIM, BANDWIDTH
                REAL*8, DIMENSION(N_DIM,M_DIM), INTENT(IN) :: A
                REAL*8, DIMENSION(N_DIM), INTENT(INOUT) :: B

                REAL*8, DIMENSION(BANDWIDTH) :: storage
                INTEGER :: i, j, a_min, a_max, b_min, b_max, s
                REAL*8 :: temp
                storage = 0

        if (BANDWIDTH > 0) THEN
        DO i=1,N_DIM
                CALL MIN_BOUND(i, BANDWIDTH, a_min)
                a_min = MAX(1, a_min)

                CALL MAX_BOUND(i, BANDWIDTH, N_DIM, a_max)
                a_max = MIN(M_DIM, a_max)

                b_min = MAX(1, i-BANDWIDTH)
                b_max = MIN(N_DIM, i+BANDWIDTH)

                s = MOD(i-1, BANDWIDTH)+1

                CALL MAT_VEC_DOT(A,B,i,a_min,a_max,b_min,b_max,temp)

        IF (i > BANDWIDTH) THEN
                B(i - BANDWIDTH) = storage(s)

        END IF
                storage(s) = temp
                
        END DO

        DO i=1,BANDWIDTH
                j = MOD(BANDWIDTH+s-i-1, BANDWIDTH)+1
                B(N_DIM-(BANDWIDTH)+i) = storage(j)
        END DO
        ELSE 
        DO i=1,N_DIM
                B(i) = A(i,1) * B(i)
        END DO
        END IF

        END SUBROUTINE
        
        SUBROUTINE BANDED_MULT_OUTPUT(A,B,N_DIM,M_DIM,BANDWIDTH,OUT)
                USE DOT
                INTEGER, INTENT(IN) :: N_DIM, M_DIM, BANDWIDTH
                REAL*8, DIMENSION(N_DIM,M_DIM), INTENT(IN) :: A
                REAL*8, DIMENSION(N_DIM), INTENT(IN) :: B
                REAL*8, DIMENSION(N_DIM), INTENT(OUT) :: OUT

                REAL*8, DIMENSION(BANDWIDTH) :: storage
                INTEGER :: i, j, a_min, a_max, b_min, b_max, s
                storage = 0

        if (BANDWIDTH > 0) THEN
        DO i=1,N_DIM
                CALL MIN_BOUND(i, BANDWIDTH, a_min)
                a_min = MAX(1, a_min)

                CALL MAX_BOUND(i, BANDWIDTH, N_DIM, a_max)
                a_max = MIN(M_DIM, a_max)

                b_min = MAX(1, i-BANDWIDTH)
                b_max = MIN(N_DIM, i+BANDWIDTH)

                CALL MAT_VEC_DOT(A,B,i,a_min,a_max,b_min,b_max,OUT(i))

        END DO
        ELSE 
        DO i=1,N_DIM
                OUT(i) = A(i,1) * B(i)
        END DO
        END IF
        END SUBROUTINE

        SUBROUTINE MIN_BOUND(i, BANDWIDTH, j_min)
                INTEGER, INTENT(IN) :: i, BANDWIDTH
                INTEGER, INTENT(OUT) :: j_min

                j_min = 2 + BANDWIDTH - i

        END SUBROUTINE

        SUBROUTINE MAX_BOUND(i, BANDWIDTH, N_DIM, j_max)
                INTEGER, INTENT(IN) :: i, BANDWIDTH, N_DIM
                INTEGER, INTENT(OUT) :: j_max

                j_max = BANDWIDTH + N_DIM + 1 - i
        END SUBROUTINE
