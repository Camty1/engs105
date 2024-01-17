        PROGRAM test
                REAL*8, DIMENSION(6,5) :: A
                REAL*8, DIMENSION(6,1) :: C
                REAL*8, DIMENSION(6) :: B, OUT
                INTEGER :: i
                C = 2
                A(:,1) = -3
                A(:,2) = 2
                A(:,3) = 1
                A(:,4) = -2
                A(:,5) = 3
        DO i=1,6
                B(i) = i
        END DO
                CALL BANDED_MULT_OUTPUT(A, B, 6, 5, 2, OUT)
                CALL BANDED_MULT_IN_PLACE(A, B, 6, 5, 2)
                PRINT *, INT(OUT)
                PRINT *, INT(B)

                CALL BANDED_MULT_IN_PLACE(C, B, 6, 1, 0)
                PRINT *, INT(B)

        END PROGRAM
