        PROGRAM TEST

        REAL, DIMENSION(10,3) :: A
        REAL, DIMENSION(10) :: b, x
        INTEGER :: i
        A = 1.0
        b = 1.0

        A(1,1) = 0
        A(10,3) = 0
        !A(1,2) = -1
        !A(1,3) = 1
        !A(2,1) = 1
        !A(2,2) = -1
        !A(2,3) = 2
        !A(3,1) = 1
        !A(3,2) = -1
        !A(4,1) = 1
        !A(4,2) = -1
        !A(4,3) = 3
        !A(5,2) = 2
        !A(5,3) = -1

        DO i=1,10
        A(i,1) = 0.0
        A(i,2) = i
        b(i) = i
        PRINT *, A(i,:), b(i)
        END DO
        
        CALL SOLVE(3, A, b, 10, 1, 10, 3)
         
        PRINT *, ""

        DO i=1,10
        PRINT *, b(i)
        END DO

        END PROGRAM
