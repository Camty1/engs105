        program hello
        
        REAL, DIMENSION(1200, 6) :: A
        REAL, DIMENSION(6, 1200) :: A_prime
        REAL, DIMENSION(6, 6) :: B
        REAL, DIMENSION(1200, 1200) :: C
        INTEGER :: i, j
        
        OPEN(unit = 1, file = "hw0.dat")

        DO i=1,2
                READ(1, *)
        END DO
        
        DO j=1,1200
                READ(1, *) A(j,:)
        END DO
        
        B = MATMUL(TRANSPOSE(A), A)

        PRINT *, B

        end program hello
