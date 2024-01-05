        program first
        
        REAL, DIMENSION(1200, 6) :: A
        REAL, DIMENSION(6, 6) :: B, B_sym, B_sym_norm
        REAL, DIMENSION(4) :: B_sub
        REAL, DIMENSION(1200, 1200) :: C, C_sym, C_sym_norm
        REAL, DIMENSION(1198) :: C_sub
        REAL :: mean_square, rms_1, rms_2
        INTEGER :: i, j
        
        OPEN(unit=1, file="hw0.dat")

        DO i=1,2
                READ(1, *)
        END DO
        
        DO j=1,1200
                READ(1, *) A(j,:)
        END DO
        
        CLOSE(1)
        OPEN(unit=1, file="B.dat")

        PRINT *, 'MAXIMUM VALUE OF A: ', MAXVAL(A)
        PRINT *, 'MINIMUM VALUE OF A: ', MINVAL(A)
        
        mean_square = 0

        DO i=1,1200
                DO j=1,6
                        mean_square = mean_square + A(i,j) * A(i,j)
                END DO
        END DO
        
        mean_square = mean_square / (1200 * 6)

        PRINT *, 'MEAN SQUARE VALUE OF A: ', mean_square 
        
        B = MATMUL(TRANSPOSE(A), A)
        C = MATMUL(A, TRANSPOSE(A))

        WRITE(1, *) B
        CLOSE(1)

        OPEN(unit=1, file="C.dat")
        WRITE(1, *) C
        CLOSE(1)
        
        rms_1 = 0
        rms_2 = 0

        DO i=1,6
                DO j=1,6
                        B_sym(i,j) = ABS(B(i,j)-B(j,i))
                        rms_1 = rms_1+B_sym(i,j)*B_sym(i,j)
                        B_sym_norm(i,j) = B_sym(i,j)/(2*(B(i,j)+B(j,i)))
                        rms_2 = rms_2+B_sym_norm(i,j)+B_sym_norm(i,j)
                END DO
        END DO
        
        rms_1 = SQRT(rms_1 / (6*6))
        rms_2 = SQRT(rms_2 / (6*6))
        
        PRINT *, ''
        PRINT *, 'SYMMETRY CHECK OF B'
        PRINT *, 'MEAN VALUE: ', SUM(B_sym)/(6*6)
        PRINT *, 'RMS VALUE: ', rms_1
        PRINT *, 'MAXIMUM VALUE: ', MAXVAL(B_sym)

        PRINT *, ''
        PRINT *, 'NORMALIZED SYMMETRY CHECK OF B'
        PRINT *, 'MEAN VALUE: ', SUM(B_sym_norm)/(6*6)
        PRINT *, 'RMS VALUE: ', rms_2 
        PRINT *, 'MAXIMUM VALUE: ', MAXVAL(B_sym_norm)
        
        PRINT *, ''
        PRINT *, 'THIRD ROW OF B WRITTEN TO B3.dat'

        OPEN(unit=1, file="B3.dat")
        WRITE(1, *) B(3,:)
        CLOSE(1)
        
        DO i=1,4
                B_sub(i) = B(i+2,i)
        END DO

        PRINT *, ''
        PRINT *, 'SECOND SUBDIAGONAL OF B WRITTEN TO B_sub.dat'

        OPEN(unit=1, file="B_sub.dat")
        WRITE(1, *) B_sub
        CLOSE(1)
        
        rms_1 = 0
        rms_2 = 0

        DO i=1,1200
                DO j=1,1200
                        C_sym(i,j) = ABS(C(i,j)-C(j,i))
                        rms_1 = rms_1+C_sym(i,j)*C_sym(i,j)
                        C_sym_norm(i,j) = C_sym(i,j)/(2*(C(i,j)+C(j,i)))
                        rms_2 = rms_2+C_sym_norm(i,j)+C_sym_norm(i,j)
                END DO
        END DO
        
        rms_1 = SQRT(rms_1 / (1200*1200))
        rms_2 = SQRT(rms_2 / (1200*1200))
        
        PRINT *, ''
        PRINT *, 'SYMMETRY CHECK OF C'
        PRINT *, 'MEAN VALUE: ', SUM(C_sym)/(1200*1200)
        PRINT *, 'RMS VALUE: ', rms_1
        PRINT *, 'MAXIMUM VALUE: ', MAXVAL(C_sym)

        PRINT *, ''
        PRINT *, 'NORMALIZED SYMMETRY CHECK OF C'
        PRINT *, 'MEAN VALUE: ', SUM(C_sym_norm)/(1200*1200)
        PRINT *, 'RMS VALUE: ', rms_2 
        PRINT *, 'MAXIMUM VALUE: ', MAXVAL(C_sym_norm)

        PRINT *, ''
        PRINT *, 'THIRD ROW OF C WRITTEN TO C3.dat'
        OPEN(unit=1, file="C3.dat")
        WRITE(1, *) C(3,:)
        CLOSE(1)
        
        DO i=1,1198
                C_sub(i) = C(i+2,i)
        END DO

        PRINT *, ''
        PRINT *, 'SECOND SUBDIAGONAL OF C WRITTEN TO C_sub.dat'
        OPEN(unit=1, file="C_sub.dat")
        WRITE(1, *) C_sub
        CLOSE(1)

        end program first
