        PROGRAM lapack_test
                INTEGER, PARAMETER :: m=7, b_width=2
                REAL*8, DIMENSION(3*b_width+1, m) :: A
                REAL*8, DIMENSION(m) :: b
                INTEGER, DIMENSION(m) :: IPIV

                INTEGER :: i, INFO
                
                A = 0
                A(3, :) = 1
                A(4, :) = 3
                A(5, :) = 1
                A(6, :) = 2
                A(7, :) = -1

                DO i=1,m
                        b(i) = i
                END DO
                
                CALL DGBTRF(m, m, b_width, b_width, A, 3*b_width+1, IPIV, INFO)
                PRINT *, INFO

                CALL DGBTRS('N', m, b_width, b_width, 1, A, 3*b_width+1, IPIV, b, m, INFO)
                print *, INFO

                DO i=1,m
                        PRINT *, b(i)
                END DO

        END PROGRAM
