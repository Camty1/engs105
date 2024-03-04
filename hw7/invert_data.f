        PROGRAM invert_data
                INTEGER, PARAMETER :: NUM_NODE=502, NUM_MEAS=12, BANDWIDTH=23
                
                REAL*8, DIMENSION(3*BANDWIDTH+1, NUM_NODE) :: K
                REAL*8, DIMENSION(NUM_MEAS, NUM_NODE) :: S
                REAL*8, DIMENSION(NUM_NODE, NUM_MEAS) :: S_trans, u_mat
                REAL*8, DIMENSION(NUM_MEAS, NUM_MEAS) :: cov_delta_inv, R
                REAL*8, DIMENSION(NUM_NODE, NUM_NODE) :: cov_b

                REAL*8, DIMENSION(NUM_MEAS) :: d, delta
                REAL*8, DIMENSION(NUM_NODE) :: lambda, b, u, temp

                INTEGER, DIMENSION(NUM_NODE) :: IPIV
                INTEGER, DIMENSION(NUM_MEAS) :: IPIV_R
                INTEGER :: INFO, i

                OPEN(1, file="output/LHS_lapack.dat")
                READ(1, *) K
                CLOSE(1)

                OPEN(1, file="output/S.dat")
                READ(1, *) S
                CLOSE(1)

                S_trans = TRANSPOSE(S)

                OPEN(1, file="output/cov_delta_inv.dat")
                READ(1, *) cov_delta_inv
                CLOSE(1)

                OPEN(1, file="output/cov_b.dat")
                READ(1, *) cov_b
                CLOSE(1)

                CALL DGBTRF(NUM_NODE, NUM_NODE, BANDWIDTH, BANDWIDTH, K, 3 * BANDWIDTH + 1, IPIV, INFO)

                IF (INFO .NE. 0) THEN
                        PRINT *, "Error with matrix factorization"
                END IF

                DO i=1,NUM_MEAS

                        delta = 0
                        delta(i) = 1

                        temp =  2 * MATMUL(MATMUL(S_trans, cov_delta_inv), delta)

                        CALL DGBTRS('T', NUM_NODE, BANDWIDTH, BANDWIDTH, 1, K, 3 * BANDWIDTH + 1, IPIV, temp, NUM_NODE, INFO)

                        lambda = temp

                        IF (INFO .NE. 0) THEN
                                PRINT *, "Error with lambda solve", i
                        END IF

                        b = MATMUL(cov_b, lambda) / 2.0

                        CALL DGBTRS('N', NUM_NODE, BANDWIDTH, BANDWIDTH, 1, K, 3 * BANDWIDTH + 1, IPIV, b, NUM_NODE, INFO)

                        IF (INFO .NE. 0) THEN
                                PRINT*, "Error with u solve", i
                        END IF

                        u_mat(:, i) = b
                        
                        R(:, i) = MATMUL(S, u_mat(:, i))
                        R(i, i) = R(i, i) + 1

                END DO

                CALL DGETRF(NUM_MEAS, NUM_MEAS, R, NUM_MEAS, IPIV_R, INFO)

                IF (INFO .NE. 0) THEN
                        PRINT *, "Error with R factor"
                END IF

                CALL DGETRS('N', NUM_MEAS, 1, R, NUM_MEAS, IPIV_R, d, NUM_MEAS, INFO)

                IF (INFO .NE. 0) THEN
                        PRINT *, "Error with R solve"
                END IF

                delta = d

                u = MATMUL(u_mat, delta)

                OPEN(1, file="output/u_fit.dat")
                WRITE(1, '( *(g0, ",") )') u
                CLOSE(1)

        END PROGRAM
