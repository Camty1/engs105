        PROGRAM hw2
        IMPLICIT REAL*8 (A-H, O-Z)

        INTEGER, PARAMETER :: N = 100
        REAL*8, PARAMETER :: alpha = 1.0, threshold = 1E-5
        
        REAL*8, DIMENSION(N) :: r_array, theta_array
        REAL*8, DIMENSION(N**2) :: U, U_true, Y, Y_cpy
        REAL*8, DIMENSION(N**2, 1) :: diag
        REAL*8, DIMENSION(N**2, 2*N+1) :: coeff_mat, lower, upper
        REAL*8, DIMENSION(100**2) :: U_100

        INTEGER :: i
        REAL*8 :: err

        coeff_mat = 0.0
        Y = 0.0
        U = 0.0

        OPEN(1, file="output/u100_1.dat")
        READ(1,*) U_100
        CLOSE(1)

        OPEN(1, file="output/u020_1.dat")

        DO i=1,N**2
                CALL INTERP_TRUE(i, U_100, N, 100, U_true(i))
                WRITE(1,*) U_true(i)
        END DO
        U = U_true/2
        CLOSE(1)
        
        CALL GEN_SYS_EQNS(N, alpha, coeff_mat, Y)
        Y_cpy = Y

        CALL DIAG_DOMINANCE_CHECK(coeff_mat, Y, N)
        diag = 0.0
        lower = 0.0
        upper = 1

        lower(:, 1:N) = coeff_mat(:, 1:N)
        upper(:, N+2:2*N+1) = coeff_mat(:, N+2:2*N+1)

        DO i=1,N**2
                diag(i,1) = coeff_mat(i, N+1)
        END DO
        
        CALL CALC_ERROR(U, U_true, N**2, err)

        DO WHILE (err > threshold)
                PRINT *, err
                CALL JACOBI(diag, lower, upper, Y, N, U)
                CALL CALC_ERROR(U, U_true, N**2, err)
        END DO

        END PROGRAM hw2

        SUBROUTINE GEN_SYS_EQNS(N, alpha, COEFF, Y)
                INTEGER, INTENT(IN) :: N
                REAL*8, INTENT(IN) :: alpha
                REAL*8, DIMENSION(N**2, 2*N+1), INTENT(INOUT) :: COEFF
                REAL*8, DIMENSION(N**2), INTENT(INOUT) :: Y

                INTEGER :: i, j, k, row, col, new_col
                REAL*8 :: r_min, r_max, theta_min, theta_max, dr, dtheta
                REAL*8 :: I0, sigma
                REAL*8 :: A, B, C, D, E, f, g
                REAL*8, DIMENSION(N) :: r_array, theta_array 
                CHARACTER (LEN=30) :: fname

                r_min = 0.1
                r_max = 1.0
                theta_min = 0
                theta_max = 1.5707963268
                I0 = 1.0
                sigma = 1.0
                k = 3

                CALL linspace(r_min, r_max, N, r_array)
                CALL linspace(theta_min, theta_max, N, theta_array)

                dr = r_array(2) - r_array(1)
                dtheta = theta_array(2) - theta_array(1)
                
                COEFF = 0.0
                Y = 0.0
                WRITE(fname,'("output/r",I3.3,".dat")')N
                OPEN(1, file=fname)
                WRITE(fname,'("output/theta",I3.3,".dat")')N
                OPEN(2, file=fname)
                WRITE(fname,'("output/x",I3.3,".dat")')N
                OPEN(3, file=fname)
                WRITE(fname,'("output/y",I3.3,".dat")')N
                OPEN(4, file=fname)

                DO j=1,N
                DO i=1,N

                WRITE(1,*) r_array(i)
                WRITE(2,*) theta_array(j)
                WRITE(3,*) r_array(i) * COS(theta_array(j))
                WRITE(4,*) r_array(i) * SIN(theta_array(j))

                CALL get_A(r_array(i), dr, dtheta, A)
                CALL get_B(r_array(i), dr, dtheta, B)
                CALL get_C(r_array(i), dr, dtheta, C)
                CALL get_D(r_array(i), D)

                row = (j-1)*N+i
                col = (j-1)*N+i
                CALL mode2_index_map(row, col, N, new_col)
                COEFF(row, new_col) = C

                IF (j == 1) THEN
                col = j*N+i
                CALL mode2_index_map(row, col, N, new_col)
                COEFF(row, new_col) = 2*D

                ELSE IF (j == N) THEN
                CALL get_E(r_array(i), dtheta, alpha, E)
                CALL get_g(r_array(i), dtheta, k, I0, sigma, g)

                col = (j-1)*N+i
                CALL mode2_index_map(row, col, N, new_col)
                COEFF(row, new_col) = C+E*D

                col = (j-2)*N+i
                CALL mode2_index_map(row, col, N, new_col)
                COEFF(row, new_col) = 2*D

                Y(row) = -g*D

                ELSE 
                col = (j-2)*N+i
                CALL mode2_index_map(row, col, N, new_col)
                COEFF(row, new_col) = D

                col = j*N+i
                CALL mode2_index_map(row, col, N, new_col)
                COEFF(row, new_col) = D

                END IF 
                
                IF (i == 1) THEN
                col = (j-1)*N+i+1
                CALL mode2_index_map(row, col, N, new_col)
                COEFF(row, new_col) = A+B

                ELSE IF (i == N) THEN
                col = (j-1)*N+i
                CALL mode2_index_map(row, col, N, new_col)
                COEFF(row, :) = 0
                COEFF(row, new_col) = 1

                CALL bc(k,I0,sigma,r_max,theta_array(j), f)
                Y(row) = f

                ELSE
                col = (j-1)*N+i-1
                CALL mode2_index_map(row, col, N, new_col)
                COEFF(row, new_col) = B

                col = (j-1)*N+i+1
                CALL mode2_index_map(row, col, N, new_col)
                COEFF(row, new_col) = A

                END IF 

                END DO
                END DO

                CLOSE(1)
                CLOSE(2)
                CLOSE(3)
                CLOSE(4)

                WRITE(fname, '( "output/coeff", I3.3, ".dat" )' ) N
                OPEN(1, file=fname)
                WRITE(fname, '( "output/b", I3.3, ".dat" )' ) N
                OPEN(2, file=fname)

                DO i=1,N**2
                WRITE(1,*) COEFF(i,:)
                WRITE(2,*) Y(i)

                END DO

                CLOSE(1)
                CLOSE(2)
                
        END SUBROUTINE

        SUBROUTINE bc(k, I0, sigma, R, theta, f)
                INTEGER, INTENT(IN) :: k
                REAL*8, INTENT(IN) :: I0, sigma, R, theta
                REAL*8, INTENT(OUT) :: f

                f = -I0 / sigma * R * COS(k * theta)
        END SUBROUTINE

        SUBROUTINE get_A(r, delta_r, delta_theta, A)
                REAL*8, INTENT(IN) :: r, delta_r, delta_theta
                REAL*8, INTENT(OUT) :: A

                A = delta_theta**2*(1/(2*r*delta_r) + 1/delta_r**2)
        END SUBROUTINE

        SUBROUTINE get_B(r, delta_r, delta_theta, B)
                REAL*8, INTENT(IN) :: r, delta_r, delta_theta
                REAL*8, INTENT(OUT) :: B

                B = delta_theta**2*(-1/(2*r*delta_r) + 1/delta_r**2)
        END SUBROUTINE

        SUBROUTINE get_C(r, delta_r, delta_theta, C)
                REAL*8, INTENT(IN) :: r, delta_r, delta_theta
                REAL*8, INTENT(OUT) :: C
                
                C = -2*(delta_theta**2/delta_r**2+1/r**2)
        END SUBROUTINE

        SUBROUTINE get_D(r, D)
                REAL*8, INTENT(IN) :: r
                REAL*8, INTENT(OUT) :: D

                D = 1/r**2
        END SUBROUTINE

        SUBROUTINE get_E(r, delta_theta, alpha, E)
                REAL*8, INTENT(IN) :: r, delta_theta, alpha
                REAL*8, INTENT(OUT) :: E

                E = -2*r*delta_theta*alpha
        END SUBROUTINE

        SUBROUTINE get_g(r, delta_theta, k, I0, sigma, g)
                REAL*8, INTENT(IN) :: r, delta_theta, I0, sigma
                INTEGER, INTENT(IN) :: k
                REAL*8, INTENT(OUT) :: g

                g = -2*r*delta_theta*k*I0/sigma
        END SUBROUTINE
        
        SUBROUTINE linspace(x_min, x_max, num_points, arr)
                REAL*8, INTENT(IN) :: x_min, x_max
                INTEGER, INTENT(IN) :: num_points
                REAL*8, DIMENSION(num_points), INTENT(OUT) :: arr
                INTEGER :: i

                DO i=1,num_points
                        arr(i) = (i-1) * (x_max-x_min)/(num_points-1)
                        arr(i) = arr(i) + x_min
                END DO
        END SUBROUTINE

        SUBROUTINE polar_2_cart(r, theta, x, y)
                REAL*8, INTENT(IN) :: r, theta
                REAL*8, INTENT(OUT) :: x, y

                x = r * COS(theta)
                y = r * SIN(theta)
        END SUBROUTINE 
        
        SUBROUTINE mode2_index_map(row, col, width, new_col)
                INTEGER, INTENT(IN) :: row, col, width
                INTEGER, INTENT(OUT) :: new_col

                new_col = (width + 1) + (col-row)
        END SUBROUTINE

        SUBROUTINE JACOBI(diag, lower, upper, Y, N, U)
                REAL*8, DIMENSION(N**2,1), INTENT(IN) :: diag
                INTEGER, INTENT(IN) :: N
                REAL*8, DIMENSION(N**2,2*N+1),INTENT(IN) :: lower, upper
                REAL*8, DIMENSION(N**2), INTENT(IN) :: Y
                REAL*8, DIMENSION(N**2), INTENT(INOUT) :: U

                CALL BANDED_MULT_IN_PLACE(-(lower+upper),U,N**2,2*N+1,N)

                U = U + Y

                CALL DSOLVE(3, diag, U, N, 0, N, 1)

        END SUBROUTINE

        SUBROUTINE GAUSS_SIDEL(diag, lower, upper, Y, N, U)
                REAL*8, DIMENSION(N**2,1), INTENT(IN) :: diag
                INTEGER, INTENT(IN) :: N
                REAL*8, DIMENSION(N**2,2*N+1),INTENT(IN) :: lower, upper
                REAL*8, DIMENSION(N**2), INTENT(IN) :: Y
                REAL*8, DIMENSION(N**2), INTENT(INOUT) :: U
                REAL*8, DIMENSION(N**2,2*N+1) :: lhs

                CALL BANDED_MULT_IN_PLACE(-upper,U,N**2,2*N+1,N)

                U = U + Y
                lhs = lower
                DO i=1,N**2
                        lhs(i, N+1) = diag(i,1)
                END DO

                CALL DSOLVE(3, lhs, U, N, 0, N, 1)

        END SUBROUTINE

        SUBROUTINE CALC_ERROR(U, U_TRUE, N, err)
                INTEGER, INTENT(IN) :: N
                REAL*8, DIMENSION(N), INTENT(IN) :: U, U_TRUE
                REAL*8, INTENT(OUT) :: err

                err = MAXVAL(ABS(U-U_TRUE))
        END SUBROUTINE

        SUBROUTINE INTERP_TRUE(i, U_true, N, N_true, val)
                INTEGER, INTENT(IN) :: i, N, N_true
                REAL*8, DIMENSION(N), INTENT(IN) :: U_true
                REAL*8, INTENT(OUT) :: val

                INTEGER :: group, item, i_cpy, a, group_flag, item_flag
                INTEGER :: idx
                REAL*8 :: ratio, r_i, r_g, alpha, beta, gamma, delta

                ratio = DBLE(N_true-1)/DBLE(N-1)

                group = (i-1)/N + 1
                item = MOD(i-1, N)+1

                r_i = ratio * DBLE(item-1)
                r_i = r_i - DBLE(INT(r_i))
                r_g = ratio * DBLE(group-1)
                r_g = r_g - DBLE(INT(r_g))

                group = FLOOR(DBLE(group-1)*ratio)+1
                item = FLOOR(DBLE(item-1)*ratio)+1

                idx = N_true * (group - 1) + item
                alpha = (1 - r_i) * (1 - r_g) * U_true(idx)

                idx = N_true * (group - 1) + item + 1
                beta = r_i * (1 - r_g) * U_true(idx)

                idx = N_true * group + item
                gamma = (1 - r_i) * r_g * U_true(idx)

                idx = N_true * group + item + 1
                delta = r_i * r_g * U_true(idx)

                val = alpha + beta + gamma + delta

        END SUBROUTINE

        SUBROUTINE DIAG_DOMINANCE_CHECK(COEFF, Y, N)
                INTEGER, INTENT(IN) :: N
                REAL*8, DIMENSION(N**2, 2*N+1), INTENT(INOUT) :: COEFF
                REAL*8, DIMENSION(N**2), INTENT(INOUT) :: Y

                INTEGER :: i, j
                REAL*8 :: diag, non_diag

                DO i=1,N**2
                non_diag = 0
                DO j=1,2*N+1
                        non_diag = non_diag + ABS(COEFF(i,j))
                END DO
                diag = COEFF(i, N+1)
                non_diag = non_diag - diag

                IF (ABS(non_diag) > ABS(diag)) THEN
                        COEFF(i,:) = -COEFF(i,:)
                        Y(i) = -Y(i)
                        non_diag = 0
                        diag = 0
                        DO j=1,2*N+1
                                non_diag = non_diag + ABS(COEFF(i,j))
                        END DO
                        diag = COEFF(i, N+1)
                        non_diag = non_diag - diag
                        IF (ABS(non_diag) > ABS(diag)) THEN
                                PRINT *, i, diag, non_diag
                        END IF
                END IF
                END DO

        END SUBROUTINE

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
