        SUBROUTINE bc(k, I0, sigma, R, theta, f)
                INTEGER, INTENT(IN) :: k
                REAL, INTENT(IN) :: I0, sigma, R, theta
                REAL, INTENT(OUT) :: f

                f = -I0 / sigma * R * COS(k * theta)
        END SUBROUTINE

        SUBROUTINE get_A(r, delta_r, A)
                REAL, INTENT(IN) :: r, delta_r
                REAL, INTENT(OUT) :: A

                A = (1/r)*(1/(2*delta_r)) + (1/delta_r**2)
        END SUBROUTINE

        SUBROUTINE get_B(r, delta_r, B)
                REAL, INTENT(IN) :: r, delta_r
                REAL, INTENT(OUT) :: B

                B = -(1/r) * (1/(2*delta_r)) + (1/delta_r**2)
        END SUBROUTINE

        SUBROUTINE get_C(r, delta_r, delta_theta, C)
                REAL, INTENT(IN) :: r, delta_r, delta_theta
                REAL, INTENT(OUT) :: C
                
                C = -2*(1/delta_r**2 + 1/(r**2 * delta_theta**2))
        END SUBROUTINE

        SUBROUTINE get_D(r, delta_theta, D)
                REAL, INTENT(IN) :: r, delta_theta
                REAL, INTENT(OUT) :: D

                D = 1/(r**2 * delta_theta**2)
        END SUBROUTINE
        
        SUBROUTINE linspace(x_min, x_max, num_points, arr)
                REAL, INTENT(IN) :: x_min, x_max
                INTEGER, INTENT(IN) :: num_points
                REAL, DIMENSION(num_points), INTENT(OUT) :: arr
                INTEGER :: i

                DO i=1,num_points
                        arr(i) = (i-1) * (x_max-x_min)/(num_points-1)
                        arr(i) = arr(i) + x_min
                END DO
        END SUBROUTINE

        SUBROUTINE polar_2_cart(r, theta, x, y)
                REAL, INTENT(IN) :: r, theta
                REAL, INTENT(OUT) :: x, y

                x = r * COS(theta)
                y = r * SIN(theta)
        END SUBROUTINE 
        
        SUBROUTINE mode2_index_map(row, col, width, new_col)
                INTEGER, INTENT(IN) :: row, col, width
                INTEGER, INTENT(OUT) :: new_col

                new_col = (width + 1) + (col-row)
        END SUBROUTINE

        PROGRAM hw1
        INTEGER :: N, i, j, k, row, col, new_col, width
        REAL :: r_min, r_max, theta_min, theta_max, delta_r, delta_theta
        REAL :: I0, sigma, theta
        REAL :: A, B, C, D, f
        REAL, DIMENSION(:), ALLOCATABLE :: r_array, theta_array, U, Y
        REAL, DIMENSION(:, :), ALLOCATABLE :: coeff_mat

        r_min = 0.1
        r_max = 1.0
        theta_min = 0.0
        theta_max = 1.5707963268
        I0 = 1.0
        sigma = 1.0
        k = 3

        width = 5

        N = 50
        ALLOCATE(r_array(N+1))
        ALLOCATE(theta_array(N+1))

        CALL linspace(r_min, r_max, N+1, r_array)
        CALL linspace(theta_min, theta_max, N+1, theta_array)

        delta_r = r_array(2) - r_array(1)
        delta_theta = theta_array(2) - theta_array(1)
        
        ALLOCATE(U(N**2))
        ALLOCATE(Y(N**2))
        ALLOCATE(coeff_mat(N**2, 2*width + 1))
        coeff_mat = 0.0
        Y = 0.0
        
        DO j=1,N
        DO i=1,N
        CALL get_A(r_array(i), delta_r, A)
        CALL get_B(r_array(i), delta_r, B)
        CALL get_C(r_array(i), delta_r, delta_theta, C)
        CALL get_D(r_array(i), delta_theta, D)

        row = (j-1)*N+i
        col = (j-1)*N+i

        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = C

        IF (i == 1) THEN
        row = (j-1)*N+i
        col = (j-1)*N+i+1
        
        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = A+B

        ELSE
        row = (j-1)*N+i
        col = (j-1)*N+i-1

        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = B

        END IF 

        IF (i == N) THEN
        CALL bc(k,I0,sigma,r_max,theta_array(j), f)
        row = (j-1)*N+i
        Y(row) = -A * f

        ELSE 
        row = (j-1)*N+i
        col = (j-1)*N+i+1

        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = A

        END IF

        IF (j == 1) THEN
        row = (j-1)*N+i
        col = j*N+i

        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = 2*D

        ELSE
        row = (j-1)*N+i
        col = (j-2)*N+i

        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = D

        END IF 

        IF (j < N) THEN
        row = (j-1)*N+i
        col = j*N+i

        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = D

        END IF

        END DO
        END DO

        OPEN(1, file="coeff.dat")
        OPEN(2, file="y.dat")
        OPEN(3, file="u.dat")

        DO i=1,N**2
        PRINT *, coeff_mat(i,:), y(i)
        WRITE(1,*) coeff_mat(i,:)
        WRITE(2,*) Y(i)

        END DO

        CALL SOLVE(3, coeff_mat, Y, N**2, width, N**2, 2*width+1)

        DO i=1,N**2
        PRINT *, y(i)
        WRITE(3,*) Y(i)

        END DO

        CLOSE(1)
        CLOSE(2)
        CLOSE(3)

        END PROGRAM hw1
