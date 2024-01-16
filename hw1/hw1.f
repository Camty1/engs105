        PROGRAM hw1
        INTEGER :: N, i, j, k, row, col, new_col, width
        REAL :: r_min, r_max, theta_min, theta_max, delta_r, delta_theta
        REAL :: I0, sigma, theta
        REAL :: A, B, C, D, f
        REAL, DIMENSION(:), ALLOCATABLE :: r_array, theta_array, U, Y
        REAL, DIMENSION(:, :), ALLOCATABLE :: coeff_mat
        CHARACTER (LEN=30) :: filename

        r_min = 0.1
        r_max = 1.0
        theta_min = 0.0
        theta_max = 1.5707963268
        I0 = 1.0
        sigma = 1.0
        k = 3

        N = 100
        ALLOCATE(r_array(N+1))
        ALLOCATE(theta_array(N+1))
        width = N
        CALL linspace(r_min, r_max, N+1, r_array)
        CALL linspace(theta_min, theta_max, N+1, theta_array)

        delta_r = r_array(2) - r_array(1)
        delta_theta = theta_array(2) - theta_array(1)
        
        ALLOCATE(U(N**2))
        ALLOCATE(Y(N**2))
        ALLOCATE(coeff_mat(N**2, 2*width + 1))
        coeff_mat = 0.0
        Y = 0.0
        WRITE(filename, '( "output/r", I3.3, ".dat" )' ) N
        OPEN(1, file=filename)
        WRITE(filename, '( "output/theta", I3.3, ".dat" )' ) N
        OPEN(2, file=filename)
        WRITE(filename, '( "output/x", I3.3, ".dat" )' ) N
        OPEN(3, file=filename)
        WRITE(filename, '( "output/y", I3.3, ".dat" )' ) N
        OPEN(4, file=filename)

        DO j=1,N
        DO i=1,N

        WRITE(1,*) r_array(i)
        WRITE(2,*) theta_array(j)
        WRITE(3,*) r_array(i) * COS(theta_array(j))
        WRITE(4,*) r_array(i) * SIN(theta_array(j))

        CALL get_A(r_array(i), delta_r, delta_theta, A)
        CALL get_B(r_array(i), delta_r, delta_theta, B)
        CALL get_C(r_array(i), delta_r, delta_theta, C)
        CALL get_D(r_array(i), D)

        row = (j-1)*N+i
        col = (j-1)*N+i
        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = C

        IF (i == 1) THEN
        col = (j-1)*N+i+1
        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = A+B

        ELSE IF (i == N) THEN
        col = (j-1)*N+i-1
        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = B

        CALL bc(k,I0,sigma,r_max,theta_array(j), f)
        Y(row) = -A * f

        ELSE
        col = (j-1)*N+i-1
        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = B

        col = (j-1)*N+i+1
        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = A

        END IF 

        IF (j == 1) THEN
        col = j*N+i
        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = 2*D

        ELSE IF (j == N) THEN
        col = (j-2)*N+i
        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = D

        ELSE 
        col = (j-2)*N+i
        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = D

        col = j*N+i
        CALL mode2_index_map(row, col, width, new_col)
        coeff_mat(row, new_col) = D

        END IF 

        END DO
        END DO

        CLOSE(1)
        CLOSE(2)
        CLOSE(3)
        CLOSE(4)

        WRITE(filename, '( "output/coeff", I3.3, ".dat" )' ) N
        OPEN(1, file=filename)
        WRITE(filename, '( "output/b", I3.3, ".dat" )' ) N
        OPEN(2, file=filename)
        WRITE(filename, '( "output/u", I3.3, ".dat" )' ) N
        OPEN(3, file=filename)

        DO i=1,N**2
        WRITE(1,*) coeff_mat(i,:)
        WRITE(2,*) Y(i)

        END DO

        CALL SOLVE(3, coeff_mat, Y, N**2, width, N**2, 2*width+1)

        DO i=1,N**2
        WRITE(3,*) Y(i)

        END DO

        CLOSE(1)
        CLOSE(2)
        CLOSE(3)

        END PROGRAM hw1

        SUBROUTINE bc(k, I0, sigma, R, theta, f)
                INTEGER, INTENT(IN) :: k
                REAL, INTENT(IN) :: I0, sigma, R, theta
                REAL, INTENT(OUT) :: f

                f = -I0 / sigma * R * COS(k * theta)
        END SUBROUTINE

        SUBROUTINE get_A(r, delta_r, delta_theta, A)
                REAL, INTENT(IN) :: r, delta_r, delta_theta
                REAL, INTENT(OUT) :: A

                A = delta_theta**2/(2*r*delta_r)
                A = A + delta_theta**2/delta_r**2
        END SUBROUTINE

        SUBROUTINE get_B(r, delta_r, delta_theta, B)
                REAL, INTENT(IN) :: r, delta_r, delta_theta
                REAL, INTENT(OUT) :: B

                B = -delta_theta**2/(2*r*delta_r)+delta_theta**2/delta_r**2
        END SUBROUTINE

        SUBROUTINE get_C(r, delta_r, delta_theta, C)
                REAL, INTENT(IN) :: r, delta_r, delta_theta
                REAL, INTENT(OUT) :: C
                
                C = -2*(delta_theta**2/delta_r**2+1/r**2)
        END SUBROUTINE

        SUBROUTINE get_D(r, D)
                REAL, INTENT(IN) :: r
                REAL, INTENT(OUT) :: D

                D = 1/r**2
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
