        FUNCTION bc(k, I0, sigma, R, theta) RESULT(f)
                INTEGER, INTENT(IN) :: k
                REAL, INTENT(IN) :: I0, sigma, R, theta
                REAL :: f

                f = -I0 / sigma * R * COS(k * theta)
        END FUNCTION

        FUNCTION get_A(r, delta_r) RESULT(A)
                REAL, INTENT(IN) :: r, delta_r
                REAL :: A

                A = (1/r)*(1/(2*delta_r)) + (1/delta_r**2)
        END FUNCTION

        FUNCTION get_B(r, delta_r) RESULT(B)
                REAL, INTENT(IN) :: r, delta_r
                REAL :: B

                B = -(1/r) * (1/(2*delta_r)) + (1/delta_r**2)
        END FUNCTION

        FUNCTION get_C(r, delta_r, delta_theta) RESULT(C)
                REAL, INTENT(IN) :: r, delta_r, delta_theta
                REAL :: C
                
                C = -2*(1/delta_r**2 + 1/(r**2 * delta_theta**2))
        END FUNCTION

        FUNCTION get_D(r, delta_theta) RESULT(D)
                REAL, INTENT(IN) :: r, delta_theta
                REAL :: D

                D = 1/(r**2 * delta_theta**2)
        END FUNCTION
        
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

        SUBROUTINE meshgrid(x, y, x_mesh, y_mesh)
                REAL, DIMENSION(:), INTENT(IN) :: x, y
                REAL, DIMENSION(:,:), INTENT(OUT):: x_mesh, y_mesh
                INTEGER :: n_x, n_y, i, j

                n_x = SIZE(x)
                n_y = SIZE(y)

                DO i=1,n_x
                        DO j=1,n_y
                                x_mesh(i, j) = x(i)
                                y_mesh(i, j) = y(j)
                        END DO
                END DO
        END SUBROUTINE 

        PROGRAM hw1
        
        INTEGER :: N, i, j
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

        N = 5
        ALLOCATE(r_array(N+1))
        ALLOCATE(theta_array(N+1))

        CALL linspace(r_min, r_max, N+1, r_array)
        CALL linspace(theta_min, theta_max, N+1, theta_array)

        delta_r = r_array(2) - r_array(1)
        delta_theta = theta_array(2) - theta_array(1)
        
        ALLOCATE(U(N**2))
        ALLOCATE(Y(N**2))
        ALLOCATE(coeff_mat(N**2, N**2))

        coeff_mat = 0.0
        Y = 0.0
        
        DO j=1,N
                DO i=1,N
                        A = get_A(r_array(i), delta_r)
                        B = get_B(r_array(i), delta_r)
                        C = get_C(r_array(i), delta_r, delta_theta)
                        D = get_D(r_array(i), delta_theta)

                        coeff_mat((j-1)*N+i, (j-1)*N+i) = C

                        IF (i == 1) THEN
                                coeff_mat((j-1)*N+i,(j-1)*N+i+1) = A+B
                        ELSE
                                coeff_mat((j-1)*N+i, (j-1)*N+i-1) = B
                        END IF 

                        IF (i == N) THEN
                                theta = theta_array(j)
                                f = bc(k, I0, sigma, r_max, theta)
                                Y((j-1)*N+i) = -A * f
                        ELSE 
                                coeff_mat((j-1)*N+i,(j-1)*N+i+1) = A
                        END IF

                        IF (j == 1) THEN
                                coeff_mat((j-1)*N+i,j*N+i) = 2*D
                        ELSE
                                coeff_mat((j-1)*N+i,(j-2)*N+i) = D
                        END IF 

                        IF (j == N) THEN
                                Y((j-1)*N+i) = 0.0
                        ELSE
                                coeff_mat((j-1)*N+i,j*N+i) = D
                        END IF
                END DO
        END DO
        OPEN(1, file="output.dat")
        DO i=1,N**2
                WRITE(1,*) coeff_mat(i,:)
        END DO
        END PROGRAM hw1
