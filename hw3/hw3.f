        PROGRAM hw2
        IMPLICIT REAL*8 (A-H, O-Z)
        
        INTEGER, PARAMETER :: N = 80, k = 3
        REAL*8, PARAMETER :: r_min = 0.1, r_max = 1, phi_min = 0
        REAL*8, PARAMETER :: phi_max = 1.570796326794897
        REAL*8, PARAMETER :: sigma = 1.0, I0 = 1.0

        END PROGRAM

        SUBROUTINE get_A(r, dr, A)
                REAL*8, INTENT(IN) :: r, dr
                REAL*8, INTENT(OUT) :: A

                A = (r + dr/2)/(r * dr**2)

        END SUBROUTINE

        SUBROUTINE get_B(r, dr, B)
                REAL*8, INTENT(IN) :: r, dr
                REAL*8, INTENT(OUT) :: B

                B = (r - dr/2)/(r * dr**2)

        END SUBROUTINE

        SUBROUTINE get_C(r, dr, dphi, C)
                REAL*8, INTENT(IN) :: r, dr, dphi
                REAL*8, INTENT(OUT) :: C

                C = -2 * (1 / dr**2 + 1 / (r * dphi)**2)

        END SUBROUTINE

        SUBROUTINE get_D(r, dphi, D)
                REAL*8, INTENT(IN) :: r, dphi
                REAL*8, INTENT(OUT) :: D

                D = 1 / (r * dphi) ** 2

        END SUBROUTINE

        SUBROUTINE BC(k, I0, sigma, r_max, phi, f)
                INTEGER, INTENT(IN) :: k
                REAL*8, INTENT(IN) :: I0, sigma, r_max, phi
                REAL*8, INTENT(OUT) :: f
                
                f = - I0 / sigma * r_max * COS(k * phi)
        END SUBROUTINE

        ! TODO
        SUBROUTINE get_E(r, dphi, alpha, E)

        END SUBROUTINE

        ! TODO
        SUBROUTINE get_g(r, dphi, k, I0, sigma, g)

        END SUBROUTINE
         
        SUBROUTINE GEN_SS_EQNS_TYPE_I(N, r_arr, phi_arr, k, I0, sigma, A_mat, b_vec)
                INTEGER, INTENT(IN) :: N, k
                REAL*8, INTENT(IN) :: I0, sigma

                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
                REAL*8, DIMENSION(N**2, 2*N+1), INTENT(OUT) :: A_mat
                REAL*8, DIMENSION(N**2), INTENT(OUT) :: B_vec
                
                INTEGER :: i, j, row, col, new_col
                
                REAL*8 :: dr, dphi, A, B, C, D, f, r_max

                r_max = MAXVAL(r_arr)

                dr = r_arr(2) - r_arr(1)
                dphi = phi_arr(2) - phi_arr(1)

                A_mat = 0.0
                B_vec = 0.0

                DO j=1,N
                        DO i=1,N
                                CALL get_A(r_arr(i), dr, A)
                                CALL get_B(r_arr(i), dr, B)
                                CALL get_C(r_arr(i), dr, dphi, C)
                                CALL get_D(r_arr(i), dphi, D)
                                row = N*(j-1)+i
                                IF (i == 1 .AND. j == 1) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = C

                                        col = N*(j-1)+i+1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = A + B

                                        col = N*j+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = 2 * D
                                              
                                ELSE IF (i == 1 .AND. 1 < j .AND. j < N) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = C

                                        col = N*(j-1)+i+1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = A + B

                                        col = N*(j-2)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = D

                                        col = N*j+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = D

                                ELSE IF (1 < i .AND. i < N .AND. j == 1) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = C

                                        col = N*(j-1)+i+1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = A

                                        col = N*(j-1)+i-1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = B

                                        col = N*j+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = 2 * D
                                
                                ELSE IF (1 < i .AND. i < N .AND. 1 < j .AND. j < N) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = C

                                        col = N*(j-1)+i+1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = A

                                        col = N*(j-1)+i-1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = B

                                        col = N*(j-2)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = D

                                        col = N*j+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = D

                                ELSE IF (i == N) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        CALL BC(k, I0, sigma, r_max, phi_arr(j), f)
                                        A_mat(row, new_col) = 1
                                        B_vec(row) = f

                                ELSE IF (j == N) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = 1
                                        B_vec(row) = 0

                                ELSE 
                                        PRINT *, "You should not be here"

                                END IF
                        END DO
                END DO
        END SUBROUTINE

        ! TODO
        SUBROUTINE GEN_SS_EQNS_TYPE_III(N, alpha, r_arr, phi_arr, k, I0, sigma, A_mat, b_vec)
                INTEGER, INTENT(IN) :: N, k
                REAL*8, INTENT(IN) :: alpha, I0, sigma

                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
                REAL*8, DIMENSION(N**2, 2*N+1), INTENT(OUT) :: A_mat
                REAL*8, DIMENSION(N**2), INTENT(OUT) :: B_vec
                
                INTEGER :: i, j, row, col, new_col
                
                REAL*8 :: dr, dphi, A, B, C, D, E, f, g, r_max

                r_max = MAXVAL(r_arr)


                A_mat = 0.0
                B_vec = 0.0
                PRINT *, "BONGO"
                DO j=1,N
                        DO i=1,N
                                CALL get_A(r_arr(i), dr, A)
                                CALL get_B(r_arr(i), dr, B)
                                CALL get_C(r_arr(i), dr, dphi, C)
                                CALL get_D(r_arr(i), dphi, D)
                                row = N*(j-1)+i
                                IF (i == 1 .AND. j == 1) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = C

                                        col = N*(j-1)+i+1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = A + B

                                        col = N*j+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = 2 * D
                                              
                                ELSE IF (i == 1 .AND. 1 < j .AND. j < N) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = C

                                        col = N*(j-1)+i+1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = A + B

                                        col = N*(j-2)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = D

                                        col = N*j+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = D

                                ELSE IF (i == 1 .AND. j == N) THEN

                                ELSE IF (1 < i .AND. i < N .AND. j == 1) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = C

                                        col = N*(j-1)+i+1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = A

                                        col = N*(j-1)+i-1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = B

                                        col = N*j+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = 2 * D
                                
                                ELSE IF (1 < i .AND. i < N .AND. 1 < j .AND. j < N) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = C

                                        col = N*(j-1)+i+1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = A

                                        col = N*(j-1)+i-1
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = B

                                        col = N*(j-2)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = D

                                        col = N*j+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = D

                                ELSE IF (1 < i .AND. i < N .AND. j == N) THEN

                                ELSE IF (i == N) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        CALL BC(k, I0, sigma, r_max, phi_arr(j), f)
                                        A_mat(row, new_col) = 1
                                        B_vec(row) = f

                                ELSE IF (j == N) THEN
                                        col = N*(j-1)+i
                                        CALL mode2_index_map(row, col, N, new_col)

                                        A_mat(row, new_col) = 1
                                        B_vec(row) = 0

                                ELSE 
                                        PRINT *, "You should not be here"

                                END IF
                        END DO
                END DO
        END SUBROUTINE

        SUBROUTINE GET_LHS_TYPE_I(N, theta, dt, r_arr, phi_arr, k, I0, sigma, A_mat)
                INTEGER, INTENT(IN) :: N, k
                REAL*8, INTENT(IN) :: I0, sigma, theta, dt

                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
                REAL*8, DIMENSION(N**2, 2*N+1), INTENT(OUT) :: A_mat
                
                INTEGER :: i, j, row, col, new_col
                
                REAL*8 :: dr, dphi, A, B, C, D, f, r_max

                r_max = MAXVAL(r_arr)

                dr = r_arr(2) - r_arr(1)
                dphi = phi_arr(2) - phi_arr(1)

                A_mat = 0.0
                B_vec = 0.0

                DO j=1,N
                DO i=1,N

                CALL get_A(r_arr(i), dr, A)
                CALL get_B(r_arr(i), dr, B)
                CALL get_C(r_arr(i), dr, dphi, C)
                CALL get_D(r_arr(i), dphi, D)

                row = N*(j-1)+i

                IF (i == 1 .AND. j == 1) THEN
                        col = N*(j-1)+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -(A + B) * theta * dt

                        col = N*j+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -2 * D * theta * dt
                              
                ELSE IF (i == 1 .AND. 1 < j .AND. j < N) THEN
                        col = N*(j-1)+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -(A + B) * theta * dt

                        col = N*(j-2)+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                        col = N*j+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                ELSE IF (1 < i .AND. i < N .AND. j == 1) THEN
                        col = N*(j-1)+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -A * theta * dt

                        col = N*(j-1)+i-1
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -B * theta * dt

                        col = N*j+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -2 * D * theta * dt
                
                ELSE IF (1 < i .AND. i < N .AND. 1 < j .AND. j < N) THEN
                        col = N*(j-1)+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = 1 - C * theta * dt

                        col = N*(j-1)+i+1
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -A * theta * dt

                        col = N*(j-1)+i-1
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -B * theta * dt

                        col = N*(j-2)+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                        col = N*j+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = -D * theta * dt

                ELSE IF (i == N) THEN
                        col = N*(j-1)+i
                        CALL mode2_index_map(row, col, N, new_col)

                        CALL BC(k, I0, sigma, r_max, phi_arr(j), f)
                        A_mat(row, new_col) = 1

                ELSE IF (j == N) THEN
                        col = N*(j-1)+i
                        CALL mode2_index_map(row, col, N, new_col)

                        A_mat(row, new_col) = 1

                ELSE 
                        PRINT *, "You should not be here"

                END IF

                END DO
                END DO
        END SUBROUTINE

        ! TODO
!        SUBROUTINE GET_RHS_TYPE_I(N, theta, dt, r_arr, phi_arr, k, I0, sigma, U, U_new)
!                INTEGER, INTENT(IN) :: N, k
!                REAL*8, INTENT(IN) :: I0, sigma, theta, dt
!
!                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
!                REAL*8, DIMENSION(N**2), INTENT(IN) :: U
!                REAL*8, DIMENSION(N**2), INTENT(OUT) :: U_new
!                
!                INTEGER :: i, j, idx
!                
!                REAL*8 :: dr, dphi, A, B, C, D, f
!
!                dr = r_arr(2) - r_arr(1)
!                dphi = phi_arr(2) - phi_arr(1)
!
!                A_mat = 0.0
!                B_vec = 0.0
!
!                DO j=1,N
!                DO i=1,N
!
!                CALL get_A(r_arr(i), dr, A)
!                CALL get_B(r_arr(i), dr, B)
!                CALL get_C(r_arr(i), dr, dphi, C)
!                CALL get_D(r_arr(i), dphi, D)
!
!                row = N*(j-1)+i
!
!                IF (i == 1 .AND. j == 1) THEN
!                        col = N*(j-1)+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = 1 - C * theta * dt
!
!                        col = N*(j-1)+i+1
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -(A + B) * theta * dt
!
!                        col = N*j+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -2 * D * theta * dt
!                              
!                ELSE IF (i == 1 .AND. 1 < j .AND. j < N) THEN
!                        col = N*(j-1)+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = 1 - C * theta * dt
!
!                        col = N*(j-1)+i+1
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -(A + B) * theta * dt
!
!                        col = N*(j-2)+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -D * theta * dt
!
!                        col = N*j+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -D * theta * dt
!
!                ELSE IF (1 < i .AND. i < N .AND. j == 1) THEN
!                        col = N*(j-1)+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = 1 - C * theta * dt
!
!                        col = N*(j-1)+i+1
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -A * theta * dt
!
!                        col = N*(j-1)+i-1
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -B * theta * dt
!
!                        col = N*j+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -2 * D * theta * dt
!                
!                ELSE IF (1 < i .AND. i < N .AND. 1 < j .AND. j < N) THEN
!                        col = N*(j-1)+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = 1 - C * theta * dt
!
!                        col = N*(j-1)+i+1
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -A * theta * dt
!
!                        col = N*(j-1)+i-1
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -B * theta * dt
!
!                        col = N*(j-2)+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -D * theta * dt
!
!                        col = N*j+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = -D * theta * dt
!
!                ELSE IF (i == N) THEN
!                        col = N*(j-1)+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        CALL BC(k, I0, sigma, r_max, phi_arr(j), f)
!                        A_mat(row, new_col) = 1
!
!                ELSE IF (j == N) THEN
!                        col = N*(j-1)+i
!                        CALL mode2_index_map(row, col, N, new_col)
!
!                        A_mat(row, new_col) = 1
!
!                ELSE 
!                        PRINT *, "You should not be here"
!
!                END IF
!
!                END DO
!                END DO
!
!        END SUBROUTINE

        ! TODO
!        SUBROUTINE GET_LHS_TYPE_III(N, alpha, theta, dt, r_arr, phi_arr, k, I0, sigma, A_mat)
!                INTEGER, INTENT(IN) :: N, k
!                REAL*8, INTENT(IN) :: I0, sigma, alpha, theta, dt
!
!                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
!                REAL*8, DIMENSION(N**2, 2*N+1), INTENT(OUT) :: A_mat
!                
!                INTEGER :: i, j, row, col, new_col
!                
!                REAL*8 :: dr, dphi, A, B, C, D, E, f, g
!
!        END SUBROUTINE

        ! TODO
!        SUBROUTINE GET_RHS_TYPE_III(N, alpha, theta, dt, r_arr, phi_arr, k, I0, sigma, U)
!                INTEGER, INTENT(IN) :: N, k
!                REAL*8, INTENT(IN) :: I0, sigma, alpha, theta, dt
!
!                REAL*8, DIMENSION(N), INTENT(IN) :: r_arr, phi_arr
!                REAL*8, DIMENSION(N**2), INTENT(INOUT) :: U
!                
!                INTEGER :: i, j, row, col, new_col
!                
!                REAL*8 :: dr, dphi, A, B, C, D, E, f, g
!
!        END SUBROUTINE

        SUBROUTINE mode2_index_map(row, col, width, new_col)
                INTEGER, INTENT(IN) :: row, col, width
                INTEGER, INTENT(OUT) :: new_col

                new_col = (width + 1) + (col-row)
        END SUBROUTINE
