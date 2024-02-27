        PROGRAM HW5
                INTEGER, PARAMETER :: NUM_ELEM=1872, NUM_NODE=969, NUM_BC=64, NUM_MAT=8, TIMESTEPS=5000
                REAL*8, PARAMETER :: theta=0.5, dt=1

                INTEGER, DIMENSION(6, NUM_ELEM) :: elem_list_input
                INTEGER, DIMENSION(4, NUM_ELEM) :: elem_list
                INTEGER, DIMENSION(NUM_ELEM) :: material_list

                REAL*8, DIMENSION(3, NUM_NODE) :: node_pos_input
                REAL*8, DIMENSION(2, NUM_NODE) :: node_pos

                REAL*8, DIMENSION(7, NUM_BC) :: bc_list_input
                INTEGER, DIMENSION(NUM_BC) :: bc_nodes
                INTEGER, DIMENSION(NUM_BC) :: bc_prev_nodes
                INTEGER, DIMENSION(NUM_BC) :: bc_next_nodes
                REAL*8, DIMENSION(NUM_BC) :: bc_h
                REAL*8, DIMENSION(NUM_BC) :: bc_ta

                REAL*8, DIMENSION(3, NUM_ELEM) :: heat_rate_list_input
                REAL*8, DIMENSION(NUM_ELEM) :: heat_rate_list

                REAL*8, DIMENSION(4, NUM_MAT - 2) :: mat_properties_input
                REAL*8, DIMENSION(3, NUM_MAT) :: mat_properties

                REAL*8, DIMENSION(:,:), ALLOCATABLE :: A_mat, B_mat
                REAL*8, DIMENSION(NUM_NODE) :: RHS, U

                REAL*8, DIMENSION(3,3) :: elem_mat_3
                REAL*8, DIMENSION(4,4) :: elem_mat_4
                REAL*8, DIMENSION(4) :: tri_quad_points
                REAL*8, DIMENSION(3) :: b_vec_3
                REAL*8, DIMENSION(4) :: b_vec_4

                REAL*8, DIMENSION(4, NUM_ELEM) :: x_elem, y_elem, sigmaE
                REAL*8, DIMENSION(3) :: dx, dy
                REAL*8 :: A, kappa, m, c_rho

                INTEGER :: i, j, k, l, bandwidth, width, row, col, new_col

                OPEN(1, file="epeltr4_tri.dat")
                READ(1, *) elem_list_input
                CLOSE(1)

                OPEN(1, file="npeltr4.dat")
                READ(1, *) node_pos_input
                CLOSE(1)

                OPEN(1, file="bpeltr4.dat")
                READ(1, *) bc_list_input
                CLOSE(1)

                OPEN(1, file="bpelt4_tri.dat")
                READ(1, *) heat_rate_list_input
                CLOSE(1)

                OPEN(1, file="material_properties.dat")
                READ(1, *) mat_properties_input
                CLOSE(1)

                elem_list = elem_list_input(2:5, :)
                material_list = elem_list_input(6, :)

                node_pos = node_pos_input(2:3, :)

                bc_nodes = INT(bc_list_input(2, :))
                bc_prev_nodes = INT(bc_list_input(4, :))
                bc_next_nodes = INT(bc_list_input(5, :))
                bc_h = bc_list_input(6,:)
                bc_ta = bc_list_input(7,:)

                heat_rate_list = heat_rate_list_input(3,:)

                mat_properties(:, 3:8) = mat_properties_input(2:4, :)

                bandwidth = 0

                DO i=1,NUM_ELEM

                x_elem(1,i) = node_pos(1,elem_list(1,i))
                x_elem(2,i) = node_pos(1,elem_list(2,i))
                x_elem(3,i) = node_pos(1,elem_list(3,i))
                x_elem(4,i) = node_pos(1,elem_list(4,i))

                y_elem(1,i) = node_pos(2,elem_list(1,i))
                y_elem(2,i) = node_pos(2,elem_list(2,i))
                y_elem(3,i) = node_pos(2,elem_list(3,i))
                y_elem(4,i) = node_pos(2,elem_list(4,i))

                DO j=1,4
                DO k=j,4
                width = ABS(elem_list(j, i) - elem_list(k, i))

                IF (width > bandwidth) THEN
                bandwidth = width

                END IF

                END DO
                END DO
                END DO

                width = 2 * bandwidth + 1

                ALLOCATE(A_mat(NUM_NODE, width))
                ALLOCATE(B_mat(NUM_NODE, width))

                DO l=1,NUM_ELEM
                ! PRINT *,l
                kappa = mat_properties(1, material_list(l))
                c_rho = mat_properties(2, material_list(l))
                m = mat_properties(3, material_list(l))

                CALL GET_TRI_DELTAS(x_elem(1:3, l), y_elem(1:3, l), dx, dy)
                CALL GET_TRI_AREA(x_elem(1:3, l), dy, A)
                CALL GET_A_ELEM_TRI(dx, dy, A, c_rho, kappa, m, theta, dt, elem_mat_3)

                DO i=1,3
                DO j=1,3

                row = elem_list(i, l)
                col = elem_list(j, l)

                CALL mode2_index_map(row, col, bandwidth, new_col)

                A_mat(row, new_col) = A_mat(row, new_col) + elem_mat_3(i,j)

                END DO
                END DO

                CALL GET_B_ELEM_TRI(dx, dy, A, c_rho, kappa, m, theta, dt, elem_mat_3)

                DO i=1,3
                DO j=1,3

                row = elem_list(i, l)
                col = elem_list(j, l)

                CALL mode2_index_map(row, col, bandwidth, new_col)

                B_mat(row, new_col) = B_mat(row, new_col) + elem_mat_3(i,j)

                END DO
                END DO
                END DO

                OPEN(1, file="output/A_1_tri.dat")
                DO i=1,NUM_NODE
                WRITE(1, '( *(g0, ",") )') A_mat(i, :)

                END DO
                CLOSE(1)

                OPEN(1, file="output/B_1_tri.dat")
                DO i=1,NUM_NODE
                WRITE(1, '( *(g0, ",") )') B_mat(i, :)

                END DO
                CLOSE(1)

                DO i=1,NUM_BC

                row = bc_nodes(i)
                A_mat(row, :) = 0
                A_mat(row, bandwidth + 1) = 1

                END DO

                CALL DSOLVE(1, A_mat, U, NUM_NODE, bandwidth, NUM_NODE, 2*bandwidth+1)

                U = 0
                OPEN(1, file="output/u_1_tri_dt.dat")

                DO k=1,TIMESTEPS

                WRITE(1, '( *(g0, ",") )') U

                CALL BANDED_MULT_OUTPUT(B_mat, U, NUM_NODE, 2*bandwidth+1, bandwidth, RHS)

                DO l=1,NUM_ELEM

                CALL GET_RHS_ELEM_TRI(A, heat_rate_list(l), b_vec_3)
                b_vec_3 = b_vec_3 * dt
                DO i=1,3

                row = elem_list(i,l)
                RHS(row) = RHS(row) + b_vec_3(i)

                END DO
                
                IF (k == 1) THEN
                        OPEN(2, file="output/RHS_1_tri_dt.dat")
                        WRITE(2, '( *(g0, "," ) )') RHS
                        CLOSE(2)
                END IF

                DO i=1,NUM_BC

                row = bc_nodes(i)
                RHS(row) = 0

                END DO
                END DO

                CALL DSOLVE(2, A_mat, RHS, NUM_NODE, bandwidth, NUM_NODE, 2*bandwidth + 1)

                U = RHS
                END DO

                WRITE(1, '( *(g0, ",") )') U
                CLOSE(1)

                !CALL DSOLVE(3, LHS, RHS, NUM_NODE, bandwidth, NUM_NODE, 2*bandwidth + 1)

        END PROGRAM

        SUBROUTINE ELEMENT_MAT_RECT(x, y, kappa, m, sigmaE, elem_mat, elem_b)
                REAL*8, DIMENSION(4), INTENT(IN) :: x, y, sigmaE
                REAL*8, INTENT(IN) :: kappa, m 
                REAL*8, DIMENSION(4,4), INTENT(OUT) :: elem_mat
                REAL*8, DIMENSION(4), INTENT(OUT) :: elem_b

                REAL*8, DIMENSION(2,2) :: jac, jac_inv
                REAL*8, DIMENSION(2,4) :: quad_points

                REAL*8 :: jac_det, xi, nu, phi_i, phi_j, dphi_i_x, dphi_i_y, dphi_j_x, dphi_j_y 
                INTEGER :: i, j, k

                quad_points = RESHAPE((/ -1/SQRT(3.0), -1/SQRT(3.0), 1/SQRT(3.0), -1/SQRT(3.0), 1/SQRT(3.0), 1/SQRT(3.0), 1/SQRT(3.0), -1/SQRT(3.0) /), SHAPE(quad_points))

                elem_mat = 0
                elem_b = 0

                DO i=1,4
                DO j=1,4
                DO k=1,4
                xi = quad_points(1, k)
                nu = quad_points(2, k)

                CALL GET_JAC(x, y, xi, nu, jac, jac_inv, jac_det)
                CALL GET_QUAD(xi, nu, jac_inv, i, phi_i, dphi_i_x, dphi_i_y)
                CALL GET_QUAD(xi, nu, jac_inv, j, phi_j, dphi_j_x, dphi_j_y)

                elem_mat(i,j) = elem_mat(i,j) + (kappa * (dphi_i_x * dphi_j_x + dphi_i_y * dphi_j_y) + m * phi_i * phi_j) * jac_det

                elem_b(i) = elem_b(i) + sigmaE(j) * phi_j * phi_i * jac_det

                END DO
                END DO
                END DO

        END SUBROUTINE

        SUBROUTINE GET_QUAD(xi, nu, jac_inv, basis, phi, dphi_x, dphi_y)
                REAL*8, INTENT(IN) :: xi, nu
                REAL*8, DIMENSION(2,2), INTENT(IN) :: jac_inv
                INTEGER, INTENT(IN) :: basis

                REAL*8, INTENT(OUT) :: phi, dphi_x, dphi_y
                REAL*8 :: dphi_xi, dphi_nu

                IF (basis == 1) THEN
                phi = (1 - xi) * (1 - nu) / 4.0
                dphi_xi = -(1 - nu) / 4.0
                dphi_nu = -(1 - xi) / 4.0

                ELSE IF (basis == 2) THEN
                phi = (1 - xi) * (1 + nu) / 4.0
                dphi_xi = -(1 + nu) / 4.0
                dphi_nu = (1 - xi) / 4.0

                ELSE IF (basis == 3) THEN
                phi = (1 + xi) * (1 + nu) / 4.0
                dphi_xi = (1 + nu) / 4.0
                dphi_nu = (1 + xi) / 4.0

                ELSE
                phi = (1 + xi) * (1 - nu) / 4.0
                dphi_xi = (1 - nu) / 4.0
                dphi_nu = -(1 + xi) / 4.0

                END IF

                dphi_x = jac_inv(1,1) * dphi_xi + jac_inv(1,2) * dphi_nu
                dphi_y = jac_inv(2,1) * dphi_xi + jac_inv(2,2) * dphi_nu

        END SUBROUTINE

        SUBROUTINE GET_JAC(x, y, xi, nu, jac, jac_inv, jac_det)
                REAL*8, DIMENSION(4), INTENT(IN) :: x, y
                REAL*8, INTENT(IN) :: xi, nu

                REAL*8, DIMENSION(2,2), INTENT(OUT) :: jac, jac_inv
                REAL*8, INTENT(OUT) :: jac_det

                REAL*8 :: dxdxi, dxdnu, dydxi, dydnu, temp
                INTEGER :: i

                dxdxi = 0
                dxdnu = 0
                dydxi = 0
                dydnu = 0

                DO i=1,4
                CALL XI_DERIV(i, nu, temp)
                dxdxi = dxdxi + x(i) * temp
                dydxi = dydxi + y(i) * temp

                CALL NU_DERIV(i, xi, temp)
                dxdnu = dxdnu + x(i) * temp
                dydnu = dydnu + y(i) * temp
                END DO

                jac(1,1) = dxdxi
                jac(2,1) = dxdnu
                jac(1,2) = dydxi
                jac(2,2) = dydnu

                jac_det = dxdxi * dydnu - dxdnu * dydxi

                jac_inv(1,1) = dydnu / jac_det
                jac_inv(2,1) = -dxdnu / jac_det
                jac_inv(1,2) = -dydxi / jac_det
                jac_inv(2,2) = dxdxi / jac_det

        END SUBROUTINE

        SUBROUTINE XI_DERIV(basis_num, nu, dphidxi)
                INTEGER, INTENT(IN) :: basis_num
                REAL*8, INTENT(IN) :: nu
                REAL*8, INTENT(OUT) :: dphidxi

                IF (basis_num == 1) THEN
                dphidxi = -(1 - nu) / 4.0

                ELSE IF (basis_num == 2) THEN
                dphidxi = -(1 + nu) / 4.0

                ELSE IF (basis_num == 3) THEN
                dphidxi = (1 + nu) / 4.0

                ELSE
                dphidxi = (1 - nu) / 4.0

                END IF

        END SUBROUTINE

        SUBROUTINE NU_DERIV(basis_num, xi, dphidnu)
                INTEGER, INTENT(IN) :: basis_num
                REAL*8, INTENT(IN) :: xi
                REAL*8, INTENT(OUT) :: dphidnu

                IF (basis_num == 1) THEN
                dphidnu = -(1 - xi) / 4.0

                ELSE IF (basis_num == 2) THEN
                dphidnu = (1 - xi) / 4.0

                ELSE IF (basis_num == 3) THEN
                dphidnu = (1 + xi) / 4.0

                ELSE
                dphidnu = -(1 + xi) / 4.0

                END IF

        END SUBROUTINE

        SUBROUTINE ELEMENT_MAT_TRI(dx, dy, A, kappa, m, sigmaE, elem_mat, elem_b)
                REAL*8, DIMENSION(3), INTENT(IN) :: dx, dy, sigmaE
                REAL*8, INTENT(IN) :: kappa, m, A

                REAL*8, DIMENSION(3, 3), INTENT(OUT) :: elem_mat
                REAL*8, DIMENSION(3), INTENT(OUT) :: elem_b

                elem_mat = 0
                elem_b = 0

                elem_mat(1,1) = kappa / (4 * A) * (dx(1) * dx(1) + dy(1) * dy(1)) + m * A / 6.0
                elem_mat(1,2) = kappa / (4 * A) * (dx(1) * dx(2) + dy(1) * dy(2)) + m * A / 12.0
                elem_mat(1,3) = kappa / (4 * A) * (dx(1) * dx(3) + dy(1) * dy(3)) + m * A / 12.0

                elem_mat(2,1) = kappa / (4 * A) * (dx(2) * dx(1) + dy(2) * dy(1)) + m * A / 12.0
                elem_mat(2,2) = kappa / (4 * A) * (dx(2) * dx(2) + dy(2) * dy(2)) + m * A / 6.0
                elem_mat(2,3) = kappa / (4 * A) * (dx(2) * dx(3) + dy(2) * dy(3)) + m * A / 12.0

                elem_mat(3,1) = kappa / (4 * A) * (dx(3) * dx(1) + dy(3) * dy(1)) + m * A / 12.0
                elem_mat(3,2) = kappa / (4 * A) * (dx(3) * dx(2) + dy(3) * dy(2)) + m * A / 12.0
                elem_mat(3,3) = kappa / (4 * A) * (dx(3) * dx(3) + dy(3) * dy(3)) + m * A / 6.0

                elem_b(1) = (2 * sigmaE(1) + sigmaE(2) + sigmaE(3)) * A / 12.0
                elem_b(2) = (sigmaE(1) + 2 * sigmaE(2) + sigmaE(3)) * A / 12.0
                elem_b(3) = (sigmaE(1) + sigmaE(2) + 2 * sigmaE(3)) * A / 12.0 

        END SUBROUTINE

        SUBROUTINE GET_TRI_AREA(x, dy, A)
                REAL*8, DIMENSION(3), INTENT(IN) :: x, dy
                REAL*8, INTENT(OUT) :: A

                A = 0.5 * (x(1) * dy(1) + x(2) * dy(2) + x(3) * dy(3))
                END SUBROUTINE

                SUBROUTINE GET_TRI_DELTAS(x, y, dx, dy)
                REAL*8, DIMENSION(3), INTENT(IN) :: x, y
                REAL*8, DIMENSION(3), INTENT(OUT) :: dx, dy

                dx(1) = x(2) - x(3)
                dx(2) = x(3) - x(1)
                dx(3) = x(1) - x(2)

                dy(1) = y(2) - y(3)
                dy(2) = y(3) - y(1)
                dy(3) = y(1) - y(2)

        END SUBROUTINE

        SUBROUTINE GET_M_ELEM_TRI(A, c_rho, M_elem)
                REAL*8, INTENT(IN) :: A, c_rho
                REAL*8, DIMENSION(3,3), INTENT(OUT) :: M_elem
                INTEGER :: i, j

                M_elem = 0

                DO i=1,3
                DO j=1,3

                IF (i == j) THEN
                M_elem(i,j) = c_rho * A / 6.0
                ELSE 
                M_elem(i,j) = c_rho * A / 12.0
                END IF

                END DO
                END DO

        END SUBROUTINE

        SUBROUTINE GET_K_ELEM_TRI(dx, dy, A, kappa, m, K_elem)
                REAL*8, DIMENSION(3), INTENT(IN) :: dx, dy
                REAL*8, INTENT(IN) :: A, kappa, m 
                REAL*8, DIMENSION(3,3), INTENT(OUT) :: K_elem
                INTEGER :: i, j

                K_elem = 0

                DO i=1,3
                DO j=1,3

                K_elem(i,j) = kappa * (dx(i) * dx(j) + dy(i) * dy(j)) / (4 * A)

                IF (i == j) THEN
                K_elem(i,j) = K_elem(i,j) + m * A / 6.0

                ELSE 
                K_elem(i,j) = K_elem(i,j) + m * A / 12.0

                END IF

                END DO
                END DO

        END SUBROUTINE

        SUBROUTINE GET_A_ELEM_TRI(dx, dy, A, c_rho, kappa, m, theta, dt, A_elem)
                REAL*8, DIMENSION(3), INTENT(IN) :: dx, dy
                REAL*8, INTENT(IN) :: A, c_rho, kappa, m, theta, dt
                REAL*8, DIMENSION(3,3), INTENT(OUT) :: A_elem

                REAL*8, DIMENSION(3,3) :: temp_mat

                CALL GET_M_ELEM_TRI(A, c_rho, temp_mat)
                A_elem = temp_mat

                CALL GET_K_ELEM_TRI(dx, dy, A, kappa, m, temp_mat)
                A_elem = A_elem + theta * dt * temp_mat

        END SUBROUTINE

        SUBROUTINE GET_B_ELEM_TRI(dx, dy, A, c_rho, kappa, m, theta, dt, B_elem)
                REAL*8, DIMENSION(3), INTENT(IN) :: dx, dy
                REAL*8, INTENT(IN) :: A, c_rho, kappa, m, theta, dt
                REAL*8, DIMENSION(3,3), INTENT(OUT) :: B_elem

                REAL*8, DIMENSION(3,3) :: temp_mat

                CALL GET_M_ELEM_TRI(A, c_rho, temp_mat)
                B_elem = temp_mat

                CALL GET_K_ELEM_TRI(dx, dy, A, kappa, m, temp_mat)
                B_elem = B_elem - (1 - theta) * dt * temp_mat

        END SUBROUTINE

        SUBROUTINE GET_RHS_ELEM_TRI(A, sigmaE, rhs_elem)
                REAL*8, INTENT(IN) :: A
                REAL*8, INTENT(IN) :: sigmaE
                REAL*8, DIMENSION(3), INTENT(OUT) :: rhs_elem

                INTEGER :: i

                rhs_elem = sigmaE * A / 3.0

        END SUBROUTINE
