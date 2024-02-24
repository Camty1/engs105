      PROGRAM HW5
            INTEGER, PARAMETER :: NUM_ELEM=1872, NUM_NODE=969, NUM_BC=64, NUM_MAT=8

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

            REAL*8, DIMENSION(3, NUM_NODE) :: heat_rate_list_input
            REAL*8, DIMENSION(NUM_NODE) :: heat_rate_list
            
            REAL*8, DIMENSION(4, NUM_MAT - 2) :: mat_properties_input
            REAL*8, DIMENSION(3, NUM_MAT) :: mat_properties

            REAL*8, DIMENSION(:,:), ALLOCATABLE :: LHS
            REAL*8, DIMENSION(NUM_NODE) :: RHS

            REAL*8, DIMENSION(3,3) :: elem_mat_3
            REAL*8, DIMENSION(4,4) :: elem_mat_4
            REAL*8, DIMENSION(4) :: tri_quad_points
            REAL*8, DIMENSION(3) :: b_vec_3
            REAL*8, DIMENSION(4) :: b_vec_4

            REAL*8, DIMENSION(4, NUM_ELEM) :: x_elem, y_elem, sigmaE
            REAL*8, DIMENSION(3) :: dx, dy
            REAL*8 :: A, kappa, m, h, ta, L1, L2

            INTEGER :: i, j, k, l, bandwidth, width, row, col, new_col, node, prev_node, next_node

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
            bc_prev_nodes = INT(bc_list_input(5, :))
            bc_next_nodes = INT(bc_list_input(4, :))
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

            sigmaE(1,i) = heat_rate_list(elem_list(1,i))
            sigmaE(2,i) = heat_rate_list(elem_list(2,i))
            sigmaE(3,i) = heat_rate_list(elem_list(3,i))
            sigmaE(4,i) = heat_rate_list(elem_list(4,i))

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
            
            ALLOCATE(LHS(NUM_NODE, width))

            DO l=1,NUM_ELEM
                  ! PRINT *,l
            kappa = mat_properties(1, material_list(l))
            m = mat_properties(3, material_list(l))

            ! Tri
            IF (elem_list(3,l) == elem_list(4,l)) THEN

            CALL GET_TRI_DELTAS(x_elem(1:3, l), y_elem(1:3, l), dx, dy)
            CALL GET_TRI_AREA(x_elem(1:3, l), dy, A)
            CALL ELEMENT_MAT_TRI(dx, dy, A, kappa, m, sigmaE(1:3, l), elem_mat_3, b_vec_3)

            DO i=1,3
            DO j=1,3

            row = elem_list(i, l)
            col = elem_list(j, l)

            CALL mode2_index_map(row, col, bandwidth, new_col)

            LHS(row, new_col) = LHS(row, new_col) + elem_mat_3(i,j)
            END DO

            RHS(row) = RHS(row) + b_vec_3(i)
            END DO

            ELSE

            CALL ELEMENT_MAT_RECT(x_elem(:,l), y_elem(:,l), kappa, m, sigmaE(:, l), elem_mat_4, b_vec_4)

            DO i=1,4
            DO j=1,4

            row = elem_list(i, l)
            col = elem_list(j, l)

            CALL mode2_index_map(row, col, bandwidth, new_col)

            LHS(row, new_col) = LHS(row, new_col) + elem_mat_4(i,j)
            END DO

            RHS(row) = RHS(row) + b_vec_4(i)

            END DO
            END IF
            END DO

            DO i=1,NUM_BC

            node = bc_nodes(i)
            prev_node = bc_prev_nodes(i)
            next_node = bc_next_nodes(i)
            h = bc_h(i)
            ta = bc_ta(i)
            
            L1 = sqrt((node_pos(1, node) - node_pos(1, prev_node)) ** 2 + (node_pos(2, node) - node_pos(2, prev_node)) ** 2)
            L2 = sqrt((node_pos(1, node) - node_pos(1, next_node)) ** 2 + (node_pos(2, node) - node_pos(2, next_node)) ** 2)
            

            row = node
            col = node

            CALL MODE2_INDEX_MAP(row, col, bandwidth, new_col)

            LHS(row, new_col) = LHS(row, new_col) + 2.0 * h / 6.0 * (L1 + L2)

            col = prev_node
            CALL MODE2_INDEX_MAP(row, col, bandwidth, new_col)

            LHS(row, new_col) = LHS(row, new_col) + h / 6.0 * L1

            col = next_node
            CALL MODE2_INDEX_MAP(row, col, bandwidth, new_col)

            LHS(row, new_col) = LHS(row, new_col) + h / 6.0 * L2
            
            RHS(row) = RHS(row) + h * ta / 2.0 * (L1 + L2)

            END DO

            OPEN(1, file="output/LHS_3_tri.dat")
            DO i=1,NUM_NODE
                  WRITE(1, '( *(g0, ",") )') LHS(i, :)

            END DO
            CLOSE(1)

            OPEN(1, file="output/RHS_3_tri.dat")
            DO i=1,NUM_NODE
                  WRITE(1, '( *(g0, ",") )') RHS(i)

            END DO
            CLOSE(1)

            CALL DSOLVE(3, LHS, RHS, NUM_NODE, bandwidth, NUM_NODE, 2*bandwidth + 1)
            
            OPEN(1, file="output/u_3_tri.dat")
            DO i=1,NUM_NODE
                  WRITE(1, '( *(g0, ",") )') RHS(i)

            END DO
            CLOSE(1)
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
