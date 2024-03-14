        PROGRAM transient
                INTEGER, PARAMETER :: NUM_NODE=969, NUM_ELEM=1089, NUM_BC=64, NUM_MAT=8, NUM_STEPS=250
                REAL*8, PARAMETER :: dt=200, theta=0.5
                REAL*8, DIMENSION(3, NUM_NODE) :: node_input
                REAL*8, DIMENSION(2, NUM_NODE) :: node

                REAL*8, DIMENSION(3, NUM_NODE) :: sigma_input
                REAL*8, DIMENSION(NUM_NODE) :: sigma

                INTEGER, DIMENSION(6, NUM_ELEM) :: elem_input 
                INTEGER, DIMENSION(4, NUM_ELEM) :: elem
                INTEGER, DIMENSION(NUM_ELEM) :: elem_material

                REAL*8, DIMENSION(4, NUM_MAT) :: material_input
                REAL*8, DIMENSION(3, NUM_MAT) :: material

                REAL*8, DIMENSION(7, NUM_BC) :: bc_input
                INTEGER, DIMENSION(NUM_BC) :: bc_node, bc_next, bc_prev
                REAL*8, DIMENSION(NUM_BC) :: bc_h, bc_Ta

                REAL*8, DIMENSION(:,:), ALLOCATABLE :: A_packed, M_packed, LHS_mat
                REAL*8, DIMENSION(NUM_NODE, NUM_NODE) :: A_unpacked, M_unpacked, RHS_mat
                REAL*8, DIMENSION(NUM_NODE) :: RHS, RHS_cpy
                REAL*8, DIMENSION(NUM_NODE, NUM_STEPS+1) :: u

                REAL*8, DIMENSION(4, NUM_ELEM) :: x, y
                REAL*8, DIMENSION(3, 3) :: A_3, M_3
                REAL*8, DIMENSION(3) :: b_3
                REAL*8, DIMENSION(4, 4) :: A_4, M_4
                REAL*8, DIMENSION(4) :: b_4

                INTEGER, DIMENSION(NUM_NODE) :: IPIV
                INTEGER :: INFO

                REAL*8 :: kappa, m, rho, h, Ta, L1, L2
                INTEGER :: i, j, k, l, width, bandwidth
                INTEGER :: row, col, new_row, new_col, first3, first4

                first3 = 0
                first4 = 0
                IPIV = 0

                OPEN(1, file="npeltr4.dat")
                READ(1, *) node_input
                CLOSE(1)

                node = node_input(2:3, :)

                OPEN(1, file="ppelt4.dat")
                READ(1, *) sigma_input

                sigma = sigma_input(3, :)

                OPEN(1, file="epeltr4.dat")
                READ(1, *) elem_input
                CLOSE(1)

                elem = elem_input(2:5, :)
                elem_material = elem_input(6, :)

                OPEN(1, file="bpeltr4.dat")
                READ(1, *) bc_input
                CLOSE(1)
                
                bc_node = INT(bc_input(2, :))
                bc_next = INT(bc_input(4, :))
                bc_prev = INT(bc_input(5, :))
                bc_h = bc_input(6, :)
                bc_Ta = bc_input(7, :)

                OPEN(1, file="material_properties.dat")
                READ(1, *) material_input
                CLOSE(1)

                material = material_input(2:4, :);

                bandwidth = -1

                DO l=1,NUM_ELEM
                DO i=1,4

                        x(i, l) = node(1, elem(i, l))
                        y(i, l) = node(2, elem(i, l))

                DO j=i+1,4

                        width = ABS(elem(i, l) - elem(j, l))

                        IF (width > bandwidth) THEN
                                bandwidth = width
                        END IF

                END DO
                END DO
                END DO

                ALLOCATE(A_packed(NUM_NODE, 2*bandwidth+1))
                ALLOCATE(M_packed(NUM_NODE, 2*bandwidth+1))
                ALLOCATE(LHS_mat(NUM_NODE, 2*bandwidth+1))
                A_packed = 0
                M_packed = 0
                A_unpacked = 0
                M_unpacked = 0
                RHS = 0

                DO l=1,NUM_ELEM

                kappa = material(1, elem_material(l))
                rho = material(2, elem_material(l))
                m = material(3, elem_material(l))

                IF (elem(3, l) == elem(4, l)) THEN

                        CALL ELEM_MAT_3(x(1:3, l), y(1:3, l), kappa, m, rho, sigma(l), A_3, M_3, b_3)

                DO i=1,3
                DO j=1,3

                        row = elem(i, l)
                        col = elem(j, l)

                        A_unpacked(row, col) = A_unpacked(row, col) + A_3(i, j) 
                        M_unpacked(row, col) = M_unpacked(row, col) + M_3(i, j)

                        CALL MODE_2_INDEX_MAP(row, col, bandwidth, new_col)
                        A_packed(row, new_col) = A_packed(row, new_col) + A_3(i, j)
                        M_packed(row, new_col) = M_packed(row, new_col) + M_3(i, j)

                END DO
                END DO

                DO i=1,3

                        row = elem(i, l)
                        RHS(row) = RHS(row) + b_3(i)

                END DO

                ELSE

                        CALL ELEM_MAT_4(x(:, l), y(:, l), kappa, m, rho, sigma(l), A_4, M_4, b_4)

                DO i=1,4
                DO j=1,4

                        row = elem(i, l)
                        col = elem(j, l)

                        A_unpacked(row, col) = A_unpacked(row, col) + A_4(i, j)
                        M_unpacked(row, col) = M_unpacked(row, col) + M_4(i, j)

                        CALL MODE_2_INDEX_MAP(row, col, bandwidth, new_col)
                        A_packed(row, new_col) = A_packed(row, new_col) + A_4(i, j)
                        M_packed(row, new_col) = M_packed(row, new_col) + M_4(i, j)

                END DO
                END DO

                DO i=1,4

                        row = elem(i, l)
                        RHS(row) = RHS(row) + b_4(i)

                END DO

                END IF

                END DO

                DO i=1,NUM_BC
                        h = bc_h(i)
                        Ta = bc_Ta(i)

                        CALL GET_LENGTH(x(:, bc_prev(i)), y(:, bc_prev(i)), x(:, bc_node(i)), y(:, bc_node(i)), L1)
                        CALL GET_LENGTH(x(:, bc_next(i)), y(:, bc_next(i)), x(:, bc_node(i)), y(:, bc_node(i)), L2)

                        row = bc_node(i)
                        col = bc_node(i)

                        RHS(row) = RHS(row) + h * Ta * (L1 + L2) / 2.0

                        CALL MODE_2_INDEX_MAP(row, col, bandwidth, new_col)

                        A_packed(row, new_col) = A_packed(row, new_col) + h * (L1 + L2) / 3.0
                        A_unpacked(row, col) = A_unpacked(row, col) + h * (L1 + L2) / 3.0

                        col = bc_prev(i)

                        CALL MODE_2_INDEX_MAP(row, col, bandwidth, new_col)

                        A_packed(row, new_col) = A_packed(row, new_col) + h * L1 / 6.0
                        A_unpacked(row, col) = A_unpacked(row, col) + h * L1 / 6.0

                        col = bc_next(i)

                        CALL MODE_2_INDEX_MAP(row, col, bandwidth, new_col)

                        A_packed(row, new_col) = A_packed(row, new_col) + h * L2 / 6.0
                        A_unpacked(row, col) = A_unpacked(row, col) + h * L2 / 6.0

                END DO

                LHS_mat = 0
                RHS_mat = 0
                
                LHS_mat = M_packed + theta * dt * A_packed
                RHS_mat = M_unpacked - (1 - theta) * dt * A_unpacked
                RHS = RHS * dt
                RHS_cpy = RHS

                CALL DSOLVE(1, LHS_mat, RHS, NUM_NODE, bandwidth, NUM_NODE, 2*bandwidth+1)

                DO k=1,NUM_STEPS

                        RHS = RHS_cpy + MATMUL(RHS_mat, u(:, k))

                        CALL DSOLVE(2, LHS_mat, RHS, NUM_NODE, bandwidth, NUM_NODE, 2*bandwith+1)

                        u(:, k+1) = RHS

                END DO

                OPEN(1, file="output/u_transient.dat")
                DO k=1,NUM_STEPS+1

                        WRITE(1, '(*(g0, ","))') u(:, k)

                END DO
                CLOSE(1)

                RHS = RHS_cpy / dt

                CALL DSOLVE(3, A_packed, RHS, NUM_NODE, bandwidth, NUM_NODE, 2*bandwidth+1)

                OPEN(1, file="output/u_test.dat")

                DO i=1,NUM_NODE

                        WRITE(1, '(*(g0, ","))') RHS(i)

                END DO

                CLOSE(1)

        END PROGRAM

        SUBROUTINE ELEM_MAT_3(x, y, kappa, m, rho, sigma, A_3, M_3, b_3)
                REAL*8, DIMENSION(3), INTENT(IN) :: x ,y
                REAL*8, INTENT(IN) :: kappa, m, rho, sigma
                REAL*8, DIMENSION(3,3), INTENT(OUT) :: A_3, M_3
                REAL*8, DIMENSION(3), INTENT(OUT) :: b_3

                REAL*8, DIMENSION(3) :: dx, dy
                REAL*8 :: A
                INTEGER :: i, j

                CALL GET_DELTAS(x, y, dx, dy)
                CALL GET_AREA(x, dy, A)

                DO i=1,3
                DO j=1,3

                IF (i == j) THEN

                        A_3(i, j) = kappa / (4 * A) * (dx(i)*dx(j) + dy(i)*dy(j)) + m * A / 6.0

                        M_3(i, j) = rho * A / 6.0

                ELSE

                        A_3(i, j) = kappa / (4 * A) * (dx(i)*dx(j) + dy(i)*dy(j)) + m * A / 12.0

                        M_3(i, j) = rho * A / 12.0

                END IF

                END DO
                END DO

                DO i=1,3
                
                        b_3(i) = sigma * A / 3.0

                END DO

        END SUBROUTINE

        SUBROUTINE GET_DELTAS(x, y, dx, dy)
                REAL*8, DIMENSION(3), INTENT(IN) :: x, y
                REAL*8, DIMENSION(3), INTENT(OUT) :: dx, dy

                dx(1) = x(2) - x(3)
                dx(2) = x(3) - x(1)
                dx(3) = x(1) - x(2)

                dy(1) = y(2) - y(3)
                dy(2) = y(3) - y(1)
                dy(3) = y(1) - y(2)

        END SUBROUTINE

        SUBROUTINE GET_AREA(x, dy, A)
                REAL*8, DIMENSION(3), INTENT(IN) :: x, dy
                REAL*8, INTENT(OUT) :: A

                INTEGER :: i

                A = 0

                DO i=1,3
                        A = A + x(i) * dy(i)
                END DO

                A = A / 2.0

        END SUBROUTINE

        SUBROUTINE ELEM_MAT_4(x, y, kappa, m, rho, sigma, A_4, M_4, b_4)
                REAL*8, DIMENSION(4), INTENT(IN) :: x, y
                REAL*8, INTENT(IN) :: kappa, m, rho, sigma
                REAL*8, DIMENSION(4,4), INTENT(OUT) :: A_4, M_4
                REAL*8, DIMENSION(4), INTENT(OUT) :: b_4

                REAL*8, DIMENSION(4) :: zeta_k, eta_k
                REAL*8 :: phi_i, phi_j
                REAL*8 :: dx_dzeta, dx_deta, dy_dzeta, dy_deta, det
                REAL*8 :: dphi_i_dzeta, dphi_i_deta, dphi_j_dzeta, dphi_j_deta
                REAL*8 :: dphi_i_dx, dphi_i_dy, dphi_j_dx, dphi_j_dy
                INTEGER :: i, j, k
                
                zeta_k = (/-1, -1, 1, 1/) * 1 / SQRT(3.0)
                eta_k = (/-1, 1, -1, 1/) * 1 / SQRT(3.0)

                A_4 = 0
                M_4 = 0
                b_4 = 0

                DO i=1,4
                DO j=1,4
                DO k=1,4

                        CALL PHI(i, zeta_k(k), eta_k(k), phi_i)
                        CALL PHI(j, zeta_k(k), eta_k(k), phi_j)

                        CALL DPHI_LOCAL(i, zeta_k(k), eta_k(k), dphi_i_dzeta, dphi_i_deta)
                        CALL DPHI_LOCAL(j, zeta_k(k), eta_k(k), dphi_j_dzeta, dphi_j_deta)

                        CALL JACOBIAN(x, y, zeta_k(k), eta_k(k), dx_dzeta, dx_deta, dy_dzeta, dy_deta, det)

                        CALL DPHI_GLOBAL(dphi_i_dzeta, dphi_i_deta, dx_dzeta, dx_deta, dy_dzeta, dy_deta, dphi_i_dx, dphi_i_dy, det)
                        CALL DPHI_GLOBAL(dphi_j_dzeta, dphi_j_deta, dx_dzeta, dx_deta, dy_dzeta, dy_deta, dphi_j_dx, dphi_j_dy, det)

                        A_4(i, j) = A_4(i, j) + (kappa * (dphi_i_dx * dphi_j_dx + dphi_i_dy * dphi_j_dy) + m * phi_i * phi_j) * det

                        M_4(i, j) = M_4(i, j) + rho * phi_i * phi_j * det

                END DO
                END DO
                END DO

                DO i=1,4
                DO k=1,4

                        CALL PHI(i, zeta_k(k), eta_k(k), phi_i)
                        CALL JACOBIAN(x, y, zeta_k(k), eta_k(k), dx_dzeta, dx_deta, dy_dzeta, dy_deta, det)

                        b_4(i) = b_4(i) + sigma * phi_i * det

                END DO
                END DO

        END SUBROUTINE

        SUBROUTINE JACOBIAN(x, y, zeta, eta, dx_dzeta, dx_deta, dy_dzeta, dy_deta, det)
                REAL*8, DIMENSION(4), INTENT(IN) :: x, y
                REAL*8, INTENT(IN) :: zeta, eta
                REAL*8, INTENT(OUT) :: dx_dzeta, dx_deta, dy_dzeta, dy_deta, det
                
                REAL*8 :: dphi_dzeta, dphi_deta
                INTEGER :: l

                dx_dzeta = 0
                dx_deta = 0
                dy_dzeta = 0
                dy_deta = 0

                DO l=1,4

                        CALL DPHI_LOCAL(l, zeta, eta, dphi_dzeta, dphi_deta)

                        dx_dzeta = dx_dzeta + x(l) * dphi_dzeta
                        dx_deta = dx_deta + x(l) * dphi_deta

                        dy_dzeta = dy_dzeta + y(l) * dphi_dzeta
                        dy_deta = dy_deta + y(l) * dphi_deta

                END DO

                det = dx_dzeta * dy_deta - dx_deta * dy_dzeta

        END SUBROUTINE

        SUBROUTINE PHI(l, zeta, eta, phi_l)
                INTEGER, INTENT(IN) :: l
                REAL*8, INTENT(IN) :: zeta, eta
                REAL*8, INTENT(OUT) :: phi_l

                REAL*8, DIMENSION(4) :: zeta_l, eta_l

                zeta_l = (/-1, 1, 1, -1/)
                eta_l = (/-1, -1, 1, 1/)

                phi_l = (1 + zeta_l(l) * zeta) * (1 + eta_l(l) * eta) / 4.0

        END SUBROUTINE

        SUBROUTINE DPHI_LOCAL(l, zeta, eta, dphi_dzeta, dphi_deta)
                INTEGER, INTENT(IN) :: l
                REAL*8, INTENT(IN) :: zeta, eta
                REAL*8, INTENT(OUT) :: dphi_dzeta, dphi_deta

                REAL*8, DIMENSION(4) :: zeta_l, eta_l

                zeta_l = (/-1, 1, 1, -1/)
                eta_l = (/-1, -1, 1, 1/)

                dphi_dzeta = zeta_l(l) * (1 + eta_l(l) * eta) / 4.0
                dphi_deta = eta_l(l) * (1 + zeta_l(l) * zeta) / 4.0

        END SUBROUTINE

        SUBROUTINE DPHI_GLOBAL(dphi_dzeta, dphi_deta, dx_dzeta, dx_deta, dy_dzeta, dy_deta, dphi_dx, dphi_dy, det)
                REAL*8, INTENT(IN) :: dphi_dzeta, dphi_deta, dx_dzeta, dx_deta, dy_dzeta, dy_deta, det
                REAL*8, INTENT(OUT) :: dphi_dx, dphi_dy

                dphi_dx = (dy_deta * dphi_dzeta - dy_dzeta * dphi_deta) / det
                dphi_dy = (dx_dzeta * dphi_deta - dx_deta * dphi_dzeta) / det

        END SUBROUTINE

        SUBROUTINE GET_LENGTH(x1, y1, x2, y2, L)
                REAL*8, DIMENSION(3), INTENT(IN) :: x1, y1, x2, y2
                REAL*8, INTENT(OUT) :: L

                L = 0

                DO i=1,3

                        L = L + (x1(i) - x2(i))**2 + (y1(i) - y2(i))**2

                END DO

                L = SQRT(L)

        END SUBROUTINE
