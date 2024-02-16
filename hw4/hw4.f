        PROGRAM hw4
                INTEGER, PARAMETER :: NUM_NODES = 28, NUM_ELEMS = 38, PART = 1
                INTEGER, PARAMETER :: NUM_BC_NODES = 16, BANDWIDTH = 5
                REAL*8, DIMENSION(2, NUM_NODES) :: node_pos
                REAL*8, DIMENSION(3, NUM_NODES) :: node_pos_input

                INTEGER, DIMENSION(3, NUM_ELEMS) :: elem_tab
                INTEGER, DIMENSION(5, NUM_ELEMS) :: elem_tab_input
                REAL*8, DIMENSION(3, NUM_ELEMS) :: x_elem, y_elem, dx_elem, dy_elem
                REAL*8, DIMENSION(NUM_ELEMS) :: A_elem

                REAL*8, DIMENSION(NUM_NODES, NUM_NODES) :: coeff
                REAL*8, DIMENSION(NUM_NODES, 2*BANDWIDTH+1) :: band_coeff
                REAL*8, DIMENSION(NUM_NODES) :: RHS
                REAL*8, DIMENSION(3,3) :: elem_mat

                REAL*8, DIMENSION(4, NUM_BC_NODES) :: bc_mat
                
                REAL*8 :: garbage1, garbage2, u_coeff, du_coeff, const
                INTEGER :: i, j, k, garbage3

                OPEN(1, file="hw4.ele1")
                READ(1,*) elem_tab_input
                CLOSE(1)
                elem_tab = elem_tab_input(2:4, :)

                OPEN(1, file="hw4.nod")
                READ(1,*) node_pos_input
                CLOSE(1)
                node_pos = node_pos_input(2:3, :)

                IF (PART == 1) THEN
                        OPEN(1, file="hw4a.bc")
                        READ(1, *) bc_mat
                ELSE IF (PART == 2) THEN
                        OPEN(1, file="hw4b.bc")
                        READ(1, *) bc_mat

                ELSE 
                        OPEN(1, file="hw4c.bc")
                        READ(1, *) bc_mat
                END IF
                CLOSE(1)

                coeff = 0
                DO i=1,NUM_ELEMS
                        x_elem(1,i) = node_pos(1,elem_tab(1,i))
                        x_elem(2,i) = node_pos(1,elem_tab(2,i))
                        x_elem(3,i) = node_pos(1,elem_tab(3,i))

                        y_elem(1,i) = node_pos(2,elem_tab(1,i))
                        y_elem(2,i) = node_pos(2,elem_tab(2,i))
                        y_elem(3,i) = node_pos(2,elem_tab(3,i))

                        CALL GET_DELTAS_VEC(x_elem(:,i), y_elem(:,i), dx_elem(:,i), dy_elem(:,i))

                        CALL GET_AREA_VEC(x_elem(:,i), dy_elem(:,i), A_elem(i))

                        CALL ELEMENT_MATRIX_VEC(dx_elem(:,i), dy_elem(:,i), A_elem(i), elem_mat)

                        DO j=1,3
                        DO k=1,3

                        coeff(elem_tab(j, i), elem_tab(k, i)) = coeff(elem_tab(j, i), elem_tab(k, i)) + elem_mat(j,k)

                        END DO
                        END DO
                        
                END DO

                OPEN(1, file="output/coeff.dat")
                OPEN(2, file="output/rhs.dat")

                RHS = 0
                DO i=1,NUM_NODES
                        CALL GET_BC(i, NUM_BC_NODES, bc_mat, u_coeff, du_coeff, const)

                        ! Interior node
                        IF (u_coeff == 0 .AND. du_coeff == 0 .AND. const == 0) THEN
                                RHS(i) = RHS(i)

                        ! Type I BC
                        ELSE IF (du_coeff == 0) THEN
                                RHS(i) = const
                                coeff(i, :) = 0
                                coeff(i, i) = u_coeff
                        
                        ! Type II or III BC
                        ELSE
                                RHS(i) = const / du_coeff
                                coeff(i, :) = -u_coeff / du_coeff
                        END IF

                        WRITE(1, '( *(g0, ",") )') coeff(i,:)
                        WRITE(2, *) RHS(i)
                END DO

                CLOSE(1)
                CLOSE(2)

        END PROGRAM

        SUBROUTINE GET_BC(node, num_bcs, bc_mat, u_coeff, du_coeff, const)
                INTEGER, INTENT(IN) :: node, num_bcs
                REAL*8, DIMENSION(4, num_bcs), INTENT(IN) :: bc_mat
                REAL*8, INTENT(OUT) :: u_coeff, du_coeff, const

                INTEGER :: i

                u_coeff = 0
                du_coeff = 0
                const = 0

                DO i=1,num_bcs
                        IF (node == bc_mat(1, i)) THEN
                                u_coeff = bc_mat(2, i)
                                du_coeff = bc_mat(3, i)
                                const = bc_mat(4, i)
                        END IF
                END DO

        END SUBROUTINE

        SUBROUTINE ELEMENT_MATRIX(dx_1, dx_2, dx_3, dy_1, dy_2, dy_3, A, Ae)
                REAL*8, INTENT(IN) :: dx_1, dx_2, dx_3, dy_1, dy_2, dy_3, A
                REAL*8, DIMENSION(3,3), INTENT(OUT) :: Ae

                Ae(1,1) = -(dy_1 * dy_1 + dx_1 * dx_1) / (4 * A)
                Ae(1,2) = -(dy_1 * dy_2 + dx_1 * dx_2) / (4 * A)
                Ae(1,3) = -(dy_1 * dy_3 + dx_1 * dx_3) / (4 * A)

                Ae(2,1) = -(dy_2 * dy_1 + dx_2 * dx_1) / (4 * A)
                Ae(2,2) = -(dy_2 * dy_2 + dx_2 * dx_2) / (4 * A)
                Ae(2,3) = -(dy_2 * dy_3 + dx_2 * dx_3) / (4 * A)

                Ae(3,1) = -(dy_3 * dy_1 + dx_3 * dx_1) / (4 * A)
                Ae(3,2) = -(dy_3 * dy_2 + dx_3 * dx_2) / (4 * A)
                Ae(3,3) = -(dy_3 * dy_3 + dx_3 * dx_3) / (4 * A)

        END SUBROUTINE

        SUBROUTINE ELEMENT_MATRIX_VEC(dx, dy, A, Ae)
                REAL*8, DIMENSION(3), INTENT(IN) :: dx, dy
                REAL*8, INTENT(IN) :: A
                REAL*8, DIMENSION(3,3), INTENT(OUT) :: Ae

                Ae(1,1) = -(dy(1) * dy(1) + dx(1) * dx(1)) / (4 * A)
                Ae(1,2) = -(dy(1) * dy(2) + dx(1) * dx(2)) / (4 * A)
                Ae(1,3) = -(dy(1) * dy(3) + dx(1) * dx(3)) / (4 * A)

                Ae(2,1) = -(dy(2) * dy(1) + dx(2) * dx(1)) / (4 * A)
                Ae(2,2) = -(dy(2) * dy(2) + dx(2) * dx(2)) / (4 * A)
                Ae(2,3) = -(dy(2) * dy(3) + dx(2) * dx(3)) / (4 * A)

                Ae(3,1) = -(dy(3) * dy(1) + dx(3) * dx(1)) / (4 * A)
                Ae(3,2) = -(dy(3) * dy(2) + dx(3) * dx(2)) / (4 * A)
                Ae(3,3) = -(dy(3) * dy(3) + dx(3) * dx(3)) / (4 * A)

        END SUBROUTINE

        SUBROUTINE GET_AREA(x_1, x_2, x_3, dy_1, dy_2, dy_3, A)
                REAL*8, INTENT(IN) :: x_1, x_2, x_3, dy_1, dy_2, dy_3
                REAL*8, INTENT(OUT) :: A

                A = 0.5 * (x_1 * dy_1 + x_2 * dy_2 + x_3 * dy_3)

        END SUBROUTINE

        SUBROUTINE GET_DELTAS(x_1, x_2, x_3, y_1, y_2, y_3, dx, dy)
                REAL*8, INTENT(IN) :: x_1, x_2, x_3, y_1, y_2, y_3
                REAL*8, DIMENSION(3), INTENT(OUT) :: dx, dy

                dx(1) = x_2 - x_3
                dx(2) = x_3 - x_1
                dx(3) = x_1 - x_2

                dy(1) = y_2 - y_3
                dy(2) = y_3 - y_1
                dy(3) = y_1 - y_2

        END SUBROUTINE

        SUBROUTINE GET_AREA_VEC(x, dy, A)
                REAL*8, DIMENSION(3), INTENT(IN) :: x, dy
                REAL*8, INTENT(OUT) :: A

                A = 0.5 * (x(1) * dy(1) + x(2) * dy(2) + x(3) * dy(3))
        
        END SUBROUTINE

        SUBROUTINE GET_DELTAS_VEC(x, y, dx, dy)
                REAL*8, DIMENSION(3), INTENT(IN) :: x, y
                REAL*8, DIMENSION(3), INTENT(OUT) :: dx, dy

                dx(1) = x(2) - x(3)
                dx(2) = x(3) - x(1)
                dx(3) = x(1) - x(2)

                dy(1) = y(2) - y(3)
                dy(2) = y(3) - y(1)
                dy(3) = y(1) - y(2)

        END SUBROUTINE
