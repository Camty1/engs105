        PROGRAM hw6
                INTEGER, PARAMETER :: NUM_NODE=505, NUM_ELEM=917, NUM_BC=33
                REAL*8, DIMENSION(:, :), ALLOCATABLE :: LHS
                REAL*8, DIMENSION(NUM_NODE) :: RHS

                INTEGER, DIMENSION(5, NUM_ELEM) :: elem_input
                INTEGER, DIMENSION(3, NUM_ELEM) :: elem_list

                REAL*8, DIMENSION(3, NUM_NODE) :: node_input 
                REAL*8, DIMENSION(2, NUM_NODE) :: node

                REAL*8, DIMENSION(2, NUM_BC) :: bc

                REAL*8, DIMENSION(3) :: phi_point_source
                
                REAL*8, DIMENSION(3, NUM_ELEM) :: x, y, dx, dy
                REAL*8, DIMENSION(NUM_ELEM) :: A

                REAL*8, DIMENSION(3,3) :: A_elem

                INTEGER :: i, j, k, l, width, bandwidth

                OPEN(1, file="hw44.ele")
                READ(1, *) elem_input
                CLOSE(1)
                
                elem_list = elem_input(2:4, :)

                OPEN(1, file="hw44.nod")
                READ(1, *) node_input
                CLOSE(1)

                node = node_input(2:3, :)

                OPEN(1, file="hw44.dnd")
                READ(1, *) bc
                CLOSE(1)

                bandwidth = 0

                DO l=1,NUM_ELEM
                DO i=1,3
                        x(i, l) = node(1, elem_list(i, l))
                        y(i, l) = node(2, elem_list(i, l))
                        CALL GET_DELTAS(x(:, l), y(:, l), dx(:, l), dy(:, l))
                        CALL GET_AREA(x(:, l), dy(:, l), A(l))

                DO j=i+1,3
                        if (ABS(elem_list(i, l) - elem_list_(j, l)) > bandwidth) THEN

                        bandwidth = ABS(elem_list(i, l) - elem_list_(j, l))
                        END IF

                END DO
                END DO
                END DO

                PRINT *, bandwidth

        END PROGRAM

        SUBROUTINE GET_ELEM_MAT(x, y, dx, dy, A, A_elem)

        END SUBROUTINE

        SUBROUTINE GET_DELTAS(x, y, dx, dy)

        END SUBROUTINE

        SUBROUTINE GET_AREA(dx, y, A)

        END SUBROUTINE
