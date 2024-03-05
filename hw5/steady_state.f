        PROGRAM hw5
                INTEGER, PARAMETER :: NUM_NODE=969, NUM_ELEM=1089, NUM_BC=64, NUM_MAT=8
                REAL*8, DIMENSION(3, NUM_NODE) :: node_input
                REAL*8, DIMENSION(2, NUM_NODE) :: node

                REAL*8, DIMENSION(3, NUM_NODE) :: heat_input
                REAL*8, DIMENSION(NUM_NODE) :: heat

                INTEGER, DIMENSION(6, NUM_ELEM) :: elem_input 
                INTEGER, DIMENSION(4, NUM_ELEM) :: elem
                INTEGER, DIMENSION(NUM_ELEM) :: elem_material

                REAL*8, DIMENSION(7, NUM_BC) :: bc_input
                INTEGER, DIMENSION(NUM_BC) :: bc_node, bc_next, bc_prev
                REAL*8, DIMENSION(NUM_BC) :: h, Ta

                REAL*8, DIMENSION(:,:), ALLOCATABLE :: LHS
                REAL*8, DIMENSION(NUM_NODE, NUM_NODE) :: LHS_unpacked
                REAL*8, DIMENSION(NUM_NODE) :: RHS

                REAL*8, DIMENSION(4, NUM_ELEM) :: x, y
                REAL*8, DIMENSION(3, 3) :: A_3
                REAL*8, DIMENSION(4, 4) :: A_4
                INTEGER :: i, j, k, l, width, bandwidth

                OPEN(1, file="npeltr4.dat")
                READ(1, *) node_input
                CLOSE(1)

                node = node_input(2:3, :)

                OPEN(1, file="ppelt4.dat")
                READ(1, *) heat_input

                heat = heat_input(3, :)

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
                h = bc_input(6, :)
                Ta = bc_input(7, :)

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

                ALLOCATE(LHS(3*bandwidth+1, NUM_NODE))
                LHS = 0
                LHS_unpacked = 0
                RHS = 0

                DO l=1,NUM_ELEM

                IF (elem(3, l) == elem(4, l)) THEN

                        CALL ELEM_MAT_3(x(1:3, l), y(1:3, l), kappa, m, A_3)

                DO i=1,3
                DO j=1,3

                        row = elem(i, l)
                        col = elem(j, l)

                        CALL LAPACK_INDEX_MAP(row, col, bandwidth, new_row)

                END DO
                END DO

                ELSE

                END IF

                END DO

        END PROGRAM

        
