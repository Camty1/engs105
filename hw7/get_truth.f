        PROGRAM get_truth
                INTEGER, PARAMETER :: NUM_INPUT_NODE=505, NUM_NODE=502, NUM_ELEM=917, NUM_BC=33, SOURCE_ELEM=288
                REAL*8, PARAMETER :: CURRENT_COEFF=1.0, REMOVE_FACTOR=0.5
                REAL*8, DIMENSION(:, :), ALLOCATABLE :: LHS, LHS_lapack
                REAL*8, DIMENSION(NUM_NODE) :: RHS, RHS_lapack
                REAL*8, DIMENSION(NUM_ELEM, 2) :: current

                INTEGER, DIMENSION(5, NUM_ELEM) :: elem_input
                INTEGER, DIMENSION(3, NUM_ELEM) :: elem_list

                REAL*8, DIMENSION(3, NUM_INPUT_NODE) :: node_input 
                REAL*8, DIMENSION(2, NUM_NODE) :: node

                REAL*8, DIMENSION(2, NUM_BC) :: bc

                REAL*8, DIMENSION(3) :: phi_point_source
                
                REAL*8, DIMENSION(3, NUM_ELEM) :: x, y, dx, dy
                REAL*8, DIMENSION(NUM_ELEM) :: A

                REAL*8, DIMENSION(3,3) :: A_elem

                INTEGER :: i, j, k, l, bandwidth, row, col, new_col

                OPEN(1, file="hw44.ele")
                READ(1, *) elem_input
                CLOSE(1)
                
                elem_list = elem_input(2:4, :)

                OPEN(1, file="hw44.nod")
                READ(1, *) node_input
                CLOSE(1)

                node = node_input(2:3, 1:NUM_NODE)

                OPEN(1, file="hw44.dnd")
                READ(1, *) bc
                CLOSE(1)

                OPEN(1, file="point_source_local_coords.txt")
                READ(1, *) phi_point_source
                CLOSE(1)

                bandwidth = 0

                DO l=1,NUM_ELEM
                DO i=1,3

                        x(i, l) = node(1, elem_list(i, l))
                        y(i, l) = node(2, elem_list(i, l))

                DO j=i+1,3

                        if (ABS(elem_list(i, l) - elem_list(j, l)) > bandwidth) THEN

                                bandwidth = ABS(elem_list(i, l) - elem_list(j, l))

                        END IF

                END DO
                END DO

                        CALL GET_DELTAS(x(:, l), y(:, l), dx(:, l), dy(:, l))
                        CALL GET_AREA(x(:, l), dy(:, l), A(l))

                END DO

                ALLOCATE(LHS(NUM_NODE, 2 * bandwidth + 1))
                ALLOCATE(LHS_lapack(3 * bandwidth + 1, NUM_NODE))

                LHS = 0.0
                RHS = 0.0
                LHS_lapack = 0.0

                DO l=1,NUM_ELEM

                        CALL GET_ELEM_MAT(x(:, l), y(:, l), dx(:, l), dy(:, l), A(l), A_elem)

                        DO i=1,3
                        DO j=1,3

                                row = elem_list(i, l)
                                col = elem_list(j, l)

                                CALL MODE2_INDEX_MAP(row, col, bandwidth, new_col)

                                LHS(row, new_col) = LHS(row, new_col) + A_elem(i,j)

                        END DO
                        END DO

                        IF (l == SOURCE_ELEM) THEN
                        DO i=1,3

                                row = elem_list(i, l)
                                RHS(row) = CURRENT_COEFF * phi_point_source(i)

                        END DO
                        END IF 

                END DO

                RHS(492) = -CURRENT_COEFF * REMOVE_FACTOR / 4.0
                RHS(493) = -CURRENT_COEFF * REMOVE_FACTOR / 4.0

                DO l=1,NUM_BC

                        row = INT(bc(1, l))

                        LHS(row, :) = 0
                        LHS(row, bandwidth + 1) = 1

                        RHS(row) = bc(2, l)

                END DO

                RHS_lapack = RHS

                LHS_lapack(bandwidth + 1:3 * bandwidth + 1, :) = TRANSPOSE(LHS)

                OPEN(1, file="output/LHS.dat")

                DO i=1,NUM_NODE

                        WRITE(1, '( *(g0, ",") )') LHS(i,:)

                END DO

                CLOSE(1)

                OPEN(1, file="output/LHS_lapack.dat")

                DO i=1,NUM_NODE

                        WRITE(1, '( *(g0, ",") )') LHS_lapack(:, i)

                END DO

                CLOSE(1)

                CALL CAM_DSOLVE(bandwidth, LHS_lapack, RHS_lapack, NUM_NODE)

                OPEN(1, file="output/truth.dat")

                DO i=1,NUM_NODE

                        WRITE(1, '( *(g0, ",") )') RHS_lapack(i)

                END DO

                CLOSE(1)
                
        END PROGRAM

        SUBROUTINE GET_ELEM_MAT(x, y, dx, dy, A, A_elem)
                REAL*8, DIMENSION(3), INTENT(IN) :: x, y, dx, dy
                REAL*8, INTENT(IN) :: A

                REAL*8, DIMENSION(3,3), INTENT(OUT) :: A_elem

                REAL*8 :: temp
                INTEGER :: i, j, k
                

                DO i=1,3
                DO j=1,3

                        A_elem(i, j) = (dx(i) * dx(j) + dy(i) * dy(j)) /(4 * A)

                        temp = 0

                        DO k=1,3

                                temp = temp + x(k) + 2 * y(k)

                        END DO

                        A_elem(i,j) = A_elem(i,j) * (1 + A * temp / 30.0)

                END DO
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

        SUBROUTINE CAM_DSOLVE(BWIDTH, LHS, RHS, NUM_EQN)
                INTEGER, INTENT(IN) :: NUM_EQN, BWIDTH
                REAL*8, DIMENSION(3*BWIDTH+1, NUM_EQN), INTENT(IN) :: LHS
                REAL*8, DIMENSION(NUM_EQN), INTENT(INOUT) :: RHS

                INTEGER, DIMENSION(NUM_EQN) :: IPIV
                INTEGER :: INFO

                CALL DGBTRF(NUM_EQN, NUM_EQN, BWIDTH, BWIDTH, LHS, 3 * BWIDTH + 1, IPIV, INFO)

                IF (INFO .NE. 0) THEN
                        PRINT *, "Error with factorization"
                ELSE
                        CALL DGBTRS('N', NUM_EQN, BWIDTH, BWIDTH, 1, LHS, 3 * BWIDTH + 1, IPIV, RHS, NUM_EQN, INFO)
                
                        IF (INFO .NE. 0) THEN
                                PRINT *, "Error with solve"
                        END IF
                END IF 

        END SUBROUTINE
