        PROGRAM hw7
                INTEGER, PARAMETER :: NUM_ELEM=917 , NUM_NODE=502, NUM_NODE_IN=505, NUM_DND=33, NUM_MEAS=12
                
                INTEGER, DIMENSION(5, NUM_ELEM) :: elem_input
                INTEGER, DIMENSION(3, NUM_ELEM) :: elem
                REAL*8, DIMENSION(3, NUM_NODE) :: node_input
                REAL*8, DIMENSION(2, NUM_NODE) :: node
                REAL*8, DIMENSION(2, NUM_DND) :: dnd_input
                INTEGER, DIMENSION(NUM_DND) :: dnd
                REAL*8, DIMENSION(2, NUM_MEAS) :: meas
                REAL*8, DIMENSION(NUM_NODE) :: truth_u

                REAL*8, DIMENSION(NUM_MEAS, NUM_NODE) :: S
                REAL*8, DIMENSION(NUM_MEAS) :: d
                REAL*8, DIMENSION(3, NUM_ELEM) :: x, y, dx, dy
                REAL*8, DIMENSION(NUM_ELEM) :: A
                INTEGER :: i, j, k, l

                OPEN(1, file="hw44.ele")
                READ(1, *) elem_input
                CLOSE(1)

                elem = elem_input(2:4,:)

                OPEN(1, file="hw44.nod")
                READ(1, *) node_input
                CLOSE(1)

                node = node_input(2:3,:)
                
                OPEN(1, file="hw44.dnd")
                READ(1, *) dnd_input
                CLOSE(1)

                OPEN(1, file="sample_points.dat")
                READ(1, *) meas
                CLOSE(1)

                OPEN(1, file="truth.dat")
                READ(1, *) truth_u
                CLOSE(1)

                DO l=1,NUM_ELEM
                DO i=1,3

                        x(i, l) = node(1, elem(i,l))
                        y(i, l) = node(2, elem(i,l))

                        CALL GET_DELTAS(x(:,l), y(:,l), dx(:,l), dy(:,l))
                        CALL GET_AREA(x(:,l), dy(:,l), A(:,l))

                END DO
                END DO

                
        END PROGRAM

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

        END SUBROUTINE
