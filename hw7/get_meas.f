        PROGRAM hw7
                INTEGER, PARAMETER :: NUM_ELEM=917 , NUM_NODE=502, NUM_NODE_IN=505, NUM_DND=33, NUM_MEAS=12
                
                INTEGER, DIMENSION(5, NUM_ELEM) :: elem_input
                INTEGER, DIMENSION(3, NUM_ELEM) :: elem
                REAL*8, DIMENSION(3, NUM_NODE) :: node_input
                REAL*8, DIMENSION(2, NUM_NODE) :: node
                REAL*8, DIMENSION(2, NUM_DND) :: dnd_input
                INTEGER, DIMENSION(NUM_DND) :: dnd
                REAL*8, DIMENSION(2, NUM_MEAS) :: meas
                REAL*8, DIMENSION(NUM_MEAS) :: min_meas_dist
                INTEGER, DIMENSION(NUM_MEAS) :: min_meas_elem
                REAL*8, DIMENSION(NUM_NODE) :: truth_u

                REAL*8, DIMENSION(NUM_MEAS, NUM_NODE) :: S
                REAL*8, DIMENSION(NUM_MEAS) :: d
                REAL*8, DIMENSION(3, NUM_ELEM) :: x, y, dx, dy, x_temp, y_temp
                REAL*8, DIMENSION(NUM_ELEM) :: A
                REAL*8, DIMENSION(2) :: delta
                REAL*8, DIMENSION(3) :: local
                REAL*8 :: dist
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

                OPEN(1, file="output/truth.dat")
                READ(1, *) truth_u
                CLOSE(1)

                DO l=1,NUM_ELEM
                DO i=1,3

                        x(i, l) = node(1, elem(i,l))
                        y(i, l) = node(2, elem(i,l))

                        CALL GET_DELTAS(x(:,l), y(:,l), dx(:,l), dy(:,l))
                        CALL GET_AREA(x(:,l), dy(:,l), A(l))

                END DO
                END DO

                min_meas_dist = 100
                min_meas_elem = -1

                DO i=1,NUM_MEAS
                DO l=1,NUM_ELEM

                dist = 0
                DO j=1,3
                        delta = meas(:, i) - node(:, elem(j, l))
                        dist = dist + SQRT(delta(1)**2 + delta(2)**2)
                END DO

                IF (dist < min_meas_dist(i)) THEN
                        min_meas_dist(i) = dist
                        min_meas_elem(i) = l
                END IF

                END DO
                END DO

                OPEN(1, file="output/meas_elems.dat")
                WRITE(1, '( *(g0, ",") )') min_meas_elem
                CLOSE(1)

                S = 0
                d = 0

                OPEN(1, file="output/S.dat")
                OPEN(2, file="output/d.dat")

                DO i=1,NUM_MEAS

                CALL GET_AREA_COORDS(meas(:, i), x(:, min_meas_elem(i)), y(:, min_meas_elem(i)), A(min_meas_elem(i)), local)
                
                DO j=1,3

                        S(i, elem(j, min_meas_elem(i))) = local(j)
                        d(i) = d(i) + local(j) * truth_u(elem(j, min_meas_elem(i)))

                END DO

                WRITE(2, '( *(g0, ",") )') d(i)

                END DO

                DO i=1,NUM_NODE
                        WRITE(1, '( *(g0, ",") )') S(:, i)
                END DO

                CLOSE(1)
                CLOSE(2)

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

                A = A / 2.0

        END SUBROUTINE

        SUBROUTINE GET_AREA_XY(x, y, A)
                REAL*8, DIMENSION(3), INTENT(IN) :: x, y
                REAL*8, INTENT(OUT) :: A

                A = ABS(x(2)*y(3) - x(3)*y(2) - x(1)*y(3) + x(3)*y(1) + x(1)*y(2) - x(2)*y(1)) / 2.0

        END SUBROUTINE

        SUBROUTINE GET_AREA_COORDS(meas_pos, x, y, A, L)
                REAL*8, DIMENSION(2), INTENT(IN) :: meas_pos
                REAL*8, DIMENSION(3), INTENT(IN) :: x, y
                REAL*8, INTENT(IN) :: A
                REAL*8, DIMENSION(3), INTENT(OUT) :: L
                
                REAL*8, DIMENSION(3) :: temp_x, temp_y
                REAL*8 :: A_prime
                INTEGER :: i

                DO i=1,3
                        temp_x = x
                        temp_y = y
                        temp_x(i) = meas_pos(1)
                        temp_y(i) = meas_pos(2)

                        CALL GET_AREA_XY(temp_x, temp_y, A_prime)

                        L(i) = A_prime / A
                END DO

        END SUBROUTINE
