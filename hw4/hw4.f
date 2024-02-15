        PROGRAM hw4
                

        END PROGRAM

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
