        PROGRAM TEST

                INTEGER, PARAMETER :: N = 4, N_big = 6
                REAL*8, DIMENSION(N_big**2) :: U, U_true
                REAL*8 :: err, val
                INTEGER :: i

                err = 5.0
                U = 0
                U_true = 1
                DO i=1,N_big**2
                        U_true(i) = i
                END DO

                CALL CALC_ERROR(U, U_true, N_big**2, err)

                DO i=1,N**2
                        PRINT*,i
                CALL INTERP_TRUE(i, U_true, N, N_big, val)
                        PRINT*, val
                END DO

        END PROGRAM

        SUBROUTINE CALC_ERROR(U, U_TRUE, N, err)
                INTEGER, INTENT(IN) :: N
                REAL*8, DIMENSION(N), INTENT(IN) :: U, U_TRUE
                REAL*8, INTENT(OUT) :: err

                err = MAXVAL(ABS(U-U_TRUE))
        END SUBROUTINE

        SUBROUTINE INTERP_TRUE(i, U_true, N, N_true, val)
                INTEGER, INTENT(IN) :: i, N, N_true
                REAL*8, DIMENSION(N), INTENT(IN) :: U_true
                REAL*8, INTENT(OUT) :: val

                INTEGER :: group, item, i_cpy, a, group_flag, item_flag
                INTEGER :: idx
                REAL*8 :: ratio, r_i, r_g, alpha, beta, gamma, delta

                ratio = DBLE(N_true-1)/DBLE(N-1)

                group = (i-1)/N + 1
                item = MOD(i-1, N)+1

                PRINT *, group, item

                r_i = ratio * DBLE(item-1)
                r_i = r_i - DBLE(INT(r_i))
                r_g = ratio * DBLE(group-1)
                r_g = r_g - DBLE(INT(r_g))

                group = FLOOR(DBLE(group-1)*ratio)+1
                item = FLOOR(DBLE(item-1)*ratio)+1
                PRINT *, group, item


                idx = N_true * (group - 1) + item
                alpha = (1 - r_i) * (1 - r_g) * U_true(idx)

                idx = N_true * (group - 1) + item + 1
                beta = r_i * (1 - r_g) * U_true(idx)

                idx = N_true * group + item
                gamma = (1 - r_i) * r_g * U_true(idx)

                idx = N_true * group + item + 1
                delta = r_i * r_g * U_true(idx)

                val = alpha + beta + gamma + delta

        END SUBROUTINE
