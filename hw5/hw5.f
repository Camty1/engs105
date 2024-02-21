      PROGRAM HW5
            INTEGER, PARAMETER :: NUM_ELEM=1089, NUM_NODE=969, NUM_BC=64, NUM_MAT=8

            REAL*8, DIMENSION(6, NUM_ELEM) :: elem_list_input
            REAL*8, DIMENSION(4, NUM_ELEM) :: elem_list
            REAL*8, DIMENSION(NUM_ELEM) :: material_list

            REAL*8, DIMENSION(3, NUM_NODE) :: node_pos_input
            REAL*8, DIMENSION(2, NUM_NODE) :: node_pos

            REAL*8, DIMENSION(7, NUM_BC) :: bc_list_input
            REAL*8, DIMENSION(5, NUM_BC) :: bc_list

            REAL*8, DIMENSION(3, NUM_NODE) :: heat_rate_list_input
            REAL*8, DIMENSION(NUM_NODE) :: heat_rate_list
            
            REAL*8, DIMENSION(4, NUM_MAT - 2) :: mat_properties_input
            REAL*8, DIMENSION(3, NUM_MAT) :: mat_properties

            REAL*8, DIMENSION(:,:), ALLOCATABLE :: LHS
            REAL*8, DIMENSION(NUM_NODE) :: RHS

            INTEGER :: i, j, k, bandwidth, width

            OPEN(1, file="epeltr4.dat")
            READ(1, *) elem_list_input
            CLOSE(1)
            
            OPEN(1, file="npeltr4.dat")
            READ(1, *) node_pos_input
            CLOSE(1)

            OPEN(1, file="bpeltr4.dat")
            READ(1, *) bc_list_input
            CLOSE(1)

            OPEN(1, file="bpelt4.dat")
            READ(1, *) heat_rate_list_input
            CLOSE(1)

            OPEN(1, file="material_properties.dat")
            READ(1, *) mat_properties_input
            CLOSE(1)

            elem_list = elem_list_input(2:5, :)
            material_list = elem_list_input(6, :)

            node_pos = node_pos_input(2:3, :)

            bc_list(1, :) = bc_list_input(2, :)
            bc_list(2:5, :) = bc_list_input(4:7, :)

            heat_rate_list = heat_rate_list_input(3,:)

            mat_properties(:, 3:8) = mat_properties_input(2:4, :)

            bandwidth = 0

            DO i=1,NUM_ELEM
            DO j=1,4
            DO k=j,4
                  width = ABS(elem_list(j, i) - elem_list(k, i)) - 1

                  IF (width > bandwidth) THEN
                        bandwidth = width

                  END IF

            END DO
            END DO
            END DO
            
            ALLOCATE(LHS(NUM_NODE, 2*bandwith+1))
      END PROGRAM

      SUBROUTINE ELEMENT_MAT_RECT()

      END SUBROUTINE

      SUBROUTINE ELEMENT_MAT_TRI(dx, dy, kappa, m, A, order, num_quad, quad_points, quad_weights, elem_mat)
            REAL*8, DIMENSION(3), INTENT(IN) :: dx, dy
            REAL*8, INTENT(IN) :: kappa, m, A
            INTEGER, INTENT(IN) :: order, num_quad
            REAL*8, DIMENSION(4, num_quad), INTENT(IN) :: quad_points
            REAL*8, DIMENSION(num_quad), INTENT(IN) :: quad_weights

            REAL*8, DIMENSION(3, 3), INTENT(OUT) :: elem_mat
            
            REAL*8, DIMENSION(6,3) :: multiplicity_mat
            REAL*8 :: dphi, coeff
            INTEGER :: i, j, k, l

            elem_mat = 0
            
            IF (order == 1) THEN
                  
            DO i=1,3
            DO j=1,3
                  coeff = 0
                  IF (j == 1) THEN
                        dphi = dy(1) - dx(1)

                  ELSE IF (j == 2) THEN
                        dphi = dy(2) - dx(2)

                  ELSE 
                        dphi = -dy(1) - dy(2) + dx(1) + dx(2)
               
                  END IF

                  DO k=1,num_quad
                  multiplicity = quad_points(4, k)
                  
                  ! Don't need a function for this, can just layout the
               ! matrix with ((1,2,3), (3,1,2), (2,3,1), ...) where you
               ! get the multiplicity already loaded
                  CALL GET_MULTIPLICITY(multiplicity, multiplicity_mat)

                  DO l=1,multiplicity

                  END DO
                  END DO
            END DO
            END DO

            ELSE
                  WRITE *, "Higher order triangular elements have not been implemented."

            END IF

      END SUBROUTINE
