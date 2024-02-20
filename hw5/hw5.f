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
            
            REAL*8, DIMENSION(4, NUM_MAT) :: mat_properties_input
            REAL*8, DIMENSION(NUM_MAT) :: kappa, c_rho, m

            INTEGER :: i

            OPEN(1, file="epeltr4.dat")
            READ(1, *) elem_list_input
            CLOSE(1)
            
            elem_list = elem_list_input(2:5, :)
            material_list = elem_list_input(6, :)

            OPEN(1, file="npeltr4.dat")
            READ(1, *) node_pos_input
            CLOSE(1)

            node_pos = node_pos_input(2:3, :)

            OPEN(1, file="bpeltr4.dat")
            READ(1, *) bc_list_input
            CLOSE(1)

            bc_list(1, :) = bc_list_input(2, :)
            bc_list(2:5, :) = bc_list_input(4:7, :)

            OPEN(1, file="bpelt4.dat")
            READ(1, *) heat_rate_list_input
            CLOSE(1)

            heat_rate_list = heat_rate_list_input(3,:)

            OPEN(1, file="material_properties.dat")
            READ(1, *) mat_properties_input
            CLOSE(1)

      END PROGRAM
