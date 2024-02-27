      PROGRAM TEST

         INTEGER, DIMENSION(2,2) :: A

         A(1,1) = 1
         A(1,2) = 2
         A(2,1) = 3
         A(2,2) = 4

         PRINT *, A

         A = A * 2

         PRINT *, A

      END PROGRAM
