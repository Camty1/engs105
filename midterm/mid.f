      PROGRAM mid
      REAL*8, PARAMETER :: L = 10, D = 0.5, dx = 0.1, dt = 0.005
      REAL*8, PARAMETER :: x_0 = 5, sigma = 0.1
      INTEGER, PARAMETER :: TIMESTEPS = 100

      REAL*8, DIMENSION(:), ALLOCATABLE :: x, U, U_old, U_new, b
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: A

      REAL*8 :: r
      INTEGER :: N, i
      CHARACTER(LEN = 50) :: filename

      r = D * dt / dx ** 2

      PRINT *, r

      N = CEILING(L / dx) + 1

      ALLOCATE(x(N), U(N), U_old(N), U_new(N), A(N,3), b(N))

      DO i=1,N
            x(i) = dx * (i-1)
      END DO

      ! OUTPUT x file
      OPEN(1, file="output/x.dat")
      WRITE(1, '(*(g0, "," ))') x 
      CLOSE(1)

      CALL IC(x_0, sigma, N, x, U)
      U_old = U
      
      WRITE(filename, '("output/A_", ES7.1, "_", ES7.1, "_", ES7.1, "_", ES7.1, ".dat")') D, sigma, dx, dt
      OPEN(1, file=filename)

      DO i=1,TIMESTEPS
            CALL TIME_STEP_A(r, N, U, U_old, U_new)
            U_old = U
            WRITE(1, '(*(g0, "," ))') U
            U = U_new
      END DO

      WRITE(1, '(*(g0, "," ))') U
      CLOSE(1)

      CALL IC(x_0, sigma, N, x, U)
      U_old = U

      CALL LHS_B(r, N, A)
      
      CALL DSOLVE(1, A, U, N, 1, N, 3)

      WRITE(filename, '("output/B_", ES7.1, "_", ES7.1, "_", ES7.1, "_", ES7.1, ".dat")') D, sigma, dx, dt
      OPEN(1, file=filename)

      DO i=1,TIMESTEPS
            WRITE(1, '(*(g0, "," ))') U
            CALL RHS_B(r, N, U, U_old, b)
            CALL DSOLVE(2, A, b, N, 1, N, 3)
            U_old = U
            U = b
      END DO

      WRITE(1, '(*(g0, "," ))') U
      CLOSE(1)

      END PROGRAM

      SUBROUTINE TIME_STEP_A(r, N, U, U_old, U_new)
      REAL*8, INTENT(IN) :: r
      INTEGER, INTENT(IN) :: N
      REAL*8, DIMENSION(N), INTENT(IN) :: U, U_old
      REAL*8, DIMENSION(N), INTENT(OUT) :: U_new

      INTEGER :: i
      
      U_new = 0.0

      DO i=2,N-1
            U_new(i) = U(i) + r / 2.0 * (3 * (U(i+1) - 2 * U(i) + U(i-1)) - (U_old(i+1) - 2 * U_old(i) + U_old(i-1)))
      END DO

      END SUBROUTINE

      SUBROUTINE LHS_B(r, N, A)
      REAL*8, INTENT(IN) :: r
      INTEGER, INTENT(IN) :: N
      REAL*8, DIMENSION(N,3), INTENT(OUT) :: A

      A = 0.0
      A(1,2) = 1.0
      A(N,2) = 1

      A(2:N-1, 1) = -5 * r / 12.0
      A(2:N-1, 2) = 1 + 10 * r / 12.0
      A(2:N-1, 3) = -5 * r / 12.0

      END SUBROUTINE

      SUBROUTINE RHS_B(r, N, U, U_old, b)
      REAL*8, INTENT(IN) :: r
      INTEGER, INTENT(IN) :: N
      REAL*8, DIMENSION(N), INTENT(IN) :: U, U_old
      REAL*8, DIMENSION(N), INTENT(OUT) :: b

      INTEGER :: i
      
      b = 0.0

      DO i=2,N-1
            b(i) = U(i) + r / 12.0 * (8 * (U(i+1) - 2 * U(i) + U(i-1)) - (U_old(i+1) - 2 * U_old(i) + U_old(i-1)))
      END DO

      END SUBROUTINE

      SUBROUTINE IC(x_0, sigma, N, x, U)
      REAL*8, INTENT(IN) :: x_0, sigma
      INTEGER, INTENT(IN) :: N
      REAL*8, DIMENSION(N), INTENT(IN) :: x
      REAL*8, DIMENSION(N), INTENT(OUT) :: U

      INTEGER :: i

      U = 0.0

      DO i=2,N-1
            U(i) = EXP(-((x(i) - x_0) ** 2 / (2 * sigma ** 2)))
      END DO

      END SUBROUTINE 
