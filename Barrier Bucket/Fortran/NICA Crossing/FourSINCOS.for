      program testrun

      INTEGER :: m
      REAL*8:: n,i,j,k,r(1000),p(1000),x(1000),A(1000),B(1000),twopi

      OPEN (UNIT=1,file='My Distribution.dat',STATUS='unknown')
      OPEN (UNIT=2,file='SINharmonic.dat',STATUS='unknown')

      twopi = 6.28318d0
      m = 5
      n = 32.0d0

      DO j=1,n,1
      OPEN (UNIT=1,file='My Distribution.dat',STATUS='unknown')
      READ(1,*) r(j),p(j),x(j)
      END DO
      CLOSE(1)


      A = 0.0
      B = 0.0
      call ht(x,A,B,n,m)

      DO k=0,n/2
      WRITE(2,*) k,A(k+1),B(k+1)
      END DO

      stop
      end

      subroutine ht(x,A,B,n,m)
      REAL*8, INTENT (IN) :: n
      INTEGER, INTENT (IN) :: m
      INTEGER :: i,j,k
      REAL(16) :: twopi
      REAL(8), INTENT (INOUT), DIMENSION (2**m) :: x
      REAL(8), INTENT (INOUT), DIMENSION (2**(m-1)) :: A,B
      twopi = 6.28318d0

      k=0
      DO k=0,n/2
      DO i=0,n-1
      IF (k.eq.0) then
      A(k+1) = A(k+1) + (1/n) * x(i+1) * cos(twopi*k*i/n)
      END IF
      IF (k.eq.n/2) then
      A(k+1) = A(k+1) + (1/n) * x(i+1) * cos(twopi*k*i/n)
      END IF
      IF ((k.gt.0).and.(k.lt.n/2)) then
      A(k+1) = A(k+1) + (2/n) * x(i+1) * cos(twopi*k*i/n)
      END IF
      B(k+1) = B(k+1) + (2/n) * x(i+1) * sin(twopi*k*i/n)
      END DO
      END DO

      RETURN
      END SUBROUTINE ht
