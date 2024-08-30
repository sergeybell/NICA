      program testrun

      INTEGER::i,j
      REAL*8:: m,k,n,r(1000),p(1000),x(1000),A(1000),B(1000),
     *twopi,Distr(1000),tau(1000),num(1000),Np,tau_bin(1000),
     *bin_count(1000),particles,x1,x2

      OPEN (UNIT=1,file='My Distribution.dat',STATUS='unknown')
      OPEN (UNIT=2,file='SINharmonic.dat',STATUS='unknown')
      OPEN (UNIT=3,file='InitPhaseSpace.dat',STATUS='unknown')

      twopi = 6.28318d0
      m = 5.0d0
      n = 32.0d0
      Np = 128.0d0
      k = 32.0d0
      x1 = -465.0d0
      x2 = +465.0d0
      particles = Np

      write(*,*) INT(m,1)
c      DO j=1,n,1
c      OPEN (UNIT=1,file='My Distribution.dat',STATUS='unknown')
c      READ(1,*) r(j),p(j),x(j)
c      END DO
c      CLOSE(1)

      DO j=1,Np,1
      OPEN (UNIT=3,file='results/InitPhaseSpace.dat',STATUS='unknown')
      READ(3,*) num(j),tau(j)
c      WRITE(*,*) num(j),tau(j)
      END DO
      CLOSE(3)

      A = 0.0
      B = 0.0

      write(*,*) Np,k,particles,x1,x2
      call hist(tau,Np,k,particles,x1,x2,Distr,
     *tau_bin,bin_count)

      DO j=1,Np
      WRITE(*,*) tau_bin(j)
      END DO

      DO j=1,k
      WRITE(*,*) j,Distr(j),bin_count(j)
      END DO

      call ht(Distr,A,B,n)

      DO k=0,n/2
      WRITE(2,*) k,A(k+1),B(k+1)
      END DO

      stop
      end

      subroutine ht(x,A,B,n)
      REAL*8, INTENT (IN) :: n
      INTEGER :: i,j,k
      REAL(16) :: twopi
      REAL(8), INTENT (INOUT), DIMENSION (INT(n,1)) :: x
      REAL(8), INTENT (INOUT), DIMENSION (INT(n/2,1)) :: A,B
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

      subroutine hist(tau,Np,k,particles,x1,x2,Distr,
     *tau_bin,bin_count)
      REAL*8, INTENT (IN) :: Np,k,particles,x1,x2
      INTEGER :: i,j,count
      REAL*8 :: h,ik
      REAL(8), INTENT (INOUT), DIMENSION (INT(Np,1)) :: tau
      REAL(8), INTENT (INOUT), DIMENSION (INT(k,1)) :: Distr
      REAL*8:: tau_bin(INT(Np,1))
      REAL*8:: bin_count(INT(k,1))

      h = (x2-x1)/k
      write(*,*) Np,k,particles,x1,x2
      write(*,*) h

      DO i=1,k
      count = 0
      DO j=1,Np
      IF((tau(j).lt.x1).or.(tau(j).gt.x2)) then
      tau_bin(j) = -1
      cycle
      END IF
      IF(((x1 +(I-1)*h).le.tau(j)).and.(tau(j).lt.(x1+(I-1)*h+h)))
     *then
      count = count+1
      tau_bin(J) = I
      END IF
      END DO
      bin_count(I) = count
      Distr(I) = count / particles
      END DO

      RETURN
      END SUBROUTINE hist
