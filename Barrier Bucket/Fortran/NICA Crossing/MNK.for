      program tesrun
c***** Integrating motion in separatrix****
      INTEGER :: J
      REAL*8 :: n_bin(1000),bin_count(1000),Distr(1000),phi(1000),
     *x(5),f(5),A(5,6),c(5)
      real*8 :: M,e

      N = 32
      DO J=1,N,1
      OPEN (UNIT=1,file='results/Distr/8InDistribution.dat')
      READ(1,*) n_bin(J),bin_count(J),Distr(J),phi(J)
      END DO
      CLOSE(1)

      x = (/1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)
      f = (/1.0d0,3.0d0,6.0d0,8.0d0,10.0d0/)

      M=4

      A=0.0d0
      call ex(2.0,4,result)
      CALL Gram(32,4,n_bin,bin_count,A,ex)
      CALL Gauss(4,A,c)
      write(*,*) c

      stop
      end

      subroutine ex(a,n,result)
      REAL*8, INTENT (IN) :: a
      INTEGER, INTENT (IN) :: n
      INTEGER :: i
      REAL*8 :: result,e

      e=1.0
      DO i=1,n
      e = e*a
      END DO
      result = e

      END subroutine ex

      subroutine Gram(n,m,x,f,A,rhs)
      INTEGER, INTENT (IN) :: n,m
      INTEGER :: i,j
      REAL*8 :: p,q,r,s,d,result
      REAL*8, INTENT (INOUT), DIMENSION (n) :: x,f
      REAL*8, INTENT (INOUT), DIMENSION (m+1,m+2)::A
      A=0.0d0

      DO j=0,m,1
      s=0.0
      r=0.0
      q=0.0
      DO i=0,n-1,1
      call rhs(x(i+1),j,p)
      s=s+p
      r=r+p*f(i+1)
      call rhs(x(i+1),m,d)
      q=q+p*d
      END DO
      A(1,j+1)=s
      A(j+1,m+1)=q
      A(j+1,m+1+1)=r
      END DO

      DO i=1,m,1
      DO j=0,m-1
      A(i+1,j+1)=A(i,j+2)
      END DO
      END DO

      RETURN
      END subroutine Gram

      subroutine Gauss(m,A,c)
      INTEGER, INTENT (IN) :: m
      INTEGER :: i,j,j1,k,l,k1,n1
      real*8 :: r,s
      REAL*8, INTENT (INOUT), DIMENSION (m+1,m+2)::A
      REAL*8, INTENT (INOUT), DIMENSION (m+1)::c
      c=0.0

      n1=m+1
      DO k=0,m
      k1=k+1
      s=A(k+1,k+1)
      DO j=k1,n1
      A(k+1,j+1)=A(k+1,j+1)/s
      END DO

      DO i=k1,m
      r=A(i+1,k+1)
      DO j=k1,n1
      A(i+1,j+1)=A(i+1,j+1)-A(k+1,j+1)*r
      END DO

      END DO
      END DO

      DO i=m,0,-1
      s=a(i+1,n1+1)
      DO j=i+1,m
      s=s-A(i+1,j+1)*c(j+1)
      ENDDO
      c(i+1)=s
      ENDDO

      RETURN
      END subroutine Gauss
