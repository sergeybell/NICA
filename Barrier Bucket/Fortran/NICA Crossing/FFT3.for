      program testrun

      INTEGER::M
      REAL*8::free,ftcoef(1000,2),I(1000),X(1000),F(1000),
     *AI(1000)
      OPEN (UNIT=2,file='harmonic.dat',STATUS='unknown')

      M=32
      NN=5
      DO J=0,M-1,1
      OPEN (UNIT=1,file='My Distribution.dat',STATUS='unknown')
      READ(1,*) I(J+1),X(J+1),F(J+1)
      AI(J) = 1.0d0
c      WRITE (*,*) I(J),X(J),F(J)
      END DO
      CLOSE(1)

c      WRITE(*,*) ISHFT(32,-1)

      CALL LFFT(1,NN,M,F,AI)

c      CALL FFT2(F,AI,M,NN)

      DO J=0,M-1
      WRITE(2,*) I(J+1),F(J+1),AI(J+1)
      END DO

      stop
      end

      SUBROUTINE FFT(DATA,NN,ISIGN)                                     FFT00010
C *** This is the Danielson and Lanczos implementation of the fast      FFT00020
C *** Fourier transform as described in Numerical Recipes, Press et     FFT00030
C *** al in section 12.2.  It has been tested by comparing with         FFT00040
C *** THE ORIGINAL COOLEY-TUKEY TRANSFORM, which is a fortran 4         FFT00050
C *** implementation of the same code.                                  FFT00060
C ***  TRANSFORM(K)=SUM(DATA(J)*EXP(ISIGN*                              FFT00070
C ***  2PI*SQRT(-1)*(J-1)*(K-1)/NN)). SUMMED OVER ALL J                 FFT00080
C *** AND K FROM 1 TO NN. DATA IS IN A ONE-DIMENSIONAL                  FFT00090
C *** COMPLEX ARRAY (I.E.,THE REAL AND IMAGINARY                        FFT00100
C *** PARTS ARE ADJACENT IN STORAGE ,SUCH AS FORTRAN IV                 FFT00110
C *** PLACES THEM) WHOSE LENGTH NN=2**K, K.GE.0 (IF NECESSARY           FFT00120
C *** APPEND ZEROES TO THE DATA). ISIGN IS +1 OR -1. IF A -1            FFT00130
C *** TRANSFORM IS FOLLOWED BY A +1 ONE (OR A +1 BY A -1) THE           FFT00140
C *** ORIGINAL DATA REAPPEAR, MULTIPLIED BY NN. TRANSFORM               FFT00150
C *** VALUES ARE RETURNED IN ARRAY DATA, REPLACING THE INPUT.           FFT00160
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(2*NN)
      N=2*NN
      J=1
      DO I=1,N,2
      IF(J.GT.I)THEN
      TEMPR=DATA(J)
      TEMPI=DATA(J+1)
      DATA(J)=DATA(I)
      DATA(J+1)=DATA(I+1)
      DATA(I)=TEMPR
      DATA(I+1)=TEMPI
      END IF
      M=N/2
    1 IF((M.GE.2).AND.(J.GT.M))THEN
      J=J-M
      M=M/2
      GOTO 1
      END IF
      J=J+M
      END DO
C *** Here begins the Danielson-Lanczos section (outer loop executed
C *** Log2 (NN) times
      MMAX=2
    2 IF(N.GT.MMAX)THEN
      ISTEP=2*MMAX
      THETA=6.28318530717959D0/(ISIGN*MMAX)
      WPR=-2*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1
      WI=0
      DO M=1,MMAX,2
      DO I=M,N,ISTEP
      J=I+MMAX
      TEMPR=WR*DATA(J)-WI*DATA(J+1)
      TEMPI=WR*DATA(J+1)+WI*DATA(J)
      DATA(J)=DATA(I)-TEMPR
      DATA(J+1)=DATA(I+1)-TEMPI
      DATA(I)=DATA(I)+TEMP
      DATA(I+1)=DATA(I+1)+TEMPI
      END DO
      WTEMP=WR
      WR=WR*WPR-WI*WPI+WR
      WI=WI*WPR+WTEMP*WPI+WI
      END DO
      MMAX=ISTEP
      GO TO 2
      END IF
      RETURN
      END

      SUBROUTINE FFT2(AR,AI,N,M)
!
! An example of the fast Fourier transform subroutine with N = 2**M.
! AR and AI are the real and imaginary part of data in the input and
! corresponding Fourier coefficients in the output.
! Copyright (c) Tao Pang 1997.
!
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: N,M
      INTEGER :: N1,N2,I,J,K,L,L1,L2
      REAL(8) :: PI,A1,A2,Q,U,V
      REAL(8), INTENT (INOUT), DIMENSION (N) :: AR,AI
      !
      PI = 4.D0*ATAN(1.D0)
      N2 = N/2
      !
      N1 = 2**M
      IF(N1.NE.N) STOP 'Indices do not match'
      !
      ! Rearrange the data to the bit reversed order
      !
      L = 1
      DO K = 1, N-1
      IF (K.LT.L) THEN
      A1    = AR(L)
      A2    = AI(L)
      AR(L) = AR(K)
      AR(K) = A1
      AI(L) = AI(K)
      AI(K) = A2
      END IF
      J   = N2
      DO WHILE (J.LT.L)
      L = L-J
      J = J/2
      END DO
      L = L+J
      END DO
      !
      ! Perform additions at all levels with reordered data
      !
      L2 = 1
      DO L = 1, M
      Q  =  0.D0
      L1 =  L2
      L2 =  2*L1
      DO K = 1, L1
      U   =  DCOS(Q)
      V   = -DSIN(Q)
      Q   =  Q + PI/L1
      DO J = K, N, L2
      I     =  J + L1
      A1    =  AR(I)*U-AI(I)*V
      A2    =  AR(I)*V+AI(I)*U
      AR(I) =  AR(J)-A1
      AR(J) =  AR(J)+A1
      AI(I) =  AI(J)-A2
      AI(J) =  AI(J)+A2
      END DO
      END DO
      END DO
      END SUBROUTINE FFT2

      SUBROUTINE LFFT(DIR,M,N,X,Y)
      INTEGER, INTENT (IN) :: dir,M,N
      INTEGER :: i,i1,j,k,i2,l,l1,l2
      REAL(16) :: c1,c2,tx,ty,t1,t2,u1,u2,z
      REAL(8), INTENT (INOUT), DIMENSION (N) :: X,Y

      i2 = ISHFT(n,-1)
      j=0
      DO i=0,n-2,1
      IF (i.lt.j) THEN
      tx = x(i+1)
      ty = y(i+1)
      x(i+1) = x(j+1)
      y(i+1) = y(j+1)
      x(j+1) = tx
      y(j+1) = ty
      END IF
      k = i2
      DO WHILE (k.le.j)
      j = j - k
      k = ISHFT(k,-1)
      END DO
      j=j+k
      END DO

      c1 = -1.0
      c2 = 0.0
      l2 = 1
      DO l=0,M-1,1
      l1 = l2
      l2 = ISHFT(l2,1)
      u1 = 1.0
      u2 = 0.0
      DO j=0,l1-1,1
      DO i=j,n-1,l2
      i1 = i + l1
      t1 = u1 * x(i1+1) - u2 * y(i1+1)
      t2 = u1 * y(i1+1) + u2 * x(i1+1)
      x(i1+1) = x(i+1) - t1
      y(i1+1) = y(i+1) - t2
      x(i+1) = x(i+1) + t1
      y(i+1) = y(i+1) + t2
      END DO
      z =  u1 * c1 - u2 * c2
      u2 = u1 * c2 + u2 * c1
      u1 = z
      END DO
      c2 = sqrt((1.0 - c1) / 2.0)
      IF (dir .eq. 1) THEN
      c2 = -c2
      END IF
      c1 = sqrt((1.0 + c1) / 2.0)
      END DO

      if (dir .eq. 1) then
      DO i=0,n-1,1
      x(i+1) = x(i+1) / n
      y(i+1) = y(i+1) / n
      END DO
      END IF

      RETURN
      END SUBROUTINE LFFT
