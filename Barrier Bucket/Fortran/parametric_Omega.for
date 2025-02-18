      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2),yy1(10000),yy2(10000),t(10000),
     *ae(10000),di(10000),ry(10000),phi(10000),ww(10000),
     *ayy(100),fi(100),ratio(100)
      COMMON /NACH/ IX
	COMMON /PARAM/ fis,Omega,am,ff
      external fsin
      OPEN (UNIT=10,file='phi_w.dat',STATUS='unknown')
      OPEN (UNIT=11,file='ratio.dat',STATUS='unknown')
      OPEN (UNIT=12,file='test.dat',STATUS='unknown')
      N=2
      pi=3.141592653589793238d0
      m=7000003
      IX=677010339
      LL=5
c-------Change of phase advance due to synchronous phase
	fis=60.d00/180.d00*pi
	Omega=dsqrt(cos(fis))
c-------Change of effective amplitude due to synchronous phase oscillation
	am=80.d00/180.d00*pi
	ff=1./dsqrt(DBESJ0(am))
c-------------------------------------------------------------------
	IN=36
	INS=19
	I1=0

      DO Q=1,100
      call random_number(p)
      write(*,*) p
      END DO

      Xmax=30./180.d00*pi
	I1=0
      do 3 II=1,IN+1
      I1=I1+1
      ay=Xmax*(-1.0+2.0*(II-1)/IN)
		if (II.eq.INS) ay=0.01
	ayy(II)=ay*180.d00/pi

      by=0.0

      J1=0
      DO 1 I=1,N
      NJ=100
      DO 2 J=1,NJ
      J1=J1+1
      di(J1)=J1
      eps=1.D-10
      kr=100
      kmax=1000000
      b=1.d0/NJ
      xa=Omega*(J1-1)*b*2.0D0*3.141592653589793238d0
      xb=Omega*J1*b*2.0D0*3.141592653589793238d0

      y(1)=ay
      y(2)=by
      call runkut(xa,xb,y,2,fsin,eps,kr,kmax)
      ay=y(1)
      by=y(2)
      yy1(J1)=y(1)
      yy2(J1)=y(2)
    2 continue
    1 continue
      phi(I1)=yy1(J1)
	    ww(I1)=yy2(J1)
    3 continue

      DO 13 I=1,IN+1
     	fi(I)=180.d00/pi*atan(ww(I)/phi(I)/Omega)
	write (12,6) fi(I),ayy(I)
   13 continue
      DO 14 I=1,IN+1
      ratio(I)=(360.d00*N+fi(I))/(360.d00*N)
	if(abs(ayy(I)).gt.110.d00) ratio(I)=(180.d00+fi(I))/(360.d00*N)
   14 continue

            write (11,6) (ayy(I),ratio(I),I=1,IN+1)

            write (10,6) (phi(I),ww(I),I=1,IN+1)
    6       FORMAT (2E12.4)
   15 continue
      stop
      end

      FUNCTION F552RD(LL)
      COMMON/NACH/ IX
      IY=IX*65539
      IF (IY) 1,2,2
    1 IY=IY+2147483647+1
    2 YEL=IY
      B=YEL*.4656613E-9
      F552RD=2*B-1.
      IX=IY
      RETURN
      END

c     m -integer, 0<m<33554432
      double precision function rndm(m)
      implicit double precision (a-h,o-z)
      rndm=m*2.98023223876953d-8
c 1.d0/33554432=2.98023223876953d-8
      m=(m+2)*61
      n=m-1073741824
      if(n.gt.0) m=n
      n=m-536870912
      if(n.gt.0) m=n
      n=m-268435456
      if(n.gt.0) m=n
      n=m-134217728
      if(n.gt.0) m=n
      n=m-67108864
      if(n.gt.0) m=n
      n=m-33554432
      if(n.gt.0) m=n
      return
      end

      subroutine fsin(x,y,f)
      implicit double precision (a-h,o-z)
	COMMON /PARAM/ fis,Omega,am,ff
      dimension y(2),f(2)
c- equation: y"-2d*y'+�**2*(1+h cos((2�+anu)*x))=0, �=1 for normalization
      pi=3.1415962
	    nu=10

	    fisrf=fis+am*dcos(nu*x)
      h1=0.0d0
    	h2=0.0
    	anu=-0.0
    	delta=0.0
      f(1)=-y(2)
      f(2)=-2.0d00*delta*y(2)+ff**2*(dsin(y(1)+fisrf)-dsin(fisrf))*
     *(1.d00+h1*dcos((3.0*Omega+anu)*x))
      return
      end


      subroutine runkut(xa,xb,y,m,rhs,eps,kr,kmax)
      implicit double precision (a-h,o-z)
c  m.le.20
      dimension y(m),yin(20),un(20),vn(20),f(20),hk(20),yh(20),yc(20)
      logical lc
      eps=dabs(eps)
      kr=kr/2
      if(kr.lt.1) kr=1
      st=(xb-xa)/kr
      do 1 n=1,m
      yc(n)=y(n)
    1 continue
    2 x=xa
      h=(xb-xa)/(kr+kr)
      d=st/6.d0
      do 3 n=1,m
      yin(n)=y(n)
      un(n)=y(n)
      vn(n)=0.d0
    3 continue
      do 8 k=1,kr
      call rhs(x,yin,f)
      do 4 n=1,m
      yh(n)=yin(n)+h*f(n)
    4 continue
      x=xa+(k+k-1)*h
      call rhs(x,yh,hk)
      do 5 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+h*hk(n)
    5 continue
      call rhs(x,yh,hk)
      do 6 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+st*hk(n)
    6 continue
      x=xa+k*st
      call rhs(x,yh,hk)
      do 7 n=1,m
      f(n)=(f(n)+hk(n))*d
      un1=un(n)+f(n)
      vn(n)=vn(n)-((un1-un(n))-f(n))
      yin(n)=un1+vn(n)
      un(n)=un1
    7 continue
    8 continue
      kr=kr+kr
      lc=.true.
      do 9 n=1,m
      lc=lc.and.dabs(yin(n)-yc(n)).lt.eps
    9 continue
      if(lc.or.kr.gt.kmax) go to 11
      do 10 n=1,m
      yc(n)=yin(n)
   10 continue
      st=h
      go to 2
   11 do 12 n=1,m
      y(n)=yin(n)
   12 continue
      return
      end
