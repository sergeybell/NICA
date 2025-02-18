      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2),yy1(10000),yy2(10000),t(10000),
     *ae(10000),di(10000),ry(10000),phi(10000),ww(10000),
     *ayy(100),fi(100),ratio(100)
      external func
      OPEN (UNIT=10,file='phi_w.dat',STATUS='unknown')
      OPEN (UNIT=11,file='ratio.dat',STATUS='unknown')
      OPEN (UNIT=12,file='test2.dat',STATUS='unknown')
      N=2
      pi=3.141592653589793238d0
      omega = 1.0
      m=7000003
      IX=677010339
      LL=5

c      IN=5
c    	I1=0
c      do 2 II=1,IN+1
c      ay=1.0d0*(-1.0+2.0*(II-1)/IN)
c      I1=I1+1
      ay = 0.99
      by = 0.0d0
      Nturn = 2
      NJ=50
      xa = 0.0
      xb = 0.0
      J1=0
c-----turn by turn------------------
      DO 3 I1=1,Nturn
      xa=xa+2*pi
c-----solution of equations---------
      DO 1 J=1,NJ
      J1=J1+1

      eps=1.D-09
      kr=100
      kmax=1000000
      b = (1.d0/NJ)
      xa = (2.0*pi) * (J1-1) * b
      xb = (2.0*pi) * J1 * b
      x = xa + (xb-xa)/2

      y(1)=ay
      y(2)=by
      call runkut(xa,xb,y,2,func,eps,kr,kmax)
      ay=y(1)
      by=y(2)
      yy1(J1)=y(1)
      yy2(J1)=y(2)

      T1 = 0.5
      T2 = 1.0
      Vbb = 1.0
      if(-T2/2-T1.le.y(1) .and. y(1).le.-T2/2) V = (+1.0)*Vbb
      if(-T2/2.lt.y(1) .and. y(1).lt.T2/2)     V = 0.0
      if( T2/2.le.y(1) .and. y(1).le.T1+T2/2)  V = (-1.0)*Vbb
      V = +cos(3.14159*y(1)+3.14159/2)
      write (12,*)  x,yy1(J1),yy2(J1), V
    1 continue
c      write (12,*)  x,yy1(J1),yy2(J1), V
    3 continue
c    2 continue
      stop
      end

      subroutine func(x,y,f)
      implicit double precision (a-h,o-z)
      dimension y(2),f(2)
c- equation: y" + c*y'+ w**2*y = B * cos(w1 * x)

      T1 = 0.5
      T2 = 1.0
      Vbb = 1.0
      if(-T2/2-T1.le.y(1).and.y(1).le.-T2/2) V = (+1.0)*Vbb
      if(-T2/2.lt.y(1).and.y(1).lt.T2/2)     V = 0.0
      if( T2/2.le.y(1).and.y(1).le.T1+T2/2)  V = (-1.0)*Vbb

      f(1) = y(2)
c      f(2) = -y(1)
      f(2) = +cos(3.14159*y(1)+3.14159/2)
c      f(2) = V
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
