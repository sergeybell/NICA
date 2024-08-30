      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2),yy1(10000),yy2(10000),t(10000),
     *ae(10000),di(10000),ry(10000),phi(10000),ww(10000),
     *ayy(100),fi(100),ratio(100),velocity(10000)
      external func
      OPEN (UNIT=12,file='test.dat',STATUS='unknown')
      N=2
      pi=3.141592653589793238d0
      omega = 1.0
      m=7000003
      IX=677010339
      LL=5

      Np=10
      Nt=4
      NJ=50
    	P1=0
      Vmax = 1.0d0
      Xmax = 1.0d0

c------different particles w different initial parameters-----
      DO 3 I=1,Np
      P1=P1+1
      call random_number(rand)
      call random_number(rand1)
      ay= Xmax * (2*rand -1)
      by= Vmax * (rand1)

      write(*,*) ay,by
      T1 = 0
c------Turn-by-turn-----
      DO 2 II=1,Nt
      T1 = T1+1

      J1=0
c-----Integration during the cicle
      DO 1 J=1,NJ
      J1=J1+1

      eps=1.D-05
      kr=100
      kmax=1000000
      b=1.d0/NJ
      xa = (2.0*pi) * (J1-1)*b + (2.0 * pi) * (T1-1)
      xb = (2.0*pi) * J1 * b + (2.0 * pi) * (T1-1)

      y(1)=ay
      y(2)=by
      call runkut(xa,xb,y,2,func,eps,kr,kmax)
      ay=y(1)
      by=y(2)
      yy1(J1)=y(1)
      yy2(J1)=y(2)
      write (12,*)  xa,yy1(J1),yy2(J1),J1,T1,P1
    1 continue
c      write (12,*) xa,phi(T1),velocity(T1),J1,T1,P1
    2 continue
      phi(P1)=yy1(J1)
      velocity(P1)=yy2(J1)
c      write (12,*) xa,phi(P1),velocity(P1),J1,T1,P1
    3 continue
      stop
      end

      subroutine func(x,y,f)
      implicit double precision (a-h,o-z)
      dimension y(2),f(2)
c- equation: y" + c*y'+ w**2*y = B * cos(w1 * x)

      f(1)=y(2)
      f(2)= - y(1) - 0.1 * y(2)
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
