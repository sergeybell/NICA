      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2)
	COMMON/PARAM/pi,T1,T2,Vacc,Vbb,proton,alpha_0,alpha_1,Etr,
     *Cring,c
      external func
      external Voltage
      OPEN (UNIT=12,file='Voltage.dat',STATUS='unknown')
      OPEN (UNIT=13,file='PhaseSpace.dat',STATUS='unknown')
c      OPEN (UNIT=14,file='Voltage.dat',STATUS='unknown')

c------constants--------
      proton = 938.0d0
      pi = 3.14159
      c = 3*1.D8
      Cring = 503.04
c------integration details--------
      eps=1.D-10
      kr=100.0d0
      kmax=1000000.0d0

c------parameters of Ring and Beam--------
      gamma_tr = 7.086619756
      alpha_0 = 0.0199123143
      alpha_1 = -0.07
      Einj = 6000.0d0
      Etr = 5709.0d0
      Vacc = 300.0 * 1.D-6
      Vbb = 5000 * 1.D-6
      gamma = (Einj/proton)+1
      beta = dsqrt(1-(1/(gamma**2)))
      E = Einj
      T0 = Cring / (c*beta) * 1.D9
      write(*,*) T0

c------voltage--------
      V = 0.0

      T1 = 80.0d0
      T2 = (T0 - 2*T1)/2

c------parameters of turns and particles--------
      tau_max = T2/2+T1-5
      dE_max  = 0.0
      Np = 1
      Nt = 500
      Nturn = 2500

      I1 = 0
c------particles--------
      DO 3 II=1,Np
      I1 = I1+1

c      call random_number(rand1)
c      call random_number(rand2)
c      tau0 = tau_max * (2*rand1-1)
c      dE0 = dE_max * (2*rand2-1)
      delta = 1.0d0/Np
c      tau0 = T2/2 + T1
c      tau0 = tau_max - (I1-1)*delta*T1
      tau0 = 2*tau_max * (I1-1) * delta - (tau_max)
c      tau0 = tau_max * (2*I1-3)
c      tau0 = tau_max
      dE0 = dE_max

      dE_E = dE0 / E
      dp_p = dE_E / (beta**2)
      ay = tau0
      by = dE0
      write(*,*) tau0,dE0,dp_p

      I2 = 0
      xa = 0.0
      xb = 0.0
c------turn by N-turn--------
      DO 2 I=1,Nt
      I2 = I2+1

c------solution of equations--------
      b = Nturn
      xa = (I2-1)*b
      xb = I2 * b

      y(1) = ay
      y(2) = by

      call runkut(xa,xb,y,2,func,eps,kr,kmax,
     *E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p)
      ay=y(1)
      by=y(2)
      yy1=y(1)
      yy2=y(2)

      dE_E = y(2) / E
      dp_p = dE_E / (beta**2)

c      IF(MOD(I2,10).eq.0) write (13,*) x,yyy1,yyy2,dp_p,dE_E,eta,E
      if((yy1.lt.(-T1-(T2/2)-500)).or.(yy1.gt.(T1+(T2/2)+500))) exit
      write (13,*) xa,yy1,yy2,dp_p,dE_E,eta,E,V,gamma
      E = E + (Vacc * Nturn)
    2 continue
    3 continue
      stop
      end

      subroutine func(x,y,f,
     *E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p)
      implicit double precision (a-h,o-z)
	COMMON/PARAM/pi,T1,T2,Vacc,Vbb,proton,alpha_0,alpha_1,Etr,
     *Cring,c
      dimension y(2),f(2)

      gamma = (E/proton) + 1
      beta = dsqrt(1-(1/(gamma**2)))
      eta0 = alpha_0 - (1/(gamma**2))
      ans = alpha_1 + ((3*(beta**2))/(2*(gamma**2)))
      eta1 = ans - eta0/(gamma**2)
      T0 = Cring/(c*beta) * 1.D9
      T2 = (T0 - 2*T1)/2

      dE_E = y(2) / E
      dp_p = dE_E / (beta**2)

      eta = eta0+eta1*dp_p

      if(E.ge.Etr) sign = -1.0
      if(E.lt.Etr) sign = 1.0

      V=0.0
      if(-T1-T2/2.le.y(1).and.y(1).le.-T2/2)
     *V=-1.0*sign*Vbb
      if(-T2/2.lt.y(1).and.y(1).lt.T2/2)
     *V=0.0
      if( T2/2.le.y(1).and.y(1).le.T1+T2/2)
     *V=sign*Vbb

      f(1) = (eta*T0*y(2))/(beta*beta*E)
c      f(2) = Vbb*y(1)/1000
c      f(2) = Vbb*sin((2*pi*y(1))/T0)
      f(2) = V

      return
      end

      subroutine runkut(xa,xb,y,m,rhs,eps,kr,kmax,
     *E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p)
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
      call rhs(x,yin,f,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p)
      do 4 n=1,m
      yh(n)=yin(n)+h*f(n)
    4 continue
      x=xa+(k+k-1)*h
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p)
      do 5 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+h*hk(n)
    5 continue
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p)
      do 6 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+st*hk(n)
    6 continue
      x=xa+k*st
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p)
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
