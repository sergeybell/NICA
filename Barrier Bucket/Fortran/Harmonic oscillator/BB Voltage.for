      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2)
	COMMON /PARAM/ pi,T0,T1,T2,Vacc,Vbb
      external func
      external Voltage
      OPEN (UNIT=12,file='Voltage.dat',STATUS='unknown')
      OPEN (UNIT=13,file='PhaseSpace.dat',STATUS='unknown')

c------constants--------
      proton = 938.0d0
      pi = 3.14159
c------integration details--------
      eps=1.D-10
      kr=100.0d0
      kmax=1000000.0d0

c------voltage--------
      V = 0.0

      T0 = 1688.0d0
      T1 = 80.0d0
      T2 = (T0 - 2*T1)/2
c------parameters of Ring and Beam--------
      C = 503.04
      gamma_tr = 7.086619756
      alpha_0 = 0.0199123143
      alpha_1 = -0.07
      Einj = 4000.0d0
      Vacc = 0.0 * 1.D-6
      Vbb = 5000 * 1.D-6
      gamma = (Einj/proton)+1
      E = Einj

c------parameters of turns and particles--------
      tau_max = T2/2 + T1 - 79
      dE_max = 0.0
      Np = 6
      Nt = 1000
      NJ = 1
      Nturn = 500

      I1 = 0
c------particles--------
      DO 3 II=1,Np
      I1 = I1+1

c      tau0 = tau_max * (2*rand1-1)
c      dE0 = dE_max * (2*rand2-1)
      delta = 1.0d0/Np
      tau0 = 2*tau_max * (I1-1) * delta - (tau_max)
c      tau0 = tau_max * (2*I1-3)
      dE0 = dE_max
      ay = tau0
      by = dE0
      write(*,*) tau0,dE0

      I2 = 0
      xa = 0.0
      xb = 0.0
c------turn by N-turn--------
      DO 2 I=1,Nt
      I2 = I2+1
      E = Vacc*(I2) * Nturn + Einj

      gamma = (E/proton) + 1
      beta = dsqrt(1-(1/(gamma**2)))
      eta0 = alpha_0 - (1/(gamma**2))
      ans = alpha_1 + ((3*(beta**2))/(2*(gamma**2)))
      eta1 = ans - eta0/(gamma**2)
c      write(*,*) 5,I1,I2,E,gamma,beta
      J1=0
c------solution of equations--------
      J1 = J1 + 1
      b = Nturn
      xa = (J1-1)*b
      xb = J1 * b
      x = xa + (xb-xa)/2
      eta = eta0

      y(1) = ay
      y(2) = by

      call Voltage(ay,V)
      call runkut(xa,xb,y,2,func,eps,kr,kmax,V,eta,beta,E)
      ay=y(1)
      by=y(2)
      yy1=y(1)
      yy2=y(2)

      dE_E = yy2 / E
      dp_p = dE_E / (beta**2)
c      write (12,*)  x,ay,V,J1,I2,I1, eta, E, beta
c      write (13,*)  x,yy1,yy2,J1,I2,I1,gamma,E
c      write(*,*) I1,I2,J1
c      IF(MOD(I2,1000).eq.0 .and. MOD(J1,10).eq.0) write (13,*) x,yy1,yy2,gamma,eta0,eta1,eta

      yyy1 = yy1
      yyy2 = yy2
      if((yyy1.lt.(-T1-(T2/2))).or.(yyy1.gt.(T1+(T2/2)))) exit
      write (13,*)  x,yyy1,yyy2,dp_p,dE_E,eta,E,gamma
c      IF(MOD(I2,100).eq.0) write (13,*) x,yyy1,yyy2,I2,gamma,E,eta
c      IF(MOD(I2,10).eq.0) write (*,*) x,yyy1,yyy2,I2,gamma,E,eta
    2 continue
      tau = yyy1
      dE  = yyy2
c      write(13, *) x,tau,dE,I1,E,beta,eta
    3 continue
      stop
      end

      subroutine func(x,y,f,V,eta,beta,E)
      implicit double precision (a-h,o-z)
	COMMON /PARAM/ pi,T0,T1,T2,Vacc,Vbb
      dimension y(2),f(2)

      f(1) = (eta*T0)/((beta**2)*E)*y(2)
      f(2) = V
      return
      end

      subroutine Voltage(t,V)
      implicit double precision (a-h,o-z)
	COMMON /PARAM/ pi,T0,T1,T2,Vacc,Vbb

      V = 0.0d0
      if(-T1-T2/2.le.t .and. t.lt.-T2/2) V = (-1.0)*Vbb
      if(-T2/2.le.t .and. t.lt.T2/2)     V = Vacc
      if( T2/2.le.t .and. t.le.T1+T2/2)  V = Vbb
      return
      end

      subroutine runkut(xa,xb,y,m,rhs,eps,kr,kmax,V,eta,beta,E)
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
      call rhs(x,yin,f,V,eta,beta,E)
      do 4 n=1,m
      yh(n)=yin(n)+h*f(n)
    4 continue
      x=xa+(k+k-1)*h
      call rhs(x,yh,hk,V,eta,beta,E)
      do 5 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+h*hk(n)
    5 continue
      call rhs(x,yh,hk,V,eta,beta,E)
      do 6 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+st*hk(n)
    6 continue
      x=xa+k*st
      call rhs(x,yh,hk,V,eta,beta,E)
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
