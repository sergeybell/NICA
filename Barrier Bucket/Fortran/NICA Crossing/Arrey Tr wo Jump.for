      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2),turn_p(1000),tau(1000),dE(1000),
     *dp_p(1000),dE_E(1000),eta_p(1000),E_p(1000),V_p(1000),
     *gamma_p(1000),g_tr_p(1000),al_0_p(1000),
     *al_1_p(1000),eta0_p(1000),eta1_p(1000),flag(1000),
     *dphi(1000),beta_p(1000)
	COMMON/PARAM/T1,T2,Vacc,Vbb,
     *T0,dn,n_tr,dn_jump,Etr
      external func
      external Voltage
      OPEN (UNIT=13,file='results/PhaseSpace1.dat',STATUS='unknown')
      OPEN (UNIT=14,file='results/InitPhaseSpace1.dat',STATUS='unknown')
      OPEN (UNIT=15,file='results/PhaseSpaceA1.dat',STATUS='unknown')
c------constants--------
      proton = 938.0d0
      pi = 3.14159
      c = 3*1.D8
      Cring = 503.04
c------integration details--------
      eps=1.D-10
      kr=100.0d0
      kmax=1000000.0d0

c------voltage--------
      V = 0.0
      T0 = 1700.0d0
      T1 = 80.0d0
      T2 = (T0 - 2*T1)/2
      Vacc = 300.0 * 1.D-6
      Vbb = 5000 * 1.D-6

c------parameters of Ring and Beam--------
      alpha_0 = 0.0199123143
      alpha_1 = -0.06872168475
      gamma_tr = 1/sqrt(alpha_0)
      g_tr = gamma_tr

      delta_g_tr = 0.045
      dgdt = 8.5d0
      A1 = 7.6629188031276*1.D-2
      A2 =-6.1171895095566*1.D-1

      Einj = 5600.0d0
      Etr = 5709.0d0
      E = Einj

      gamma = (E/proton)+1
      beta = dsqrt(1-(1/(gamma**2)))

c------parameters of turns and particles--------
      Np = 500.0d0
      Nt = 50.0d0
      Sturn = 20000.0d0
      Nturn = Sturn

      dp_p_max = 0.005
      tau_max = T2/2+T1
      dE_max  = Einj*dp_p_max
      write(*,*) dE_max

      I1 = 0
      DO IX=1,Np
      I1 = I1+1
      call random_number(rand1)
      call random_number(rand2)

      tau0 = tau_max * (2*rand1-1)
      dE0 = dE_max * (2*rand2-1)

      tau(I1) = tau0
      dE(I1) = dE0

      write (14,*) I1,tau(I1),dE(I1)
      END DO

      xa = 0.0
      xb = 0.0
      turn=0.0

      I2 = 0
c------turn by N-turn--------
      DO 2 I=1,Nt
      I2 = I2+1
      g_tr = gamma_tr
      al_1 = alpha_1

      gamma = (E/proton) + 1
      beta = dsqrt(1-(1/(gamma**2)))

      al_0 = 1 / (g_tr**2)
      eta0 = al_0 - (1/(gamma**2))
      ans = al_1 + ((3*(beta**2))/(2*(gamma**2)))
      eta1 = ans - eta0/(gamma**2)

      I1 = 0
c------particles--------
      DO 3 J=1,Np
      I1 = I1+1
      flag(I1)=1
      if((tau(I1).lt.(-T1-(T2/2)-100)).or.
     *(tau(I1).gt.(T1+(T2/2)+100)))
     *flag(I1)=0
      if (flag(I1).eq.0.0) cycle

      turn_p(I1) = turn
      E_p(I1) = E
      gamma_p(I1) = gamma
      beta_p(I1) = beta
      g_tr_p(I1) = g_tr
      al_0_p(I1) = al_0
      al_1_p(I1) = al_1
      eta0_p(I1) = eta0
      eta1_p(I1) = eta1

      y(1) = tau(I1)
      y(2) = dE(I1)

      eta_p(I1) = eta0_p(I1) + eta1_p(I1)*dp_p(I1)
c------solution of equations--------
      b = Nturn
      xa = 0.0
      xb = b

      call runkut(xa,xb,y,2,func,eps,kr,kmax,
     *E_p(I1),eta_p(I1),V_p(I1),gamma_p(I1),beta_p(I1),
     *eta0_p(I1),eta1_p(I1),dE_E(I1),dp_p(I1),turn_p(I1))

      tau(I1) = y(1)
      dE(I1) = y(2)

      dphi(I1) = 2*3.14159*tau(I1)/(T0/2)
      dE_E(I1) = dE(I1) / E_p(I1)
      dp_p(I1) = dE_E(I1) / (beta_p(I1)*beta_p(I1))

      write (15,*) turn_p(I1),tau(I1),dE(I1),dp_p(I1),dE_E(I1),
     *eta_p(I1),E_p(I1),V_p(I1),gamma_p(I1),g_tr_p(I1),
     *al_0_p(I1),al_1_p(I1),eta0_p(I1),eta1_p(I1),flag(I1),
     *dphi(I1),I1,I2

    3 END DO
      turn = turn+Nturn
      E = E + (Vacc * (Nturn))
    2 END DO

      I1 = 0
      DO IX=1,Np
      I1 = I1+1
      write (13,*) turn_p(I1),tau(I1),dE(I1),dp_p(I1),dE_E(I1),
     *eta_p(I1),E_p(I1),V_p(I1),gamma_p(I1),g_tr_p(I1),
     *al_0_p(I1),al_1_p(I1),eta0_p(I1),eta1_p(I1),flag(I1),
     *dphi(I1),I1,I2
      END DO
      stop
      end


      subroutine func(x,y,f,
     *E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,flag)
      implicit double precision (a-h,o-z)
      COMMON/PARAM/T1,T2,Vacc,Vbb,
     *T0,dn,n_tr,dn_jump,Etr
      dimension y(2),f(2)

      sign = 0.0
      if(E.lt.Etr) sign = 1.0
      if(E.ge.Etr) sign = -1.0

      V=0.0
      if(-T1-T2/2.le.y(1).and.y(1).le.-T2/2)
     *V=-sign*Vbb
      if(-T2/2.lt.y(1).and.y(1).lt.T2/2)
     *V=0.0
      if( T2/2.le.y(1).and.y(1).le.T1+T2/2)
     *V=sign*Vbb

      f(1) = (eta*T0*y(2))/(beta*beta*E)
      f(2) = V
      return
      end

      subroutine runkut(xa,xb,y,m,rhs,eps,kr,kmax,
     *E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn)
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
      call rhs(x,yin,f,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,flag)
      do 4 n=1,m
      yh(n)=yin(n)+h*f(n)
    4 continue
      x=xa+(k+k-1)*h
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,flag)
      do 5 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+h*hk(n)
    5 continue
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,flag)
      do 6 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+st*hk(n)
    6 continue
      x=xa+k*st
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,flag)
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
