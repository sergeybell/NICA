      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2)
	COMMON/PARAM/T1,T2,Vacc,Vbb,
     *T0,dn,n_tr,dn_jump
      external func
      external Voltage
      OPEN (UNIT=13,file='PhaseSpace.dat',STATUS='unknown')
      OPEN (UNIT=14,file='InitPhaseSpace.dat',STATUS='unknown')

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

      Einj = 5500.0d0
      Etr = 5709.0d0

      gamma = (Einj/proton)+1
      beta = dsqrt(1-(1/(gamma**2)))
      E = Einj

      n_tr = INT(((g_tr-1)*proton-Einj)/Vacc)
      dn = INT(n_tr-(((g_tr-delta_g_tr-1)*proton-Einj)/Vacc))+1
      dn_jump = Real(INT((2*delta_g_tr)/(T0*dgdt*1.D-9)))
c      write(*,*) T0,n_tr,dn,dn_jump

c------parameters of turns and particles--------
      tau_max = T2/2+T1
      dE_max  = 21
      tau_max = 0.0
      dE_max  = 10.0
      Np = 1.0d0
      Nt = 200.0d0
      Sturn = 15000.0d0
c      Nt = 200.0d0
c      Sturn = 20000.0d0
      Nturn = Sturn

      I1 = 0
c------particles--------
      DO 3 II=1,Np
      I1 = I1+1
      E = Einj
      delta = 1.0d0/Np
      call random_number(rand1)
      call random_number(rand2)
      tau0 = tau_max * (2*rand1-1)
      dE0 = dE_max * (2*rand2-1)

      tau0 = tau_max
      dE0 = dE_max

      dE_E = dE0 / E
      dp_p = dE_E / (beta**2)
      ay = tau0
      by = dE0

      xa = 0.0
      xb = 0.0
      turn=0.0

      write(*,*) tau0,dE0,dp_p
      write (14,*) turn,tau0,dE0,dp_p,dE_E,eta,E,V,gamma,g_tr
c------turn by N-turn--------
      DO 2 I=1,Nt
      I2 = I2+1

      if(turn.lt.n_tr-dn) g_tr = gamma_tr
      if(turn.gt.n_tr+dn) g_tr = gamma_tr
      if((n_tr-dn.le.turn).and.(turn.le.n_tr-dn+dn_jump))
     * g_tr = gamma_tr-(dgdt*T0*1.D-9*(turn-(n_tr-dn)))
      if((n_tr-dn+dn_jump.lt.turn).and. (turn.le.n_tr+dn))
     * g_tr = gamma -  delta_g_tr

      if(turn.lt.n_tr-dn) al_1 = alpha_1
      if(turn.gt.n_tr+dn) al_1 = alpha_1
      if((n_tr-dn.le.turn).and.(turn.le.n_tr-dn+dn_jump))
     *al_1 = A1*g_tr+A2
      if((n_tr-dn+dn_jump.lt.turn).and. (turn.le.n_tr+dn))
     *al_1 = A1*g_tr+A2

c      eta0 = alpha_0 - (1/(gamma**2))
c      ans = alpha_1 + ((3*(beta**2))/(2*(gamma**2)))
c      eta1 = ans - eta0/(gamma**2)

      eta = eta0 + eta1 * dp_p

      if(turn.lt.n_tr-dn) Nturn = Sturn
      if(turn.gt.n_tr+dn) Nturn = Sturn
      if((n_tr-dn-(2*Sturn).le.turn)
     *.and.(turn.le.n_tr-dn+dn_jump+(2*Sturn)))
     * Nturn = Sturn/20
      if((n_tr-dn+dn_jump.lt.turn).and.(turn.le.n_tr+dn))
     * Nturn = Sturn


      E = E + (Vacc * Nturn)
      gamma = (E/proton) + 1
      beta = dsqrt(1-(1/(gamma**2)))

      al_0 = 1 / (g_tr**2)
      eta0 = al_0 - (1/(gamma**2))
      ans = al_1 + ((3*(beta**2))/(2*(gamma**2)))
      eta1 = ans - eta0/(gamma**2)
c------solution of equations--------
      b = Nturn
      xa = (I2-1)*b
      xb = I2 * b
      turn = turn+Nturn

      y(1) = ay
      y(2) = by

      call runkut(xa,xb,y,2,func,eps,kr,kmax,
     *E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn)
      ay=y(1)
      by=y(2)
      yy1=y(1)
      yy2=y(2)

      dE_E = y(2) / E
      dp_p = dE_E / (beta**2)

      flag=1
      if((yy1.lt.(-T1-(T2/2)-100)).or.(yy1.gt.(T1+(T2/2)+100))) flag=0
      write (13,*) turn,yy1,yy2,dp_p,dE_E,eta,E,V,gamma,g_tr,
     *al_0,al_1,eta0,eta1,flag
      if((yy1.lt.(-T1-(T2/2)-100)).or.(yy1.gt.(T1+(T2/2)+100))) exit

    2 continue
c      write (13,*) turn,yy1,yy2,dp_p,dE_E,eta,E,V,gamma,g_tr,
c     *al_0,al_1,eta0,eta1,flag
    3 continue
      stop
      end



      subroutine func(x,y,f,
     *E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn)
      implicit double precision (a-h,o-z)
      COMMON/PARAM/T1,T2,Vacc,Vbb,
     *T0,dn,n_tr,dn_jump
      dimension y(2),f(2)

      sign = 0.0
      if(turn.lt.n_tr-dn) sign=+1.0
      if(turn.gt.n_tr+dn) sign=-1.0
      if((n_tr-dn.le.turn).and.(turn.le.n_tr-dn+dn_jump))
     *sign=0.0
      if((n_tr-dn+dn_jump.lt.turn).and.(turn.le.n_tr+dn))
     *sign=-1.0

c      if(turn.lt.n_tr) sign = 1.0
c      if(turn.ge.n_tr) sign = -1.0

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
      call rhs(x,yin,f,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn)
      do 4 n=1,m
      yh(n)=yin(n)+h*f(n)
    4 continue
      x=xa+(k+k-1)*h
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn)
      do 5 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+h*hk(n)
    5 continue
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn)
      do 6 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+st*hk(n)
    6 continue
      x=xa+k*st
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn)
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
