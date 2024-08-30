      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2),turn_p(1000),tau(1000),dE(1000),
     *dp_p(1000),dE_E(1000),eta_p(1000),E_p(1000),V_p(1000),
     *gamma_p(1000),g_tr_p(1000),al_0_p(1000),
     *al_1_p(1000),eta0_p(1000),eta1_p(1000),flag(1000),
     *dphi(1000),beta_p(1000),V_sc(1000),
     *Distr(1000),n_bin(1000),bin_count(1000),bin_phi(1000)
      REAL*8 :: A(1000),B(1000)
      real*8 :: k,Np,lp,lb
	COMMON/PARAM/T1,T2,Vacc,Vbb,
     *T0,dn,n_tr,dn_jump
      external func
      external Voltage
      OPEN (UNIT=13,file='results/PhaseSpace.dat',STATUS='unknown')
      OPEN (UNIT=14,file='results/InitPhaseSpace.dat',STATUS='unknown')
      OPEN (UNIT=15,file='results/PhaseSpaceA.dat',STATUS='unknown')

      OPEN (UNIT=16,file='results/InDistribution.dat',STATUS='unknown')
      OPEN (UNIT=17,file='results/InTransform.dat',STATUS='unknown')
      OPEN (UNIT=18,file='results/Distribution.dat',STATUS='unknown')
      OPEN (UNIT=19,file='results/Transform.dat',STATUS='unknown')

      OPEN (UNIT=20,file='results/InDer.dat',STATUS='unknown')
      OPEN (UNIT=21,file='results/Der.dat',STATUS='unknown')
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

      Einj = 5100.0d0
      Etr = 5709.0d0
      E = Einj

      gamma = (E/proton)+1
      beta = dsqrt(1-(1/(gamma**2)))

      n_tr = INT(((g_tr-1)*proton-Einj)/Vacc)
      dn = INT(n_tr-(((g_tr-delta_g_tr-1)*proton-Einj)/Vacc))+1
      dn_jump = Real(INT((2*delta_g_tr)/(T0*dgdt*1.D-9)))
c------parameters of turns and particles--------
      Np = 1000.0d0
      Nt = 1.0d0
      Sturn = 20000.0d0
      Nturn = Sturn

      dp_p_max = 0.005
      tau_max = T2/2
      dE_max  = Einj*dp_p_max

      DO I=1,Np
      call random_number(rand1)
      call random_number(rand2)
      tau0 = tau_max * ((2/(Np-1))*I-(Np+1)/(Np-1))
      write(*,*) ((2/(Np-1))*I-(Np+1)/(Np-1))
      dE0 = dE_max * (2*rand2-1)

      flag(I) = 1
      E_p(I) = E
      beta_p(I) = beta
      gamma_p(I) = gamma
      tau(I) = tau0
      dE(I) = dE0
      dphi(I) = 2*pi*tau(I)/(T2+2*T1)
      dE_E(I) = dE(I) / E_p(I)
      dp_p(I) = dE_E(I) / (beta_p(I)**2)
      write (14,*) I,tau(I),dE(I)
      END DO

      k = 64.0d0
      particles = Np
      x1=-pi
      x2=+pi
      R0 = Cring/(2*pi)
      Z0 = 377
      lp = 5.0
      lb = 0.55
      g0 = 1+2*log(lp/lb)
      el = 1.6*1.d-19
      total = 1*1.D13

      const = (g0*Z0*c*el)/(2*R0) * total

      call hist(dphi,Np,k,particles,x1,x2,Distr,
     *bin_count,bin_phi)

      write(*,*) k
      DO j=1,k
      write(16,*) j,bin_count(j),Distr(j),bin_phi(j)
      END DO

      stop
      end


      subroutine func(x,y,f,
     *E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,V_sc)
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
      call rhs(x,yin,f,E,eta,V,gamma,beta,eta0,eta1,dE_E,
     *dp_p,turn,V_sc)
      do 4 n=1,m
      yh(n)=yin(n)+h*f(n)
    4 continue
      x=xa+(k+k-1)*h
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,
     *dp_p,turn,V_sc)
      do 5 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+h*hk(n)
    5 continue
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,
     *dp_p,turn,V_sc)
      do 6 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+st*hk(n)
    6 continue
      x=xa+k*st
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,
     *dp_p,turn,V_sc)
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

      subroutine hist(tau,Np,k,particles,x1,x2,Distr,
     *bin_count,bin_phi)
      REAL*8, INTENT (IN) :: Np,k,particles,x1,x2
      INTEGER :: ik,j,count
      REAL*8 :: h,count1
      REAL(8), INTENT (INOUT), DIMENSION (INT(Np,2)) :: tau
      REAL(8), INTENT (INOUT), DIMENSION (INT(k,2)) :: Distr
      REAL*8:: n_bin(INT(Np,2))
      REAL*8:: bin_count(INT(k,2))
      REAL*8:: bin_phi(INT(k,2))


      h = (x2-x1)/k
      write(*,*) Np,k,particles,x1,x2,INT(Np,2),INT(k,1)
      count1=0.0
      DO ik=1,k
      count = 0
      count1=0.0
      DO j=1,Np
      IF((tau(j).lt.x1).or.(tau(j).gt.x2)) then
      cycle
      END IF
      IF(((x1 +(ik-1)*h).le.tau(j)).and.
     *(tau(j).lt.(x1+(ik-1)*h+h))) then
      count = count+1
      END IF
      count1 = count1+1
      END DO
      bin_phi(ik) = (x1+(ik-1)*h)+h/2
      bin_count(ik) = count
      Distr(ik) = count / count1
      END DO
      write(*,*) k

      RETURN
      END SUBROUTINE hist


      subroutine hist1(tau,Np,k,particles,x1,x2,Distr,
     *n_bin,bin_count,bin_phi)
      REAL*8, INTENT (IN) :: Np,k,particles,x1,x2
      INTEGER :: ik,j,count
      REAL*8 :: h,count1
      REAL(8), INTENT (INOUT), DIMENSION (INT(Np,1)) :: tau
      REAL(8), INTENT (INOUT), DIMENSION (INT(k,1)) :: Distr
      REAL*8:: n_bin(INT(Np,1))
      REAL*8:: bin_count(INT(k,1))
      REAL*8:: bin_phi(INT(k,1))

      h = (x2-x1)/k

      count1=0.0
      DO ik=1,k
      count = 0
      count1=0.0
      DO j=1,Np
      IF((tau(j).lt.x1).or.(tau(j).gt.x2)) then
      n_bin(j) = -1
      cycle
      END IF
      IF(((x1 +(ik-1)*h).le.tau(j)).and.(tau(j).le.(x1+(ik-1)*h+h)))
     *then
      count = count+1
      n_bin(j) = ik
      END IF
      count1 = count1+1
      END DO
      bin_phi(ik) = (x1+(ik-1)*h)+h/2
      bin_count(ik) = count
      Distr(ik) = count / count1
      END DO

      RETURN
      END SUBROUTINE hist1
