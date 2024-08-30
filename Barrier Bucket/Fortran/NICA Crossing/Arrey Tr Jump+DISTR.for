      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2),turn_p(1000),tau(1000),dE(1000),
     *dp_p(1000),dE_E(1000),eta_p(1000),E_p(1000),V_p(1000),
     *gamma_p(1000),g_tr_p(1000),al_0_p(1000),
     *al_1_p(1000),eta0_p(1000),eta1_p(1000),flag(1000),
     *dphi(1000),beta_p(1000),V_sc(1000),
     *Distr(100),n_bin(100),bin_count(100),bin_phi(100),
     *A(100),B(100)
      real*8 :: k,particles,x1,x2,pi,Np,dis,deriv
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

      Einj = 5600.0d0
      Etr = 5709.0d0
      E = Einj

      gamma = (E/proton)+1
      beta = dsqrt(1-(1/(gamma**2)))

      n_tr = INT(((g_tr-1)*proton-Einj)/Vacc)
      dn = INT(n_tr-(((g_tr-delta_g_tr-1)*proton-Einj)/Vacc))+1
      dn_jump = Real(INT((2*delta_g_tr)/(T0*dgdt*1.D-9)))
c------parameters of turns and particles--------
      Np = 128.0d0
      Nt = 1.0d0
      Sturn = 20000.0d0
      Nturn = Sturn

      dp_p_max = 0.003
      tau_max = T2/2+T1
      dE_max  = Einj*dp_p_max

      DO I=1,Np
      call random_number(rand1)
      call random_number(rand2)
      tau0 = tau_max * (2*rand1-1)
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

      k = 32.0d0
      particles = Np
      x1=-pi
      x2=+pi

      call hist(dphi,Np,k,particles,x1,x2,Distr,
     *n_bin,bin_count,bin_phi)
      call ht(dphi,A,B,k)

      DO j=1,k
      write(16,*) j,bin_count(j),Distr(j),bin_phi(j)
      END DO

      DO j=1,k/2+1
      write(17,*) j,A(j),B(j)
      END DO

      DO i=1,Np
      call derivative(A,B,k,dphi(i),dis,deriv)
      write(20,*) i,dphi(i),dis,deriv
      END DO

      turn=0.0
c------turn by N-turn--------
      DO 2 I=1,Nt

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

      if(turn.lt.n_tr-dn) Nturn = Sturn
      if(turn.gt.n_tr+dn) Nturn = Sturn
      if((n_tr-dn-2*Sturn.le.turn).and.(turn.le.n_tr-dn+dn_jump+Sturn))
     * Nturn = Sturn/20
      if((n_tr-dn+dn_jump.lt.turn).and.(turn.le.n_tr+dn))
     * Nturn = Sturn

      gamma = (E/proton) + 1
      beta = dsqrt(1-(1/(gamma**2)))

      al_0 = 1 / (g_tr**2)
      eta0 = al_0 - (1/(gamma**2))
      ans = al_1 + ((3*(beta**2))/(2*(gamma**2)))
      eta1 = ans - eta0/(gamma**2)

      I1 = 0
c------particles--------
      DO 3 J=1,Np-1
c      if((tau(J).lt.(-T1-(T2/2)-100)).or.
c     *(tau(J).gt.(T1+(T2/2)+100)))
c     *flag(J)=0
c      if (flag(J).eq.0.0) cycle

      turn_p(J) = turn
      E_p(J) = E
      gamma_p(J) = gamma
      beta_p(J) = beta
      g_tr_p(J) = g_tr
      al_0_p(J) = al_0
      al_1_p(J) = al_1
      eta0_p(J) = eta0
      eta1_p(J) = eta1

      y(1) = tau(J)
      y(2) = dE(J)

      eta_p(J) = eta0_p(J)+eta1_p(J)*dp_p(J)

c------solution of equations--------
      xa = 0.0
      xb = Nturn

      call runkut(xa,xb,y,2,func,eps,kr,kmax,
     *E_p(J),eta_p(J),V_p(J),gamma_p(J),beta_p(J),
     *eta0_p(J),eta1_p(J),dE_E(J),dp_p(J),turn_p(J),V_sc(J))

      tau(J) = y(1)
      dE(J) = y(2)

      dphi(J) = 2*pi*tau(J)/(T2+2*T1)
      dE_E(J) = dE(J) / E_p(J)
      dp_p(J) = dE_E(J) / (beta_p(J)*beta_p(J))

      write (15,*) turn_p(J),tau(J),dE(J),dp_p(J),dE_E(J),
     *eta_p(J),E_p(J),V_p(J),gamma_p(J),g_tr_p(J),
     *al_0_p(J),al_1_p(J),eta0_p(J),eta1_p(J),flag(J),
     *dphi(J),V_sc(J),I,J

    3 END DO
      WRITE(*,*) I,J
      turn = turn+Nturn
      E = E + (Vacc * (Nturn))
    2 END DO

      DO J=1,Np-1
      write (13,*) turn_p(J),tau(J),dE(J),dp_p(J),dE_E(J),
     *eta_p(J),E_p(J),V_p(J),gamma_p(J),g_tr_p(J),
     *al_0_p(J),al_1_p(J),eta0_p(J),eta1_p(J),flag(J),
     *dphi(J),V_sc(J),J
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
     *E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,V_sc)
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
      call rhs(x,yin,f,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,V_sc)
      do 4 n=1,m
      yh(n)=yin(n)+h*f(n)
    4 continue
      x=xa+(k+k-1)*h
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,V_sc)
      do 5 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+h*hk(n)
    5 continue
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,V_sc)
      do 6 n=1,m
      f(n)=f(n)+(hk(n)+hk(n))
      yh(n)=yin(n)+st*hk(n)
    6 continue
      x=xa+k*st
      call rhs(x,yh,hk,E,eta,V,gamma,beta,eta0,eta1,dE_E,dp_p,turn,V_sc)
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

      subroutine ht(x,A,B,n)
      REAL*8, INTENT (IN) :: n
      INTEGER :: i,j,k
      REAL(16) :: twopi
      REAL(8), INTENT (INOUT), DIMENSION (INT(n,1)) :: x
      REAL(8), INTENT (INOUT), DIMENSION (INT(n/2,1)) :: A,B
      twopi = 6.28318d0

      k=0
      DO k=0,n/2
      DO i=0,n-1
      IF (k.eq.0) then
      A(k+1) = A(k+1) + (1/n) * x(i+1) * cos(twopi*k*i/n)
      END IF
      IF (k.eq.n/2) then
      A(k+1) = A(k+1) + (1/n) * x(i+1) * cos(twopi*k*i/n)
      END IF
      IF ((k.gt.0).and.(k.lt.n/2)) then
      A(k+1) = A(k+1) + (2/n) * x(i+1) * cos(twopi*k*i/n)
      END IF
      B(k+1) = B(k+1) + (2/n) * x(i+1) * sin(twopi*k*i/n)
      END DO
      END DO

      RETURN
      END SUBROUTINE ht

      subroutine hist(tau,Np,k,particles,x1,x2,Distr,
     *n_bin,bin_count,bin_phi)
      REAL*8, INTENT (IN) :: Np,k,particles,x1,x2
      INTEGER :: i,j,count
      REAL*8 :: h,ik
      REAL(8), INTENT (INOUT), DIMENSION (INT(Np,1)) :: tau
      REAL(8), INTENT (INOUT), DIMENSION (INT(k,1)) :: Distr
      REAL*8:: n_bin(INT(Np,1))
      REAL*8:: bin_count(INT(k,1))
      REAL*8:: bin_phi(INT(k,1))

      h = (x2-x1)/k

      DO i=1,k
      count = 0
      DO j=1,Np
      IF((tau(j).lt.x1).or.(tau(j).gt.x2)) then
      n_bin(j) = -1
      cycle
      END IF
      IF(((x1 +(I-1)*h).le.tau(j)).and.(tau(j).lt.(x1+(I-1)*h+h)))
     *then
      count = count+1
      n_bin(J) = I
      END IF
      END DO
      bin_phi(I) = (x1+(I-1)*h)+h/2
      bin_count(I) = count
      Distr(I) = count / particles
      END DO

      RETURN
      END SUBROUTINE hist

      subroutine derivative(A,B,k,phi,dis,deriv)
      REAL*8, INTENT (IN) :: k,phi
      INTEGER :: i,j
      REAL*8 :: dis,deriv,pi,twopi
      REAL(8), INTENT (INOUT), DIMENSION (INT(k,1)) :: A,B
      pi=3.14159
      twopi = 2*pi
      x = (phi/pi * (k/2)) + k/2
      dis=0.0d0
      deriv=0.0d0
c      write(*,*) phi,x

      DO i=0,k-1
      dis = dis +
     *A(i+1)*cos(twopi*i*x/k) + B(i+1)*sin(twopi*i*x/k)
      deriv = deriv +
     *(twopi*i/k)*((-1)*A(i+1)*sin(twopi*i*x/k)
     *+B(i+1)*cos(twopi*i*x/k))
      END DO

      RETURN
      END SUBROUTINE derivative
