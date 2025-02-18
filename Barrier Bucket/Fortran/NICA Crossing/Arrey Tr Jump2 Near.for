      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2),turn_p(2000),tau(2000),dE(2000),
     *dp_p(2000),dE_E(2000),eta_p(2000),E_p(2000),V_p(2000),
     *gamma_p(2000),g_tr_p(2000),al_0_p(2000),
     *al_1_p(2000),eta0_p(2000),eta1_p(2000),flag(2000),
     *dphi(2000),beta_p(2000),V_sc(2000),
     *Distr(100),bin_phi(100),bin_count(100)
      REAL*8 :: n_bin(100),A(6,7),coef(6),dcoef(5)
      REAL*8 :: k,Np,lp,lb,ll,deriv,dis
	COMMON/PARAM/T1,T2,Vacc,Vbb,
     *T0,dn,n_tr,dn_jump
      external func
      external Voltage
      OPEN (UNIT=13,file='results/Distr/21PhaseSpace.dat',
     *STATUS='unknown')
      OPEN (UNIT=14,file='results/Distr/21InitPhaseSpace.dat',
     *STATUS='unknown')
      OPEN (UNIT=15,file='results/Distr/21PhaseSpaceA.dat',
     *STATUS='unknown')

      OPEN (UNIT=16,file='results/Distr/21InDistribution.dat',
     *STATUS='unknown')
      OPEN (UNIT=17,file='results/Distr/21InTransform.dat',
     *STATUS='unknown')
      OPEN (UNIT=18,file='results/Distr/21Distribution.dat',
     *STATUS='unknown')
      OPEN (UNIT=19,file='results/Distr/21Transform.dat',
     *STATUS='unknown')

      OPEN (UNIT=20,file='results/Distr/21InDer.dat',
     *STATUS='unknown')
      OPEN (UNIT=21,file='results/Distr/21Der.dat',
     *STATUS='unknown')
      OPEN (UNIT=22,file='results/Distr/21DistributionAll.dat',
     *STATUS='unknown')

      OPEN (UNIT=23,file='results/Distr/21MNK.dat',
     *STATUS='unknown')
      OPEN (UNIT=24,file='results/Distr/21Deriv.dat',
     *STATUS='unknown')
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
      T1 = 40.0d0
      T2 = (T0 - 2*T1)/(2*2)
      Vacc = 300.0 * 1.D-6
      Vbb = 5000 * 1.D-6

c------parameters of Ring and Beam--------
      alpha_0 = 0.0199123143
      alpha_1 = -0.06371671224
      gamma_tr = 1/sqrt(alpha_0)
      g_tr = gamma_tr

      delta_g_tr = 0.045
      dgdt = 8.5d0
      A1 = 6.9756292530169*1.D-2
      A2 =-5.5818760781807*1.D-1

      Einj = 5704.0d0
      Etr = 5709.0d0
      E = Einj

      gamma = (E/proton)+1
      beta = dsqrt(1-(1/(gamma**2)))

      n_tr = INT(((g_tr-1)*proton-Einj)/Vacc)
      dn = INT(n_tr-(((g_tr-delta_g_tr-1)*proton-Einj)/Vacc))+1
      dn_jump = Real(INT((2*delta_g_tr)/(T0*dgdt*1.D-9)))
c------parameters of turns and particles--------
      Np = 2000.0d0
      Nt = 40.0d0
      Sturn = 1000.0d0
      Nturn = Sturn

      dp_p_max = 0.0033
      tau_max = T2/2
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
      write (14,*) 0,tau(I),dE(I),dp_p(I),dE_E(I),
     *0,E_p(I),0,0,0,
     *0,0,0,0,1,
     *dphi(I),0
      END DO

      k = 32.0d0
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
      harm = 2.0

      const = ((harm**2)*g0*Z0*c*el)/(2*R0) * total

      call hist(dphi,Np,k,particles,x1,x2,Distr,
     *bin_count,bin_phi)
      call ex(10.0d0,0,ll)

      turn=0.0
      I2 = 0
      I=0
c------turn by N-turn--------
      DO 2 I=1,Nt,1
      WRITE(*,*) I
      I2 = I2+1

      if(turn.lt.n_tr-dn) g_tr = gamma_tr
      if(turn.gt.n_tr+dn) g_tr = gamma_tr
      if((n_tr-dn.le.turn).and.(turn.le.n_tr-dn_jump/2))
     * g_tr = gamma + delta_g_tr
      if((n_tr+dn_jump/2.lt.turn).and. (turn.le.n_tr+dn))
     * g_tr = gamma -  delta_g_tr
      if((n_tr-dn_jump/2.lt.turn).and. (turn.le.n_tr+dn_jump/2))
     * g_tr = gamma_tr-(dgdt*T0*1.D-9*(turn-(n_tr)))

      if(turn.lt.n_tr-dn) al_1 = alpha_1
      if(turn.gt.n_tr+dn) al_1 = alpha_1
      if((n_tr-dn.le.turn).and.(turn.le.n_tr-dn_jump/2))
     *al_1 = A1*g_tr+A2
      if((n_tr+dn_jump/2.lt.turn).and. (turn.le.n_tr+dn))
     *al_1 = A1*g_tr+A2
      if((n_tr-dn_jump/2.lt.turn).and. (turn.le.n_tr+dn_jump/2))
     *al_1 = A1*g_tr+A2

      if(turn.lt.n_tr-dn) Nturn = Sturn
      if(turn.gt.n_tr+dn) Nturn = Sturn
      if((n_tr-2*Sturn.le.turn).and.(turn.le.n_tr+Sturn))
     * Nturn = Sturn

      gamma = (E/proton) + 1
      beta = dsqrt(1-(1/(gamma**2)))

      al_0 = 1 / (g_tr**2)
      eta0 = al_0 - (1/(gamma**2))
      ans = al_1 + ((3*(beta**2))/(2*(gamma**2)))
      eta1 = ans - eta0/(gamma**2)

      call hist(dphi,Np,k,particles,x1,x2,Distr,
     *bin_count,bin_phi)
      CALL Gram(32,5,bin_phi,Distr,A,ex)
      CALL Gauss(5,A,coef,dcoef)

      write(23,*) coef

      count = 0.0
      DO II=1,Np,1
      if((tau(II).lt.(-T1-(T2/2))).or.
     *(tau(II).gt.(T1+(T2/2)))) then
      count = count+1
      END IF
      END DO
      particles = Np - count

      I1 = 0
c------particles--------
      DO 3 J=1,Np,1
      I1 = I1+1

      call ddistr(dphi(I1),coef,6,dis,deriv)
      V_sc(I1) = const*deriv*(particles/Np)/(gamma_p(I1)**2)
      write(24,*) i,j,dphi(I1),dis,deriv,V_sc(I1)

      if((tau(I1).lt.(-T1-(T2/2))).or.
     *(tau(I1).gt.(T1+(T2/2)))) then
      flag(I1)=0
      END IF
      if (flag(I1).eq.0) then
      cycle
      end if

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

      eta_p(I1) = eta0_p(I1)+eta1_p(I1)*dp_p(I1)
c      eta_p(I1) = eta0_p(I1)
c------solution of equations--------
      xa = 0.0
      xb = Nturn

      call runkut(xa,xb,y,2,func,eps,kr,kmax,
     *E_p(I1),eta_p(I1),V_p(I1),gamma_p(I1),beta_p(I1),
     *eta0_p(I1),eta1_p(I1),dE_E(I1),dp_p(I1),
     *turn_p(I1),V_sc(I1))

      tau(I1) = y(1)
      dE(I1) = y(2)

      dphi(I1) = 2*pi*tau(I1)/(T2+2*T1)
      dE_E(I1) = dE(I1) / E_p(I1)
      dp_p(I1) = dE_E(I1) / (beta_p(I1)*beta_p(I1))

      write (15,*) turn_p(I1),tau(I1),dE(I1),dp_p(I1),dE_E(I1),
     *eta_p(I1),E_p(I1),V_p(I1),gamma_p(I1),g_tr_p(I1),
     *al_0_p(I1),al_1_p(I1),eta0_p(I1),eta1_p(I1),flag(I1),
     *dphi(I1),V_sc(I1)

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
     *dphi(I1),V_sc(I1)
      END DO

      call hist(dphi,Np,k,particles,x1,x2,Distr,
     *bin_count,bin_phi)
      A=0
      B=0

      DO j=1,k
      write(18,*) j,bin_count(j),Distr(j),bin_phi(j)
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
      if((n_tr-dn.le.turn).and.(turn.le.n_tr-dn_jump/2))
     *sign=1.0
      if((n_tr+dn_jump/2.lt.turn).and. (turn.le.n_tr+dn))
     *sign=-1.0
      if((n_tr-dn_jump/2.lt.turn).and. (turn.le.n_tr+dn_jump/2))
     *sign=0.0

      V=0.0
      if(-T1-T2/2.le.y(1).and.y(1).le.-T2/2)
     *V=-sign*Vbb
      if(-T2/2.lt.y(1).and.y(1).lt.T2/2)
     *V=0.0
      if( T2/2.le.y(1).and.y(1).le.T1+T2/2)
     *V=sign*Vbb

      f(1) = (2*eta*T0*y(2))/(beta*beta*E)
      f(2) = V
!      f(2) = V + (V_sc/1000000.0d0)
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

      subroutine ht(x,A,B,n)
      REAL*8, INTENT (IN) :: n
      INTEGER :: i,j,k
      REAL(16) :: twopi
      REAL(8), INTENT (INOUT), DIMENSION (INT(n,1)) :: x
      REAL(8), INTENT (INOUT), DIMENSION (INT(n/2,1)) :: A,B
      twopi = 6.28318d0

      A=0
      B=0
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

      DO i=0,k-1
      dis = dis +
     *A(i+1)*cos(twopi*i*x/k) + B(i+1)*sin(twopi*i*x/k)
      deriv = deriv +
     *(twopi*i/k)*((-1)*A(i+1)*sin(twopi*i*x/k)
     *+B(i+1)*cos(twopi*i*x/k))
      END DO

      RETURN
      END SUBROUTINE derivative

      subroutine ex(a,n,result)
      REAL*8, INTENT (IN) :: a
      INTEGER, INTENT (IN) :: n
      INTEGER :: i
      REAL*8 :: result,e

      e=1.0
      DO i=1,n
      e = e*a
      END DO
      result = e

      END subroutine ex

      subroutine Gram(n,m,x,f,A,rhs)
      INTEGER, INTENT (IN) :: n,m
      INTEGER :: i,j
      REAL*8 :: p,q,r,s,d,result
      REAL*8, INTENT (INOUT), DIMENSION (n) :: x,f
      REAL*8, INTENT (INOUT), DIMENSION (m+1,m+2)::A
      A=0.0d0

      DO j=0,m,1
      s=0.0
      r=0.0
      q=0.0
      DO i=0,n-1,1
      call rhs(x(i+1),j,p)
      s=s+p
      r=r+p*f(i+1)
      call rhs(x(i+1),m,d)
      q=q+p*d
      END DO
      A(1,j+1)=s
      A(j+1,m+1)=q
      A(j+1,m+1+1)=r
      END DO

      DO i=1,m,1
      DO j=0,m-1
      A(i+1,j+1)=A(i,j+2)
      END DO
      END DO

      RETURN
      END subroutine Gram

      subroutine Gauss(m,A,c,c1)
      INTEGER, INTENT (IN) :: m
      INTEGER :: i,j,j1,k,l,k1,n1
      real*8 :: r,s
      REAL*8, INTENT (INOUT), DIMENSION (m+1,m+2)::A
      REAL*8, INTENT (INOUT), DIMENSION (m+1)::c
      REAL*8, INTENT (INOUT), DIMENSION (m)::c1
      c=0.0

      n1=m+1
      DO k=0,m
      k1=k+1
      s=A(k+1,k+1)
      DO j=k1,n1
      A(k+1,j+1)=A(k+1,j+1)/s
      END DO

      DO i=k1,m
      r=A(i+1,k+1)
      DO j=k1,n1
      A(i+1,j+1)=A(i+1,j+1)-A(k+1,j+1)*r
      END DO

      END DO
      END DO

      DO i=m,0,-1
      s=a(i+1,n1+1)
      DO j=i+1,m
      s=s-A(i+1,j+1)*c(j+1)
      ENDDO
      c(i+1)=s
      ENDDO

      DO i=1,m
      c1(i) = i*c(i+1)
      END DO

      RETURN
      END subroutine Gauss


      subroutine ddistr(dphi,coef,m,result1,result2)
      INTEGER, INTENT (IN) :: m
      REAL*8, INTENT (IN) :: dphi
      INTEGER :: i
      REAL*8 :: result1,result2
      REAL*8, INTENT (INOUT), DIMENSION (m)::coef

      result1 = 0.0
      DO I=1,m
      result1 = result1+coef(i)*((dphi)**(i-1))
      END DO

      result2 =0.0
      DO I=1,m-1
      result2 = result2+(i)*coef(i+1)*(dphi)**(i-1)
      END DO

      RETURN
      END subroutine ddistr
