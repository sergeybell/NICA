      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2),turn_p(1000),tau(10000),dE(10000),
     *dp_p(1000),dE_E(1000),eta_p(1000),E_p(1000),V_p(1000),
     *gamma_p(1000),g_tr_p(1000),al_0_p(1000),
     *al_1_p(1000),eta0_p(1000),eta1_p(1000),flag(1000),
     *dphi(1000),beta_p(1000), tau_r(1000),
     *Distr(1000),Density(1000),Der(1000),tau_bin(1000)
	COMMON/PARAM/T1,T2,Vacc,Vbb,
     *T0,dn,n_tr,dn_jump,Etr
      external func
      external Voltage
      OPEN (UNIT=14,file='results/InPS.dat',STATUS='unknown')
      OPEN (UNIT=16,file='results/Distr.dat',STATUS='unknown')

c------constants--------
      proton = 938.0d0
      pi = 3.14159
      e=2.718
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
      Np = 10.0d0

      dp_p_max = 0.005
      tau_max = T2/2+T1
      dE_max  = Einj*dp_p_max

      delta = 1.0/(Np-1)

      ro = 50.444d0
      mo = ANINT(ro+1)-DMOD(ro,10.0d0)
      Write(*,*) Np, Delta, tau_max,mo

      count=0.0
      DO IX=1,Np
      count = count+1.0
      call random_number(rand1)
      call random_number(rand2)
      call random_number(rand1)
      call random_number(rand2)

      xi = 2*rand1-1
      yi = 2*rand2-1

      s = xi**2+yi**2
      if (s .ge. 1.0) cycle
c-----from x = -1 to 1
c      x = 2*(count-1)*delta-1
c-----равномерное распределение
c      tau0 = tau_max*(x)
c-----какое-то распределение
c      tau0 = -tau_max*(-2*(x)**2+1)
c      tau0 = tau_max*((1.0/2)*(x+1)**2-1)
c-----экспонента
c      A = 2.0/(e-e**(-1))
c      C = -A/2*(e+e**(-1))
c      tau0 = tau_max * (A*e**x+C)
c      tau0 = -tau_max*1.0/sqrt(2.0) * sqrt(x+1)
c      tau0 = tau_max * e**x / e
      z0 = xi * sqrt(-2.0*DLOG(s)/s)
      z1 = yi * sqrt(-2.0*DLOG(s)/s)

      dE0 =0.0

      tau(IX) = z1*tau_max
      dE(IX) = dE0

      write (14,*) count,X,tau(IX),dE(IX),z0,z1
      END DO

      K=INT(100)

      x1 = -T2/2-T1
      x2 = T2/2+T1

      h = (x2-x1)/k

      DO I=1,K
      count = 0
      DO J=1,Np
      IF((x1 + (I-1)*h).le.tau(J).and.tau(J).lt.(x1 + (I-1)*h+h))
     *count = count+1
      IF((x1 + (I-1)*h).le.tau(J).and.tau(J).lt.(x1 + (I-1)*h+h))
     *tau_bin(J) = I
      END DO
      tau_r(I) = x1 + (I-1)*h+h/2
      Distr(I) = count / Np
      Density(I) = count / (Np*h)
      IF (I.ge.2)
     *Der(I) = (Density(I) - Density(I-1))/h
      WRITE(16,*) I,tau_r(I), Distr(I),Density(I),Der(I),count
      END DO

      DO I=1,Np
      WRITE(*,*) I,tau(I),tau_bin(I),Distr(tau_bin(I)),
     *Der(tau_bin(I))
      END DO

      stop
      end
