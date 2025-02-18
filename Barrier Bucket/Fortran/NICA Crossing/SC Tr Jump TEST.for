      program tesrun
c***** Integrating motion in separatrix****
      implicit double precision (a-h,o-z)
      dimension y(2),tau(1000),dE(1000)
	COMMON/PARAM/T1,T2,Vacc,Vbb,
     *T0,dn,n_tr,dn_jump
      external func
      external Voltage
      OPEN (UNIT=13,file='results/PhaseSpace.dat',STATUS='unknown')
      OPEN (UNIT=14,file='results/InitPhaseSpace.dat',STATUS='unknown')
      OPEN (UNIT=15,file='results/PhaseSpaceA.dat',STATUS='unknown')
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

      Einj = 5000.0d0
      Etr = 5709.0d0
      E = Einj

      gamma = (E/proton)+1
      beta = dsqrt(1-(1/(gamma**2)))

      n_tr = INT(((g_tr-1)*proton-Einj)/Vacc)
      dn = INT(n_tr-(((g_tr-delta_g_tr-1)*proton-Einj)/Vacc))+1
      dn_jump = Real(INT((2*delta_g_tr)/(T0*dgdt*1.D-9)))

c------parameters of turns and particles--------
      Np = 10.0d0
      Nt = 60.0d0
      Sturn = 20000.0d0
      Nturn = Sturn
      write(*,*) dE_max

      dp_p_max = 0.005
      tau_max = T2/2+T1
      dE_max  = Einj*dp_p_max

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
      end
