c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: SUSY_FLAVOR.F
c     Released: 15:02:2010(J.R.)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Driver file for SUSY_FLAVOR library                             c
c     Example of MSSM parameter initialization                        c
c     Test output for SUSY spectrum and implemented rare decays       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program susy_flavor
      implicit double precision (a-h,o-z)
      dimension sll(3),slr(3),amsq(3),amsu(3),amsd(3)
      double complex asl(3),asu(3),asd(3)
      double complex slmi_l(3),slmi_r(3),slmi_lr(3,3)
      double complex sqmi_l(3),sdmi_r(3),sumi_r(3)
      double complex sdmi_lr(3,3),sumi_lr(3,3)
      double complex amg,amgg,amue
      common/sf_cont/eps,indx(3,3),iconv

c     decide if input parameters are read from file susy_flavor.in or
c     defined directly inside the program
      write(*,'(a,$)')
     $     'Read input from file susy_flavor.in (no=1,yes=2)? '
      read(*,*) input_type
      if (input_type.eq.2) then
         call sflav_input       ! Parameters read from file susy_flavor.in 
         goto 100
      end if

c     Parameters defined inside the code.  

c     Input parameters convention choice
c     iconv = 1                 ! SLHA2 input conventions
      iconv = 2                 ! hep-ph/9511250 input conventions

c     SM basic input initialization
      zm0 = 91.1876d0           ! M_Z
      wm0 = 80.398d0            ! M_W
      alpha_z = 1/127.934d0     ! alpha_em(M_Z)
      call vpar_update(zm0,wm0,alpha_z)

c     QCD parameters
      alpha_s = 0.1172d0        ! alpha_s(MZ)
      call lam_fit(alpha_s)     ! fits Lambda_QCD at 3 loop level
      call lam_fit_nlo(alpha_s) ! fits Lambda_QCD at NLO level

c     CKM matrix initialization
      alam = 0.2258d0           ! lambda
      apar = 0.808d0            ! A
      rhobar = 0.177d0          ! rho bar
      etabar = 0.360d0          ! eta bar
      call ckm_wolf(alam,apar,rhobar,etabar)

c     Fermion mass initialization, input: MSbar running quark masses
      top_scale = 163.2d0
      top = 163.2d0             ! m_t(top_scale)
      bot_scale = 4.17d0
      bot = 4.17d0              ! m_b(bot_scale)
      call init_fermion_sector(top,top_scale,bot,bot_scale)

c     Higgs sector parameters
      pm    = 200               ! M_A
      tanbe = 10                ! tan(beta)
      amue  = (200.d0,100.d0)   ! mu
      call init_higgs_sector(pm,tanbe,amue,ierr)
      if (ierr.ne.0) stop 'negative tree level Higgs mass^2?'

c     Gaugino sector parametersc. CAUTION: if M1 is set to 0 here then
c     program sets M1 and M2 GUT-related, i.e. M1 = 5/3 s_W^2/c_W^2*M2
      amgg  = 0.d0              ! M1 (bino mass)
      amg   = (200.d0,0.d0)     ! M2 (wino mass)
      amglu = 3*abs(amg)        ! M3 (gluino mass)
      call init_ino_sector(amgg,amg,amglu,amue,tanbe,ierr)
      if (ierr.ne.0) write(*,*) '-ino mass below M_Z/2?'

c     Slepton diagonal soft breaking parameters
      sll(1) = 300.d0           ! left selectron mass scale
      sll(2) = 300.d0           ! left smuon mass scale
      sll(3) = 300.d0           ! left stau mass scale
      slr(1) = 300.d0           ! right selectron mass scale
      slr(2) = 300.d0           ! right smuon mass scale
      slr(3) = 300.d0           ! right stau mass scale
c     Dimensionless (normalized to masses) slepton diagonal LR mixing  
      asl(1) = (1.d0,0.d0)      ! 1st generation 
      asl(2) = (1.d0,0.d0)      ! 2nd generation 
      asl(3) = (1.d0,0.d0)      ! 3rd generation 
c     Slepton LL and RR mass insertions (hermitian matrices)
c     slmi_x(1),slmi_x(2), slmi_x(3) are 12,23,31 entry respectively
      do i=1,3
         slmi_l(i) = (0.d0,0.d0) ! slepton LL mass insertion
         slmi_r(i) = (0.d0,0.d0) ! slepton RR mass insertion
      end do 
      slmi_l(2) = (2.d-2,1.d-2) ! example, non-vanishing LL 23 entry
c     Slepton LR mass insertions, non-hermitian in general
      do i=1,3
         do j=1,3
            slmi_lr(i,j) = (0.d0,0.d0) ! slepton LR ij mass insertion
         end do
      end do 
c     Calculate physical masses and mixing angles
      call init_slepton_sector(sll,slr,asl,ierr,slmi_l,slmi_r,slmi_lr)
      if (ierr.ne.0) stop 'negative tree level slepton mass^2?'

c     Squark diagonal soft breaking parameters
      amsq(1) = 500.d0          ! left squark mass, 1st generation
      amsq(2) = 500.d0          ! left squark mass, 2nd generation
      amsq(3) = 400.d0          ! left squark mass, 3rd generation
      amsd(1) = 550.d0          ! right down squark mass
      amsd(2) = 550.d0          ! right strange squark mass
      amsd(3) = 300.d0          ! right sbottom mass
      amsu(1) = 450.d0          ! right up squark mass
      amsu(2) = 450.d0          ! right charm squark mass
      amsu(3) = 200.d0          ! right stop mass
c     Dimensionless (normalized to masses) squark diagonal LR mixing  
      asd(1) = (1.d0,0.d0)      ! down squark LR mixing, 1st generation
      asd(2) = (1.d0,0.d0)      ! down squark LR mixing, 2nd generation
      asd(3) = (1.d0,0.d0)      ! down squark LR mixing, 3rd generation
      asu(1) = (1.d0,0.d0)      ! up squark LR mixing, 1st generation
      asu(2) = (1.d0,0.d0)      ! up squark LR mixing, 2nd generation
      asu(3) = (1.d0,0.d0)      ! up squark LR mixing, 3rd generation
c     Squark LL and RR mass insertions (hermitian matrices)
c     sqmi_l(1),sqmi_l(2), sqmi_l(3) are 12,23,31 entry respectively, etc.
      do i=1,3
         sqmi_l(i) = (0.d0,0.d0) ! squark LL mass insertion
         sumi_r(i) = (0.d0,0.d0) ! up-squark RR mass insertion
         sdmi_r(i) = (0.d0,0.d0) ! down-squark RR mass insertion
      end do 
      sqmi_l(2) = (2.d-2,-1.d-2) ! example, non-vanishing LL 23 entry
c     Squark LR mass insertions, non-hermitian in general
      do i=1,3
         do j=1,3
            sumi_lr(i,j) = (0.d0,0.d0) ! up-squark LR ij mass insertion
            sdmi_lr(i,j) = (0.d0,0.d0) ! down-squark LR ij mass insertion
         end do
      end do 
c     Calculate physical masses and mixing angles
      call init_squark_sector(amsq,amsu,amsd,asu,asd,ierr,sqmi_l,sumi_r,
     $     sdmi_r,sumi_lr,sdmi_lr)
      if (ierr.ne.0) stop 'negative tree level squark mass^2?'

c     reset status of physical Higgs mass after parameter changes
      call reset_phys_data
c     Neutral CP-even Higgs masses in the 1-loop Effective Potential
c     Approximation. Only real mu, A_t, A_b allowed - replaced x->abs(x)
      call fcorr_EPA(tanbe,pm,top,abs(amue),amsq(3),amsd(3),amsu(3)
     $     ,abs(asd(3)),abs(asu(3)),ierr)
      if (ierr.ne.0) stop 'negative 1-loop EPA CP-even Higgs mass^2?'

c     !!! End of input section !!!
 100  continue

c     Control output
      write(*,99)'************************************************'
      write(*,99)'MSSM Lagrangian parameters and tree level masses'
      write(*,99)'written on file mssm_data.txt'
      write(*,99)'************************************************'
      ifl = 1                   ! output file number
      open(ifl,file='mssm_data.txt',status='unknown')
      call print_MSSM_par(ifl)  ! Lagrangian parameters
      call print_MSSM_masses(ifl) ! tree level physical masses
      close(ifl)

c     Results for implemented observables:
      write(*,*)
      write(*,99)'Physical observables:'
      write(*,*)
      write(*,99)'Electric dipole moments:'
      write(*,99)'Electron EDM = ',edm_l(1)
      write(*,99)'Muon EDM     = ',edm_l(2)
      write(*,99)'Tau EDM      = ',edm_l(3)
      write(*,99)'Neutron EDM  = ',edm_n()
      write(*,*)

      write(*,99)'Neutrino K decays:'
      call k_pivv(br_k0,br_kp)
      write(*,99)'BR(K_L^0 -> pi^0 vv) = ',br_k0
      write(*,99)'BR(K^+   -> pi^+ vv) = ',br_kp
      write(*,*)

      write(*,99)'Leptonic B decays:'
      write(*,99)'BR(B_d -> mu^+ mu^-) = ',b_ll(3,1,2,2)
      write(*,99)'BR(B_s -> mu^+ mu^-) = ',b_ll(3,2,2,2)
      write(*,*)

      write(*,99)'B-> X_s photon decay:'
c     Physical quantities for BR(B->X_s g) calculation
      delb = 0.99d0             ! Photon energy infrared cutoff
      amiu_b= 4.8d0             ! Renormalization scale miu_b
      write(*,99)'BR(B -> X_S gamma) = ',bxg_nl(delb,amiu_b)
      write(*,*)

      write(*,99)'KK mixing:'
      call dd_kaon(eps_k,delta_mk)
      write(*,99)'eps_K = ',eps_k
      write(*,99)'Delta m_K = ',delta_mk
      write(*,*)

      write(*,99)'DD mixing:'
      call uu_dmeson(delta_md)
      write(*,99)'Delta m_D = ',delta_md
      write(*,*)

      write(*,99)'BB mixing:'
      call dd_bmeson(1,delta_mbd)
      write(*,99)'Delta m_B_d = ',delta_mbd
      call dd_bmeson(2,delta_mbs)
      write(*,99)'Delta m_B_s = ',delta_mbs
      write(*,*)

 99   format(a,1pe11.4)
      end
