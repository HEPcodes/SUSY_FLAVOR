c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: SFLAV_IO.F
c     Released: 20:02:2010(J.R.)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     SLHA2 based input routine for SUSY_FLAVOR                       c
c     Test output routines for MSSM parameters and SUSY spectrum      c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine skip_comments(ifl)
c     read the ifl file until the first character in line is blank space
      character*1 line
 10   read(ifl,'(a1)',END=100)line
      if (line(1:1).ne.' ') goto 10
      backspace(ifl)
      return
 100  stop 'Incomplete or corrupted SUSY_FLAVOR.SLHA input file'
      end

      subroutine sflav_input
c     Input read from file susy_flavor.in
      implicit double precision (a-h,o-z)
      double complex amg,amgg,amue
      double complex ls,ks,ds,es,us,ws
      double complex lms,rms,ums,dms,qms
      character*14 line
      logical sldiag,sqdiag
      dimension tmpi(3,3),tmpr(3,3)
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/vev/v1,v2
      common/sf_cont/eps,indx(3,3),iconv
c     open input file
      ifl = 1
      open(ifl,file='susy_flavor.in',status='old')
c     Find non-standard Block SOFTINP, check the status of soft parameters
 21   read(ifl,'(a14)',END=100)line
      if (line(7:13).ne.'SOFTINP') goto 21
      call skip_comments(ifl)
      read(ifl,*)k,iconv        ! soft parameters convention
c     iconv = 1:
c          sfermion input parameters in conventions of hep-ph/9511250
c     iconv = 2:
c          sfermion input parameters in SLHA2 conventions
      if ((iconv.ne.1).and.(iconv.ne.2)) stop
     $     'susy_flavor.in: incorrect iconv value'
      read(ifl,*)k,input_type   ! dimension of soft parameters 
c     input_type = 1:
c          sfermion off-diagonal terms given as mass insertions
c          LR diagonal terms given as dimensionless parameters
c     input_type = 2:
c          fermion soft terms given as absolute values
c     for more information read comments in susy_flavor.in file
      if ((input_type.ne.1).and.(input_type.ne.2)) stop
     $     'susy_flavor.in: incorrect input_type value'
c     find SM input block, read the SM data
 10   read(ifl,'(a14)',END=100)line
      if (line(7:14).ne.'SMINPUTS') goto 10
      call skip_comments(ifl)
      read(ifl,*)k,alpha_z      ! 1/alpha_em(M_Z)
      alpha_z = 1/alpha_z
      read(ifl,*)k,alpha_s      ! alpha_s(MZ)  
      read(ifl,*)k,zm0          ! M_Z
c     Fermion mass initialization, input: MSbar running quark masses
      read(ifl,*)k,bottom       ! mb(mb)
      bot_scale = bottom
      read(ifl,*)k,top          ! mt(mt)
      top_scale = top
      read(ifl,*)k,em(3)        ! m_tau (pole)
      read(ifl,*)k,em(1)        ! m_e (pole)
      read(ifl,*)k,em(2)        ! m_mu (pole)
      read(ifl,*)k,dml(1)       ! md(2 GeV)
      read(ifl,*)k,uml(1)       ! mu(2 GeV)
      read(ifl,*)k,dml(2)       ! ms(2 GeV)
      read(ifl,*)k,uml(2)       ! mc(2 GeV)
      read(ifl,*)k,wm0          ! M_W (not standard SLHA2!)
c     Electroweak and strong parameter initialization
      call vpar_update(zm0,wm0,alpha_z) ! sets electroweak parameters
      CALL lam_fit(alpha_s)     ! fits Lambda_QCD at 3 loop level 
      CALL lam_fit_NLO(alpha_s) ! fits Lambda_QCD at NLO level
      call init_fermion_sector(top,top_scale,bottom,bot_scale)
c     find V_CKM input block, read the CKM data
      rewind(ifl)
 11   read(ifl,'(a14)',END=100)line
      if (line(7:12).ne.'VCKMIN') goto 11
      call skip_comments(ifl)
      read(ifl,*)k,alam         ! lambda
      read(ifl,*)k,apar         ! A
      read(ifl,*)k,rhobar       ! rho bar
      read(ifl,*)k,etabar       ! eta bar
      call ckm_wolf(alam,apar,rhobar,etabar)
c     find EXTPAR input block, read the real Higgs and gaugino data
      rewind(ifl)
 12   read(ifl,'(a14)',END=100)line
      if (line(7:12).ne.'EXTPAR') goto 12
      call skip_comments(ifl)
      read(ifl,*)k,x1           ! M1 (bino mass, complex)
      read(ifl,*)k,x2           ! M2 (wino mass, complex)
      read(ifl,*)k,amglu        ! M3 (gluino mass)
      read(ifl,*)k,x3           ! mu (complex)
      read(ifl,*)k,tanbe        ! tan(beta)
      read(ifl,*)k,pm           ! M_A
c     find IMEXTPAR input block, read the imaginary Higgs and gaugino data
      rewind(ifl)
 121  read(ifl,'(a14)',END=100)line
      if (line(7:14).ne.'IMEXTPAR') goto 121
      call skip_comments(ifl)
      read(ifl,*)k,y1           ! M1 (bino mass, complex)
      amgg = dcmplx(x1,y1)
      read(ifl,*)k,y2           ! M2 (wino mass, complex)
      amg = dcmplx(x2,y2)
      read(ifl,*)k,y3           ! mu (complex)
      amue = dcmplx(x3,y3)
c     Higgs sector
      call init_higgs_sector(pm,tanbe,amue,ierr)
      if (ierr.ne.0) stop 'negative tree level Higgs mass^2?'
c     SUSY fermion sector 
      call init_ino_sector(amgg,amg,amglu,amue,tanbe,ierr)
      if (ierr.ne.0) write(*,*) '-ino mass below M_Z/2?'
      
c     find slepton soft input blocks, read the slepton data
c     left and right slepton mass
      rewind(ifl)
 13   read(ifl,'(a14)',END=100)line
      if (line(7:12).ne.'MSL2IN') goto 13
      do i=1,6
         read(ifl,*)k,l,tmpr(k,l) ! left slepton mass M_SL^2, real part
         if (k.eq.l) lms(k,l) = dcmplx(tmpr(k,l),0.d0)
      end do
      rewind(ifl)
 131  read(ifl,'(a14)',END=100)line
      if (line(7:14).ne.'IMMSL2IN') goto 131
      do i=1,3
         read(ifl,*)k,l,tmpi(k,l) ! left slepton mass M_SL^2, imaginary part
         lms(k,l) = dcmplx(tmpr(k,l),tmpi(k,l))
      end do
      rewind(ifl)
 14   read(ifl,'(a14)',END=100)line
      if (line(7:12).ne.'MSE2IN') goto 14
      do i=1,6
         read(ifl,*)k,l,tmpr(k,l) ! right slepton mass M_ER^2, real part
         if (k.eq.l) rms(k,l) = dcmplx(tmpr(k,l),0.d0)
      end do
      rewind(ifl)
 141  read(ifl,'(a14)',END=100)line
      if (line(7:14).ne.'IMMSE2IN') goto 141
      do i=1,3
         read(ifl,*)k,l,tmpi(k,l) ! right slepton mass M_ER^2, imaginary part
         rms(k,l) = dcmplx(tmpr(k,l),tmpi(k,l))
      end do
      if (input_type.eq.1) then
c     remove instabilities by adding tiny mass splitting
         do i=1,3
            lms(i,i) = (1 + eps*i)*lms(i,i)
            rms(i,i) = (1 - eps*i)*rms(i,i)
         end do
c     expand delta parameters to full mass entries and make h.c.
         do k=1,2
            do l=k+1,3
               lms(k,l) = lms(k,l)*sqrt(lms(k,k)*lms(l,l))
               lms(l,k) = dconjg(lms(k,l))
               rms(k,l) = rms(k,l)*sqrt(rms(k,k)*rms(l,l))
               rms(l,k) = dconjg(rms(k,l))
            end do
         end do
      end if
c     slepton LR mixing
      rewind(ifl)
 15   read(ifl,'(a14)',END=100)line
      if (line(7:10).ne.'TEIN') goto 15
      do i=1,9
         read(ifl,*)k,l,tmpr(k,l) ! slepton LR mixing A_L, real part
      end do
      rewind(ifl)
 151  read(ifl,'(a14)',END=100)line
      if (line(7:12).ne.'IMTEIN') goto 151
      do i=1,9
         read(ifl,*)k,l,tmpi(k,l) ! slepton LR mixing A_L, imaginary part
         ls(k,l) = dcmplx(tmpr(k,l),tmpi(k,l))
         ks(k,l) = (0.d0,0.d0)  ! non-holomorphic terms set to zero
      end do
      if (input_type.eq.1) then
c     calculate full LR mixing 
         do k=1,3
            do l=1,3
               if (k.ne.l) then
                  ls(k,l) = sqrt(abs(lms(k,k)*rms(l,l)))*sq2/v1*ls(k,l)
               end if  
            end do
            ls(k,k) = yl(k)*ls(k,k)*(lms(k,k)*rms(k,k))**0.25d0
         end do
      end if
c     if slepton input data in SLHA format, rewrite them to
c     hep-ph/9511250 convention
      if (iconv.eq.1) call sl_slha_to_jr
c     slepton diagonalization routine
      if (sldiag()) stop 'negative tree level slepton mass^2?'

c     find squark soft input blocks, read the squark data
c     left and right squark mass
      rewind(ifl)
 16   read(ifl,'(a14)',END=100)line
      if (line(7:12).ne.'MSQ2IN') goto 16
      do i=1,6
         read(ifl,*)k,l,tmpr(k,l) ! Left squark mass M_QL^2, real part
         if (k.eq.l) qms(k,l) = dcmplx(tmpr(k,l),0.d0)
      end do
      rewind(ifl)
 161  read(ifl,'(a14)',END=100)line
      if (line(7:14).ne.'IMMSQ2IN') goto 161
      do i=1,3
         read(ifl,*)k,l,tmpi(k,l) ! Left squark mass M_QL^2, imaginary part
         qms(k,l) = dcmplx(tmpr(k,l),tmpi(k,l))
      end do
      rewind(ifl)
 17   read(ifl,'(a14)',END=100)line
      if (line(7:12).ne.'MSU2IN') goto 17
      do i=1,6
         read(ifl,*)k,l,tmpr(k,l) ! right up-squark mass M_UR^2, real part
         if (k.eq.l) ums(k,l) = dcmplx(tmpr(k,l),0.d0)
       end do
      rewind(ifl)
 171  read(ifl,'(a14)',END=100)line
      if (line(7:14).ne.'IMMSU2IN') goto 171
      do i=1,3
         read(ifl,*)k,l,tmpi(k,l) ! right up-squark mass M_UR^2, imaginary part
         ums(k,l) = dcmplx(tmpr(k,l),tmpi(k,l))
       end do
      rewind(ifl)
 18   read(ifl,'(a14)',END=100)line
      if (line(7:12).ne.'MSD2IN') goto 18
      do i=1,6
         read(ifl,*)k,l,tmpr(k,l) ! right down-squark mass M_DR^2, real part 
         if (k.eq.l) dms(k,l) = dcmplx(tmpr(k,l),0.d0)
      end do
      rewind(ifl)
 181  read(ifl,'(a14)',END=100)line
      if (line(7:14).ne.'IMMSD2IN') goto 181
      do i=1,3
         read(ifl,*)k,l,tmpi(k,l) ! right down-squark mass M_DR^2, imaginary part
         dms(k,l) = dcmplx(tmpr(k,l),tmpi(k,l))
      end do
      if (input_type.eq.1) then
c     remove instabilities by adding tiny mass splitting
         do i=1,2
            qms(i,i) = (1 + eps*i)*dble(qms(i,i))
            ums(i,i) = (1 - eps*i)*dble(ums(i,i))
            dms(i,i) = (1 - eps*i)*dble(dms(i,i))
         end do
c     expand delta parameters to full mass entries and make h.c.
         do k=1,2
            do l=k+1,3
               qms(k,l) = qms(k,l)*sqrt(qms(k,k)*qms(l,l))
               qms(l,k) = dconjg(qms(k,l))
               ums(k,l) = ums(k,l)*sqrt(ums(k,k)*ums(l,l))
               ums(l,k) = dconjg(ums(k,l))
               dms(k,l) = dms(k,l)*sqrt(dms(k,k)*dms(l,l))
               dms(l,k) = dconjg(dms(k,l))
            end do
         end do
      end if
c     squark LR mixing
      rewind(ifl)
 19   read(ifl,'(a14)',END=100)line
      if (line(7:10).ne.'TUIN') goto 19
      do i=1,9
         read(ifl,*)k,l,tmpr(k,l) ! up-squark LR mixing A_U, real part
       end do
      rewind(ifl)
 191  read(ifl,'(a14)',END=100)line
      if (line(7:12).ne.'IMTUIN') goto 191
      do i=1,9
         read(ifl,*)k,l,tmpi(k,l) ! up-squark LR mixing A_U, imaginary part
         us(k,l) = dcmplx(tmpr(k,l),tmpi(k,l))
         ws(k,l) = (0.d0,0.d0)
      end do
      rewind(ifl)
 20   read(ifl,'(a14)',END=100)line
      if (line(7:10).ne.'TDIN') goto 20
      do i=1,9
         read(ifl,*)k,l,tmpr(k,l) ! down-squark LR mixing A_D, real part
      end do
      rewind(ifl)
 201  read(ifl,'(a14)',END=100)line
      if (line(7:12).ne.'IMTDIN') goto 201
      do i=1,9
         read(ifl,*)k,l,tmpi(k,l) ! down-squark LR mixing A_D, imaginary part
         ds(k,l) = dcmplx(tmpr(k,l),tmpi(k,l))
         es(k,l) = (0.d0,0.d0)
      end do
c     store dimensionless |A_t|, |A_b| for EPA Higgs mass calculation
      if (input_type.eq.1) then 
         yts = abs(us(3,3))
         ybs = abs(ds(3,3))
      else
         yts = abs(us(3,3))/(qms(3,3)*ums(3,3))**0.25d0
         ybs = abs(ds(3,3))/(qms(3,3)*dms(3,3))**0.25d0
      end if
c     calculate full LR mixing (non-holomorphic terms set to zero!)
      if (input_type.eq.1) then
         do k=1,3
            do l=1,3
               if (k.ne.l) then
                  ds(k,l) = sqrt(abs(qms(k,k)*dms(l,l)))*sq2/v1*ds(k,l)
                  us(k,l) = sqrt(abs(qms(k,k)*ums(l,l)))*sq2/v2*us(k,l)
               end if
            end do
            ds(k,k) = yd(k)*ds(k,k)*(qms(k,k)*dms(k,k))**0.25d0
            us(k,k) = yu(k)*us(k,k)*(qms(k,k)*ums(k,k))**0.25d0
         end do
      end if 
c     if squark input data in SLHA format, rewrite them to
c     hep-ph/9511250 convention
      if (iconv.eq.1) call sq_slha_to_jr
c     squark diagonalization routine
      if (sqdiag()) stop 'negative tree level squark mass^2?'
c     reset status of physical Higgs mass after parameter changes
      call reset_phys_data                  
c     Neutral CP-even Higgs masses in the simple 1-loop Effective
c     Potential Approximation (EPA). Only real mu, A_t, A_b allowed
      stbl = sqrt(abs(qms(3,3)))
      str = sqrt(abs(ums(3,3)))
      sbr = sqrt(abs(dms(3,3)))
      call fcorr_EPA(tanbe,pm,top,abs(amue),stbl,sbr,str,ybs,yts,ierr)
      if (ierr.ne.0) stop 'negative 1-loop EPA CP-even Higgs mass^2?'
      close(ifl)
      return
c     incorrect input file?
 100  stop 'Incomplete or corrupted SUSY_FLAVOR.SLHA input file'
      end

      subroutine print_MSSM_par(ifl)
      implicit double precision (a-h,o-z)
      double complex hmu
      double complex gm2,gm3
      double complex lms,rms,ums,dms,qms
      double complex ls,ks,ds,es,us,ws
      common/soft/ls(3,3),ks(3,3),ds(3,3),es(3,3),us(3,3),ws(3,3)
      common/msoft/lms(3,3),rms(3,3),ums(3,3),dms(3,3),qms(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hpar/hm1,hm2,hm12,hmu
      common/gmass/gm1,gm2,gm3
      write(ifl,*)
      write(ifl,99)'******* MSSM Lagrangian parameters *******'
      write(ifl,*)
      write(ifl,99)'QED coupling 1/alpha_em(M_Z) = ',1/alpha
      write(ifl,99)'Weinberg angle s_W^2 = ',st2
      write(ifl,99)'Z boson mass = ',zm
      write(ifl,99)'W boson mass = ',wm
      write(ifl,99)'QCD coupling alpha_s(M_Z) = ',alfas(zm)
      write(ifl,*)
      write(ifl,99)'Higgs mixing parameter mu (complex)   = ',hmu
      write(ifl,99)'Higgs soft mixing parameter m_{12}^2  = ',hm12
      write(ifl,99)'Higgs soft masses m_{H_1}^2,m_{H_2}^2 = ',hm1,hm2
      write(ifl,*)
      write(ifl,99)'U(1)  gaugino mass (complex) =',gm3
      write(ifl,99)'SU(2) gaugino mass (complex) =',gm2
      write(ifl,99)'SU(3) gaugino mass (real)    =',gm1
      write(ifl,*)
      write(ifl,99)'Left slepton mass matrix, real part:'
      call cr_mat_print(lms,3,ifl)
      write(ifl,*)
      write(ifl,99)'Left slepton mass matrix, imaginary part:'
      call ci_mat_print(lms,3,ifl)
      write(ifl,*)
      write(ifl,99)'Right slepton mass matrix, real part:'
      call cr_mat_print(rms,3,ifl)
      write(ifl,*)
      write(ifl,99)'Right slepton mass matrix, imaginary part:'
      call ci_mat_print(rms,3,ifl)
      write(ifl,*)
      write(ifl,99)'Slepton LR mixing matrix, real part:'
      call cr_mat_print(ls,3,ifl)
      write(ifl,*)
      write(ifl,99)'Slepton LR mixing matrix, imaginary part:'
      call ci_mat_print(ls,3,ifl)
      write(ifl,*)
      write(ifl,99)'Left squark mass matrix, real part:'
      call cr_mat_print(qms,3,ifl)
      write(ifl,*)
      write(ifl,99)'Left squark mass matrix, imaginary part:'
      call ci_mat_print(qms,3,ifl)
      write(ifl,*)
      write(ifl,99)'Right up-squark mass matrix, real part:'
      call cr_mat_print(ums,3,ifl)
      write(ifl,*)
      write(ifl,99)'Right up-squark mass matrix, imaginary part:'
      call ci_mat_print(ums,3,ifl)
      write(ifl,*)
      write(ifl,99)'Right down-squark mass matrix, real part:'
      call cr_mat_print(dms,3,ifl)
      write(ifl,*)
      write(ifl,99)'Right down-squark mass matrix, imaginary part:'
      call ci_mat_print(dms,3,ifl)
      write(ifl,*)
      write(ifl,99)'Up-squark LR mixing matrix, real part:'
      call cr_mat_print(us,3,ifl)
      write(ifl,*)
      write(ifl,99)'Up-squark LR mixing matrix, imaginary part:'
      call ci_mat_print(us,3,ifl)
      write(ifl,*)
      write(ifl,99)'Down-squark LR mixing matrix, real part:'
      call cr_mat_print(ds,3,ifl)
      write(ifl,*)
      write(ifl,99)'Down-squark LR mixing matrix, imaginary part:'
      call ci_mat_print(ds,3,ifl)
 99   format(a,10(1pe11.4,1x))
      return
      end

      subroutine print_MSSM_masses(ifl)
      implicit double precision (a-h,o-z)
      double complex gm2,gm3
      double complex zu,zd,zn,zpos,zneg,zv,zl
      common/hmass/cm(2),rm(2),ppm(2),zr(2,2),zh(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/neut/fnm(4),zn(4,4)
      common/fmass/em(3),um(3),dm(3)
      common/gmass/gm1,gm2,gm3
      common/hmass_EPA/pm_epa,hm1,hm2,sa,ca,sb,cb
      write(ifl,*)
      write(ifl,99)'******* Particle masses in GeV: *******'
      write(ifl,*)
      write(ifl,99)'** Fermion masses **'
      write(ifl,99)'Charged lepton masses                ',
     $     (em(ij),ij=1,3)
      write(ifl,99)'Running u quark masses at m_t scale ',
     $     (um(ij),ij=1,3)
      write(ifl,99)'Running d quark masses at m_t scale ',
     $     (dm(ij),ij=1,3)
      write(ifl,*)
      write(ifl,99)'** Higgs masses **'
      write(ifl,99)'Tree level (H,h,A,H+):          ',rm(1),rm(2),ppm(1)
     $     ,cm(1)
      write(ifl,99)'1-loop, EPA approximation (H,h):',hm1,hm2
      write(ifl,*)
      write(ifl,99)'** Tree level SUSY masses **'
      write(ifl,99)'Sneutrino masses   ',(vm(ij),ij=1,3)
      write(ifl,99)'Slepton masses     ',(slm(ij),ij=1,6)
      write(ifl,99)'U squark masses    ',(sum(ij),ij=1,6)
      write(ifl,99)'D squark masses    ',(sdm(ij),ij=1,6)
      write(ifl,99)'Chargino masses    ',(fcm(ij),ij=1,2)
      write(ifl,99)'Neutralino masses  ',(fnm(ij),ij=1,4)
      write(ifl,99)'Gluino mass        ',gm1
 99   format(a,10(1pe10.3,1x))
      return
      end

      subroutine cr_mat_print(a,n,ifl)
c     print real part of complex matrix
      double complex a(n,n)
      do i=1,n
         write(ifl,'(100(1pe12.5,1x))')(dble(a(i,j)),j=1,n)
      end do
      return
      end

      subroutine ci_mat_print(a,n,ifl)
c     print imaginary part of complex matrix
      double complex a(n,n)
      do i=1,n
         write(ifl,'(100(1pe12.5,1x))')(dimag(a(i,j)),j=1,n)
      end do
      return
      end
