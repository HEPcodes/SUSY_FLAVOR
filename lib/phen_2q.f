c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: PHEN_2Q.F
c     Released: 29:08:2007 (J.R.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for QCD evolution of          c
c     Wilson coefficients of the operators present in 2-d quark    c
c     mixing and for the phenomenological quantities like B->ll,   c
c     B->Xll, K/B-> pivv (to be added gradually).                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine eta_2q_evol
c     QCD evolution of Wilson coefficients of the effective
c     2-quark operators
      implicit double precision (a-h,o-z)
      logical init_eta_2q,init_alpha_susy
      dimension qm(3)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
c     put into the common below evolution coefficients...
      common/ev_mat_2q/vx(2),sx(2),tx(2),init_eta_2q
      common/fmass_high/umu(3),uml(3),amuu(3),dmu(3),dml(3),amud(3)
      common/bx_4q/bk(5),bd(5),bb(2,5),amu_k,amu_d,amu_b
c     running s,b masses at m_t scale
      call init_run_qmass
c      NLO alpha_s at M_Z
      al = alfas_nlo(zm)
c      NLO alpha_s at M_SUSY scale(s)
c      M_S_{U,D} = (M_gluino + average M_{U,D} mass)/2
      if (init_alpha_susy) call init_alpha_s_susy
c      running quark masses (c,b,t)
      qm(1) = uml(2)
      qm(2) = dml(3)
      qm(3) = umu(3)
c      evolution from M_SUSY to mu_t
c      start from M_SUSY = (M_gluino + averaged M_D)/2
      al_s = alfas_nlo(qm(3))/4/pi
      eta = g3d*g3d/al_s/16/pi/pi
c     define evolution below for V, S and T sectors:
      vx(2) = 1
      sx(2) = 1
      tx(2) = 1
c      evolution from mu_t to mu_B
      al_s = alfas_nlo(amu_b)/4/pi
      eta = alfas_nlo(qm(3))/alfas_nlo(amu_b)
c     define evolution below for V, S and T sectors:
      vx(1) = 1
      sx(1) = 1
      tx(1) = 1
c     initialization of QCD factors finished:
      init_eta_2q = .false.
      return
      end

      subroutine dl_wil_run(i,j,k,l)
      implicit double precision (a-h,o-z)
      logical init_eta_2q
      double complex dl_vll,dl_vrr,dl_vlr,dl_vrl,dl_sll,dl_srr,
     $     dl_slr,dl_srl,dl_tl,dl_tr
      double complex dls_vll,dls_vrr,dls_vlr,dls_vrl,dls_sll,dls_srr,
     $     dls_slr,dls_srl,dls_tl,dls_tr
      common/dl_wil_coeff/dls_vll,dls_vrr,dls_vlr,dls_vrl,
     $     dls_sll,dls_srr,dls_slr,dls_srl,dls_tl,dls_tr
      common/ev_mat_2q/vx(2),sx(2),tx(2),init_eta_2q
      if (init_eta_2q) call eta_2q_evol()
c     NLO evolution of the MSSM part
c     First step: from mu = M_SUSY = (M_gluino + aver M_D)/2 to mu = m_t
      dls_vll  = vx(2)*dl_vll(i,j,k,l)
      dls_vrr  = vx(2)*dl_vrr(i,j,k,l)
      dls_vlr  = vx(2)*dl_vlr(i,j,k,l)
      dls_vrl  = vx(2)*dl_vrl(i,j,k,l)
      dls_sll  = sx(2)*dl_sll(i,j,k,l)
      dls_srr  = sx(2)*dl_srr(i,j,k,l)
      dls_slr  = sx(2)*dl_slr(i,j,k,l)
      dls_srl  = sx(2)*dl_srl(i,j,k,l)
      dls_tl   = tx(2)*dl_tl(i,j,k,l)
      dls_tr   = tx(2)*dl_tr(i,j,k,l)
c     Second step: from mu = m_t to mu = m_b
      dls_vll  = vx(1)*dls_vll
      dls_vrr  = vx(1)*dls_vrr
      dls_vlr  = vx(1)*dls_vlr
      dls_vrl  = vx(1)*dls_vrl
      dls_sll  = sx(1)*dls_sll
      dls_srr  = sx(1)*dls_srr
      dls_slr  = sx(1)*dls_slr
      dls_srl  = sx(1)*dls_srl
      dls_tl   = tx(1)*dls_tl
      dls_tr   = tx(1)*dls_tr
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Br(B->l^+l^-) calculation                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function b_ll(i,j,k,l)
c     Decay of B_d(s)-> l_K l_L pair
c     i=3 denotes B decay, j=3 is \bar B decay
c     second index j,i=1 or 2 defines B meson, B_d or B_s respectively
      implicit double precision (a-h,o-z)
      double complex fs,fp,fv,fa
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/meson_data/dmk,amk,epsk,fk,dmd,amd,fd,
     $    amb(2),dmb(2),gam_b(2),fb(2)
      common/fmass/em(3),um(3),dm(3)
      double complex dls_vll,dls_vrr,dls_vlr,dls_vrl,dls_sll,dls_srr,
     $     dls_slr,dls_srl,dls_tl,dls_tr
      common/dl_wil_coeff/dls_vll,dls_vrr,dls_vlr,dls_vrl,
     $     dls_sll,dls_srr,dls_slr,dls_srl,dls_tl,dls_tr
      common/ph_units/hbar,gev
      external init_4q,init_2q,init_units
      if (max(i,j).ne.3) stop 'No b quark index in b_ll?'
      call dl_wil_run(i,j,k,l)
c     form factors
      ii = min(i,j)      
      amb2 = amb(ii)*amb(ii)
      del = em(k) - em(l)
      sum = em(k) + em(l)
      fv = dls_vrr - dls_vll + dls_vrl - dls_vlr 
      fa = dls_vll + dls_vrr - dls_vlr - dls_vrl
      fs = amb2/(dm(3) + dm(ii))*(dls_sll - dls_srr + dls_slr - dls_srl)
      fp = amb2/(dm(3) + dm(ii))*(dls_slr + dls_srl - dls_sll - dls_srr)
c     matrix element
      b_ll = (amb2 - sum*sum)*abs(fs)**2 + (amb2 - del*del)*abs(fp)**2
     $     + sum*sum*(amb2 - del*del)*abs(fa)**2
     $     + 2*sum*(amb2 - del*del)*dble(fp*dconjg(fa))

      if (k.ne.l) b_ll = b_ll + del*del*(amb2 - sum*sum)*abs(fv)**2
     $     - 2*del*(amb2 - sum*sum)*dble(fs*dconjg(fv))
c     branching ratio itself
      b_ll = b_ll*fb(ii)*fb(ii)*gam_b(ii)/128/pi/amb(ii)/hbar
     $     * sqrt(1 - sum*sum/amb2)*sqrt(1 - del*del/amb2)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Br(K_L -> pi^0 \bar v v) and  Br(K^+ -> pi^+ \bar v v) calculation  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine k_pivv(br_k0,br_kp)
c     Decays K^0_L -> pi^0 \bar v v and K^+ -> pi^+ \bar v v
c     compare hep-ph/0408142
      implicit double precision (a-h,o-z)
      double complex xx
      double complex dd_vv_l,dd_vv_r
      double complex ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/kpivv/ak0,del_ak0,akp,del_akp,pc,del_pc,alam
      common/km_mat/ckm(3,3)
      external init_2q
      xx = (0.d0,0.d0)
      do k=1,3
         do l=1,3
            xx = xx + dd_vv_l(2,1,k,l) + dd_vv_r(2,1,k,l)
         end do 
      end do
      xx = - (4*pi*st2*wm/e2)**2*xx/3/alam**5
c     Br(K^0_L -> pi^0 \bar v v)
      br_k0 = ak0*dimag(xx)**2
c     Br(K^+ -> pi^+ \bar v v)
      br_kp = akp*(dimag(xx)**2 
     $     + dble(dconjg(ckm(2,2))*ckm(1,2)*pc/alam + xx)**2)
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Data initialization block                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      block data init_2q
      implicit double precision (a-h,o-z)
      logical init_eta_2q
      common/ev_mat_2q/vx(2),sx(2),tx(2),init_eta_2q
      common/kpivv/ak0,del_ak0,akp,del_akp,pc,del_pc,alam
      data init_eta_2q/.true./
      data ak0,akp/2.231d-10,5.173d-11/
      data del_ak0,del_akp/0.013d-10,0.024d-11/
      data pc,del_pc/0.41d0,0.03d0/
      data alam/0.225d0/
      end 
