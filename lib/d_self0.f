c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM} 
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: D_SELF0.F
c     Released: 28: 1:2000 (P.Ch.,J.R.)

c     Vector, scalar and pseudoscalar d quark self energy at p^2=0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Vector self-energy (proportional to k(mu)gamma(mu) = G(k)))       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dv_sig0_1(i,j)
c     Up quark + W in loop
      implicit double precision (a-h,o-z)
      double complex b1,ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      dv_sig0_1 = 0
      do l=1,3
         dv_sig0_1 = dv_sig0_1 + e2/2/st2*dconjg(ckm(i,l))*ckm(j,l)
     $        * b1(0.d0,um(l),wm)
      end do
      return
      end
      
      double complex function dv_sig0_4(i,j)
c     Up quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b1,ckm
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/fmass/em(3),um(3),dm(3)
      dv_sig0_4 = 0
      do l=1,3
         do k=1,2
            dv_sig0_4 = dv_sig0_4 + ((zh(2,k)*yu(l))**2
     $           + yd(i)*yd(j)*zh(1,k)**2)/2
     $           * dconjg(ckm(i,l))*ckm(j,l)*b1(0.d0,um(l),cm(k))
         end do
      end do
      return
      end

      double complex function dv_sig0_7(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zn
      double complex vl_ddn,vr_ddn
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      dv_sig0_7 = 0
      do k=1,6
        do l=1,4
          dv_sig0_7 = dv_sig0_7 
     $        + (vl_ddn(i,k,l)*dconjg(vl_ddn(j,k,l))
     $        + vr_ddn(i,k,l)*dconjg(vr_ddn(j,k,l)))/2
     $        * b1(0.d0,fnm(l),sdm(k))
        end do
      end do
      return
      end

      double complex function dv_sig0_8(i,j)
c     Chargino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zpos,zneg
      double complex vl_duc,vr_duc
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      dv_sig0_8 = 0
      do k=1,6
        do l=1,2
          dv_sig0_8 = dv_sig0_8 
     $        + (vl_duc(i,k,l)*dconjg(vl_duc(j,k,l))
     $        + vr_duc(i,k,l)*dconjg(vr_duc(j,k,l)))/2
     $        * b1(0.d0,fcm(l),sum(k))
        end do
      end do
      return
      end

      double complex function dv_sig0_9(i,j)
c     Gluino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      dv_sig0_9 = 0
      if (init_alpha_susy) call init_alpha_s_susy
      al = 4*g3d*g3d/3.d0
      do k=1,6
        dv_sig0_9 = dv_sig0_9 + al*(zd(i,k)*dconjg(zd(j,k))
     $      + zd(i+3,k)*dconjg(zd(j+3,k)))*b1(0.d0,gm1,sdm(k))
      end do
      return
      end

      double complex function dv_sig0(i,j)
c      Full bare down quark self-energy, vector part
      implicit double precision (a-h,o-z)
      double complex dv_sig0_1,dv_sig0_4,dv_sig0_7,dv_sig0_8,dv_sig0_9
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      common/d_self_stat/iv,ia,is,ip
      external init_d_self
      dv_sig0 = (0.d0,0.d0)
      if (ih.eq.1) dv_sig0 = dv_sig0 + dv_sig0_1(i,j) + dv_sig0_4(i,j)
      if (in.eq.1) dv_sig0 = dv_sig0 + dv_sig0_7(i,j)
      if (ic.eq.1) dv_sig0 = dv_sig0 + dv_sig0_8(i,j)
      if (ig.eq.1) dv_sig0 = dv_sig0 + dv_sig0_9(i,j)
      dv_sig0 = iv*dv_sig0/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Axial self-energy (proportional to k(mu)gamma(mu)gamma(5)         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function da_sig0_1(i,j)
c     Up quark + W in loop
      implicit double precision (a-h,o-z)
      double complex b1,ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      da_sig0_1 = 0
      do l=1,3
         da_sig0_1 = da_sig0_1 - e2/2/st2*dconjg(ckm(i,l))*ckm(j,l)
     $        * b1(0.d0,um(l),wm)
      end do
      return
      end
      
      double complex function da_sig0_4(i,j)
c     Up quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b1,ckm
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      da_sig0_4 = 0
      do l=1,3
         do k=1,2
            da_sig0_4 = da_sig0_4 - ((zh(2,k)*yu(l))**2
     $           - yd(i)*yd(j)*zh(1,k)**2)/2
     $           * dconjg(ckm(i,l))*ckm(j,l)*b1(0.d0,um(l),cm(k))
         end do
      end do
      return
      end

      double complex function da_sig0_7(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zn
      double complex vl_ddn,vr_ddn
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      da_sig0_7 = 0
      do k=1,6
         do l=1,4
            da_sig0_7 = da_sig0_7 
     $           - (vl_ddn(i,k,l)*dconjg(vl_ddn(j,k,l))
     $           - vr_ddn(i,k,l)*dconjg(vr_ddn(j,k,l)))/2
     $           * b1(0.d0,fnm(l),sdm(k))
         end do
      end do
      return
      end

      double complex function da_sig0_8(i,j)
c     Chargino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zpos,zneg
      double complex vl_duc,vr_duc
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      da_sig0_8 = 0
      do k=1,6
         do l=1,2
            da_sig0_8 = da_sig0_8 
     $           - (vl_duc(i,k,l)*dconjg(vl_duc(j,k,l))
     $           - vr_duc(i,k,l)*dconjg(vr_duc(j,k,l)))/2
     $           * b1(0.d0,fcm(l),sum(k))
         end do
      end do
      return
      end

      double complex function da_sig0_9(i,j)
c     Gluino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      da_sig0_9 = 0
      if (init_alpha_susy) call init_alpha_s_susy
      al = 4*g3d*g3d/3.d0
      do  k=1,6
         da_sig0_9 = da_sig0_9 - al*(zd(i,k)*dconjg(zd(j,k))
     $        - zd(i+3,k)*dconjg(zd(j+3,k)))*b1(0.d0,gm1,sdm(k))
      end do
      return
      end

      double complex function da_sig0(i,j)
c     Full bare down quark self-energy, axial part
      implicit double precision (a-h,o-z)
      double complex da_sig0_1,da_sig0_4,da_sig0_7,da_sig0_8,da_sig0_9
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      common/d_self_stat/iv,ia,is,ip
      external init_d_self
      da_sig0 = (0.d0,0.d0)
      if (ih.eq.1) da_sig0 = da_sig0 + da_sig0_1(i,j) + da_sig0_4(i,j) 
      if (in.eq.1) da_sig0 = da_sig0 + da_sig0_7(i,j) 
      if (ic.eq.1) da_sig0 = da_sig0 + da_sig0_8(i,j) 
      if (ig.eq.1) da_sig0 = da_sig0 + da_sig0_9(i,j) 
      da_sig0 = ia*da_sig0/16/pi/pi
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Scalar self-energy                                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function ds_sig0_4(i,j)
c     Up quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b0,ckm
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/fmass/em(3),um(3),dm(3)
      ds_sig0_4 = 0
      do l=1,3
         do k=1,2
            ds_sig0_4 = ds_sig0_4 + zh(2,k)*zh(1,k)*um(l)*yu(l)/2
     $           * (yd(i) + yd(j))
     $           * dconjg(ckm(i,l))*ckm(j,l)*b0(0.d0,um(l),cm(k))
         end do
      end do
      return
      end

      double complex function ds_sig0_7(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zn
      double complex vl_ddn,vr_ddn
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      ds_sig0_7 = 0
      do k=1,6
        do l=1,4
          ds_sig0_7 = ds_sig0_7 
     $        - (vl_ddn(i,k,l)*dconjg(vr_ddn(j,k,l))
     $        + vr_ddn(i,k,l)*dconjg(vl_ddn(j,k,l)))*fnm(l)/2
     $        * b0(0.d0,fnm(l),sdm(k))
        end do
      end do
      return
      end

      double complex function ds_sig0_8(i,j)
c     Chargino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zpos,zneg
      double complex vl_duc,vr_duc
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      ds_sig0_8 = 0
      do k=1,6
        do l=1,2
          ds_sig0_8 = ds_sig0_8 
     $        - (vl_duc(i,k,l)*dconjg(vr_duc(j,k,l))
     $        + vr_duc(i,k,l)*dconjg(vl_duc(j,k,l)))*fcm(l)/2
     $        * b0(0.d0,fcm(l),sum(k))
        end do
      end do
      return
      end

      double complex function ds_sig0_9(i,j)
c     Gluino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      ds_sig0_9 = 0
      al = 8*g3d*g3d/3.d0
      do  k=1,6
         ds_sig0_9 = ds_sig0_9 + al*(zd(i,k)*dconjg(zd(j+3,k))
     $        + zd(i+3,k)*dconjg(zd(j,k)))*gm1*b0(0.d0,gm1,sdm(k))/2
      end do
      return
      end

      double complex function ds_sig0(i,j)
      implicit double precision (a-h,o-z)
      double complex ds_sig0_4,ds_sig0_7,ds_sig0_8,ds_sig0_9
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      common/d_self_stat/iv,ia,is,ip
      external init_d_self
      ds_sig0 = (0.d0,0.d0)
      if (ih.eq.1) ds_sig0 = ds_sig0 + ds_sig0_4(i,j) 
      if (in.eq.1) ds_sig0 = ds_sig0 + ds_sig0_7(i,j) 
      if (ic.eq.1) ds_sig0 = ds_sig0 + ds_sig0_8(i,j) 
      if (ig.eq.1) ds_sig0 = ds_sig0 + ds_sig0_9(i,j) 
      ds_sig0 = is*ds_sig0/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Pseudoscalar self-energy (proportional to G(5))                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dp_sig0_4(i,j)
c     Up quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b0,ckm
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/km_mat/ckm(3,3)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      dp_sig0_4 = 0
      do l=1,3
         do k=1,2
            dp_sig0_4 = dp_sig0_4 + zh(2,k)*zh(1,k)*um(l)*yu(l)
     $           * (yd(i) - yd(j))/2
     $           * dconjg(ckm(i,l))*ckm(j,l)*b0(0.d0,um(l),cm(k))
         end do
      end do
      return
      end

      double complex function dp_sig0_7(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zn
      double complex vl_ddn,vr_ddn
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      dp_sig0_7 = 0
      do k=1,6
        do l=1,4
          dp_sig0_7 = dp_sig0_7 
     $        + (vl_ddn(i,k,l)*dconjg(vr_ddn(j,k,l))
     $        - vr_ddn(i,k,l)*dconjg(vl_ddn(j,k,l)))*fnm(l)/2
     $        * b0(0.d0,fnm(l),sdm(k))
        end do
      end do
      return
      end

      double complex function dp_sig0_8(i,j)
c     Chargino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zpos,zneg
      double complex vl_duc,vr_duc
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      dp_sig0_8 = 0
      do k=1,6
        do l=1,2
          dp_sig0_8 = dp_sig0_8 
     $        + (vl_duc(i,k,l)*dconjg(vr_duc(j,k,l))
     $        - vr_duc(i,k,l)*dconjg(vl_duc(j,k,l)))*fcm(l)/2
     $        * b0(0.d0,fcm(l),sum(k))
        end do
      end do
      return
      end

      double complex function dp_sig0_9(i,j)
c     Gluino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      dp_sig0_9 = 0
      al = 8*g3d*g3d/3.d0
      do k=1,6
        dp_sig0_9 = dp_sig0_9 - al*(zd(i,k)*dconjg(zd(j+3,k))
     $      - zd(i+3,k)*dconjg(zd(j,k)))*gm1*b0(0.d0,gm1,sdm(k))/2
      end do
      return
      end
      
      double complex function dp_sig0(i,j)
      implicit double precision (a-h,o-z)
      double complex dp_sig0_4,dp_sig0_7,dp_sig0_8,dp_sig0_9
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      common/d_self_stat/iv,ia,is,ip
      external init_d_self
      dp_sig0 = (0.d0,0.d0)
      if (ih.eq.1) dp_sig0 = dp_sig0 + dp_sig0_4(i,j) 
      if (in.eq.1) dp_sig0 = dp_sig0 + dp_sig0_7(i,j) 
      if (ic.eq.1) dp_sig0 = dp_sig0 + dp_sig0_8(i,j) 
      if (ig.eq.1) dp_sig0 = dp_sig0 + dp_sig0_9(i,j) 
      dp_sig0 = ip*dp_sig0/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Basis change: 1,gamma(5) -> P_L,P_R                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dvl_sig0(i,j)
      implicit double precision (a-h,o-z)
      double complex dv_sig0,da_sig0
      dvl_sig0 = dv_sig0(i,j) - da_sig0(i,j)
      return
      end

      double complex function dvr_sig0(i,j)
      implicit double precision (a-h,o-z)
      double complex dv_sig0,da_sig0
      dvr_sig0 = dv_sig0(i,j) + da_sig0(i,j)
      return
      end

      double complex function dsl_sig0(i,j)
      implicit double precision (a-h,o-z)
      double complex ds_sig0,dp_sig0
      dsl_sig0 = ds_sig0(i,j) - dp_sig0(i,j)
      return
      end

      double complex function dsr_sig0(i,j)
      implicit double precision (a-h,o-z)
      double complex ds_sig0,dp_sig0
      dsr_sig0 = ds_sig0(i,j) + dp_sig0(i,j)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccc

      block data init_d_self
      implicit double precision (a-h,o-z)
      common/d_self_stat/iv,ia,is,ip
      data iv,ia,is,ip/4*1/
      end



