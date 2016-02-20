c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: DD_LL.F
c     Released: 02:08:2007 (J.R.)
c     Revised: 23.10.2008
c     Higgs penguin diagrams addded

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expression for the coefficient of the    c
c     effective Hamiltonian for the B/K->ll and B/K->Xll decays   c
c     (only box diagrams at the quark level).                     c
c     General form of the Hamiltonian is:                         c
c     H = A^V_LL H^V_LL + A^V_RR H^V_RR                           c
c       + A^V_LR H^V_LR + A^V_LR H^V_RL                           c
c       + A^S_LL H^S_LL + A^S_RR H^S_RR                           c
c       + A^S_LR H^S_LR + A^S_LR H^S_RL                           c
c       + A^T_L H^T_L   + A^T_R H^T_R                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_LL                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_vll_g(i,j,k,l)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      double complex ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      dl_vll_g = (0.d0,0.d0)
      if (k.ne.l) return
      do m=1,3
         dl_vll_g = dl_vll_g + dconjg(ckm(i,m))*ckm(j,m)
     $        * dp1(wm,wm,um(m),0.d0)
      end do
      dl_vll_g = e2*e2/4/st2/st2*dl_vll_g
      return
      end

      double complex function dl_vll_h(i,j,k,l)
c     Double Higgs contribution
      implicit double precision (a-h,o-z)
      dl_vll_h = (0.d0,0.d0)
      return
      end

      double complex function dl_vll_hg(i,j,k,l)
c     W-Higgs contribution
      implicit double precision (a-h,o-z)
      dl_vll_hg = (0.d0,0.d0)
      return
      end

      double complex function dl_vll_c(i,j,k,l)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vl_duc,zpos,zneg,zu,zd,zv,zl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_vll_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do nn=1,3
               do mm=1,6
                 dl_vll_c = dl_vll_c + vl_duc(i,mm,m)*zpos(1,n)*zv(l,nn)
     $                 * dconjg(vl_duc(j,mm,n)*zpos(1,m)*zv(k,nn))
     $                 * dp1(fcm(m),fcm(n),sum(mm),vm(nn))
               end do
            end do
         end do
      end do
      dl_vll_c = e2/st2/4*dl_vll_c
      return
      end

      double complex function dl_vll_n(i,j,k,l)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vl_lln,zn,zu,zd,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_vll_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do mm=1,6
               do nn=1,6
                  dl_vll_n = dl_vll_n + vl_ddn(i,mm,m)*vl_lln(k,nn,n)
     $                 * dconjg(vl_ddn(j,mm,n)*vl_lln(l,nn,m))/4
     $                 * dp1(fnm(m),fnm(n),sdm(mm),slm(nn))
     $                 + vl_ddn(i,mm,m)*vl_lln(k,nn,m)
     $                 * dconjg(vl_ddn(j,mm,n)*vl_lln(l,nn,n))/2
     $                 *fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(mm),slm(nn))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dl_vll_box(i,j,k,l)
c     Full box A^V_LL formfactor
      implicit double precision (a-h,o-z)
      double complex dl_vll_g,dl_vll_h,dl_vll_hg,dl_vll_c,dl_vll_n
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      dl_vll_box = (ih*dl_vll_g(i,j,k,l) + ih*dl_vll_h(i,j,k,l)
     $     + ic*dl_vll_c(i,j,k,l) + in*dl_vll_n(i,j,k,l)
     $     + ih*dl_vll_hg(i,j,k,l))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_RR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_vrr_g(i,j,k,l)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      dl_vrr_g = (0.d0,0.d0)
      return
      end

      double complex function dl_vrr_h(i,j,k,l)
c     Double Higgs contribution
      implicit double precision (a-h,o-z)
      double complex yh_eff_r
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_vrr_h = (0.d0,0.d0)
      if (k.ne.l) return
      do m=1,3
         do kk=1,2
            do ll=1,2
               dl_vrr_h = dl_vrr_h + zh(1,kk)*zh(1,ll)*yh_eff_r(i,m,kk)
     $              * dconjg(yh_eff_r(j,m,ll))
     $              * dp0(cm(kk),cm(ll),um(m),0.d0)
            end do
         end do
      end do
      dl_vrr_h = yl(k)*yl(l)/4*dl_vrr_h
      return
      end

      double complex function dl_vrr_hg(i,j,k,l)
c     W-Higgs contribution
      implicit double precision (a-h,o-z)
      dl_vrr_hg = (0.d0,0.d0)
      return
      end

      double complex function dl_vrr_c(i,j,k,l)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vr_duc,zpos,zneg,zu,zd,zv,zl
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_vrr_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do nn=1,3
               do mm=1,6
                 dl_vrr_c = dl_vrr_c + vr_duc(i,mm,m)*zneg(2,m)*zv(l,nn)
     $                 * dconjg(vr_duc(j,mm,n)*zneg(2,n)*zv(k,nn))
     $                 * dp1(fcm(m),fcm(n),sum(mm),vm(nn))
               end do
            end do
         end do
      end do
      dl_vrr_c = yl(k)*yl(l)/4*dl_vrr_c
      return
      end

      double complex function dl_vrr_n(i,j,k,l)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vr_ddn,vr_lln,zn,zu,zd,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_vrr_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do mm=1,6
               do nn=1,6
                  dl_vrr_n = dl_vrr_n + vr_ddn(i,mm,m)*vr_lln(k,nn,n)
     $                 * dconjg(vr_ddn(j,mm,n)*vr_lln(l,nn,m))/4
     $                 * dp1(fnm(m),fnm(n),sdm(mm),slm(nn))
     $                 + vr_ddn(i,mm,m)*vr_lln(k,nn,m)
     $                 * dconjg(vr_ddn(j,mm,n)*vr_lln(l,nn,n))/2
     $                 *fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(mm),slm(nn))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dl_vrr_box(i,j,k,l)
c     Full box A^V_RR formfactor
      implicit double precision (a-h,o-z)
      double complex dl_vrr_g,dl_vrr_h,dl_vrr_hg,dl_vrr_c,dl_vrr_n
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      dl_vrr_box = (ih*dl_vrr_g(i,j,k,l) + ih*dl_vrr_h(i,j,k,l)
     $     + ic*dl_vrr_c(i,j,k,l) + in*dl_vrr_n(i,j,k,l)
     $     + ih*dl_vrr_hg(i,j,k,l))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_LR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_vlr_g(i,j,k,l)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      dl_vlr_g = (0.d0,0.d0)
      return
      end

      double complex function dl_vlr_h(i,j,k,l)
c     Double Higgs contribution
      implicit double precision (a-h,o-z)
      double complex yh_eff_l
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_vlr_h = (0.d0,0.d0)
      if (k.ne.l) return
      do m=1,3
         do kk=1,2
            do ll=1,2
               dl_vlr_h = dl_vlr_h + zh(1,kk)*zh(1,ll)*yh_eff_l(i,m,kk)
     $              * dconjg(yh_eff_l(j,m,ll))
     $              * dp0(cm(kk),cm(ll),um(m),0.d0)
            end do
         end do
      end do
      dl_vlr_h = yl(k)*yl(l)/4*dl_vlr_h
      return
      end

      double complex function dl_vlr_hg(i,j,k,l)
c     W-Higgs contribution
      implicit double precision (a-h,o-z)
      dl_vlr_hg = (0.d0,0.d0)
      return
      end

      double complex function dl_vlr_c(i,j,k,l)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vl_duc,zpos,zneg,zu,zd,zv,zl
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_vlr_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do nn=1,3
               do mm=1,6
                 dl_vlr_c = dl_vlr_c + vl_duc(i,mm,m)*zneg(2,m)*zv(l,nn)
     $                 * dconjg(vl_duc(j,mm,n)*zneg(2,n)*zv(k,nn))
     $                 *fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(mm),vm(nn))
               end do
            end do
         end do
      end do
      dl_vlr_c = - yl(k)*yl(l)/2*dl_vlr_c
      return
      end

      double complex function dl_vlr_n(i,j,k,l)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_lln,zn,zu,zd,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_vlr_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do mm=1,6
               do nn=1,6
                  dl_vlr_n = dl_vlr_n - vl_ddn(i,mm,m)*vr_lln(k,nn,n)
     $                 * dconjg(vl_ddn(j,mm,n)*vr_lln(l,nn,m))/2
     $                 *fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(mm),slm(nn))
     $                 - vl_ddn(i,mm,m)*vr_lln(k,nn,m)
     $                 * dconjg(vl_ddn(j,mm,n)*vr_lln(l,nn,n))/4
     $                 * dp1(fnm(m),fnm(n),sdm(mm),slm(nn))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dl_vlr_box(i,j,k,l)
c     Full box A^V_LR formfactor
      implicit double precision (a-h,o-z)
      double complex dl_vlr_g,dl_vlr_h,dl_vlr_hg,dl_vlr_c,dl_vlr_n
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      dl_vlr_box = (ih*dl_vlr_g(i,j,k,l) + ih*dl_vlr_h(i,j,k,l)
     $     + ic*dl_vlr_c(i,j,k,l) + in*dl_vlr_n(i,j,k,l)
     $     + ih*dl_vlr_hg(i,j,k,l))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_RL                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_vrl_g(i,j,k,l)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      dl_vrl_g = (0.d0,0.d0)
      return
      end

      double complex function dl_vrl_h(i,j,k,l)
c     Double Higgs contribution
      implicit double precision (a-h,o-z)
      dl_vrl_h = (0.d0,0.d0)
      return
      end

      double complex function dl_vrl_hg(i,j,k,l)
c     W-Higgs contribution
      implicit double precision (a-h,o-z)
      dl_vrl_hg = (0.d0,0.d0)
      return
      end

      double complex function dl_vrl_c(i,j,k,l)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vr_duc,zpos,zneg,zu,zd,zv,zl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_vrl_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do nn=1,3
               do mm=1,6
                 dl_vrl_c = dl_vrl_c + vr_duc(i,mm,m)*zpos(1,n)*zv(l,nn)
     $                 * dconjg(vr_duc(j,mm,n)*zpos(1,m)*zv(k,nn))
     $                 *fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(mm),vm(nn))
               end do
            end do
         end do
      end do
      dl_vrl_c = - e2/st2/2*dl_vrl_c
      return
      end

      double complex function dl_vrl_n(i,j,k,l)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vr_ddn,vl_lln,zn,zu,zd,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_vrl_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do mm=1,6
               do nn=1,6
                  dl_vrl_n = dl_vrl_n - vr_ddn(i,mm,m)*vl_lln(k,nn,n)
     $                 * dconjg(vr_ddn(j,mm,n)*vl_lln(l,nn,m))/2
     $                 *fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(mm),slm(nn))
     $                 - vr_ddn(i,mm,m)*vl_lln(k,nn,m)
     $                 * dconjg(vr_ddn(j,mm,n)*vl_lln(l,nn,n))/4
     $                 * dp1(fnm(m),fnm(n),sdm(mm),slm(nn))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dl_vrl_box(i,j,k,l)
c     Full box A^V_RL formfactor
      implicit double precision (a-h,o-z)
      double complex dl_vrl_g,dl_vrl_h,dl_vrl_hg,dl_vrl_c,dl_vrl_n
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      dl_vrl_box = (ih*dl_vrl_g(i,j,k,l) + ih*dl_vrl_h(i,j,k,l)
     $     + ic*dl_vrl_c(i,j,k,l) + in*dl_vrl_n(i,j,k,l)
     $     + ih*dl_vrl_hg(i,j,k,l))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^S_LL                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_sll_g(i,j,k,l)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      dl_sll_g = (0.d0,0.d0)
      return
      end

      double complex function dl_sll_h(i,j,k,l)
c     Double Higgs contribution
      implicit double precision (a-h,o-z)
      dl_sll_h = (0.d0,0.d0)
      return
      end

      double complex function dl_sll_hg(i,j,k,l)
c     W-Higgs contribution
      implicit double precision (a-h,o-z)
      dl_sll_hg = (0.d0,0.d0)
      return
      end

      double complex function dl_sll_c(i,j,k,l)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd,zv,zl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_sll_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do nn=1,3
               do mm=1,6
                 dl_sll_c = dl_sll_c + vl_duc(i,mm,m)*zpos(1,n)*zv(l,nn)
     $                 * zneg(2,m)*dconjg(vr_duc(j,mm,n)*zv(k,nn))
     $                 * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(mm),vm(nn))
               end do
            end do
         end do
      end do
      dl_sll_c = - e/st*yl(l)/2*dl_sll_c
      return
      end

      double complex function dl_sll_n(i,j,k,l)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,vl_lln,vr_lln,zn,zu,zd,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_sll_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do mm=1,6
               do nn=1,6
                  dl_sll_n = dl_sll_n
     $                 - (vl_lln(k,nn,n)*dconjg(vr_lln(l,nn,m))
     $                 + vl_lln(k,nn,m)*dconjg(vr_lln(l,nn,n)))
     $                 * vl_ddn(i,mm,m)*dconjg(vr_ddn(j,mm,n))/2
     $                 *fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(mm),slm(nn))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dl_sll_box(i,j,k,l)
c     Full box A^S_LL formfactor
      implicit double precision (a-h,o-z)
      double complex dl_sll_g,dl_sll_h,dl_sll_hg,dl_sll_c,dl_sll_n
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      dl_sll_box = (ih*dl_sll_g(i,j,k,l) + ih*dl_sll_h(i,j,k,l)
     $     + ic*dl_sll_c(i,j,k,l) + in*dl_sll_n(i,j,k,l)
     $     + ih*dl_sll_hg(i,j,k,l))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^S_RR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_srr_g(i,j,k,l)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      dl_srr_g = (0.d0,0.d0)
      return
      end

      double complex function dl_srr_h(i,j,k,l)
c     Double Higgs contribution
      implicit double precision (a-h,o-z)
      dl_srr_h = (0.d0,0.d0)
      return
      end

      double complex function dl_srr_hg(i,j,k,l)
c     W-Higgs contribution
      implicit double precision (a-h,o-z)
      dl_srr_hg = (0.d0,0.d0)
      return
      end

      double complex function dl_srr_c(i,j,k,l)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd,zv,zl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_srr_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do nn=1,3
               do mm=1,6
                 dl_srr_c = dl_srr_c + vr_duc(i,mm,m)*zv(l,nn)
     $              *dconjg(vl_duc(j,mm,n)*zpos(1,m)*zneg(2,n)*zv(k,nn))
     $              * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(mm),vm(nn))
               end do
            end do
         end do
      end do
      dl_srr_c = - e/st*yl(k)/2*dl_srr_c
      return
      end

      double complex function dl_srr_n(i,j,k,l)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,vl_lln,vr_lln,zn,zu,zd,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_srr_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do mm=1,6
               do nn=1,6
                  dl_srr_n = dl_srr_n
     $                 - (vr_lln(k,nn,n)*dconjg(vl_lln(l,nn,m))
     $                 + vr_lln(k,nn,m)*dconjg(vl_lln(l,nn,n)))
     $                 * vr_ddn(i,mm,m)*dconjg(vl_ddn(j,mm,n))/2
     $                 *fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(mm),slm(nn))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dl_srr_box(i,j,k,l)
c     Full box A^S_RR formfactor
      implicit double precision (a-h,o-z)
      double complex dl_srr_g,dl_srr_h,dl_srr_hg,dl_srr_c,dl_srr_n
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      dl_srr_box = (ih*dl_srr_g(i,j,k,l) + ih*dl_srr_h(i,j,k,l)
     $     + ic*dl_srr_c(i,j,k,l) + in*dl_srr_n(i,j,k,l)
     $     + ih*dl_srr_hg(i,j,k,l))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^S_LR                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_slr_g(i,j,k,l)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      dl_slr_g = (0.d0,0.d0)
      return
      end

      double complex function dl_slr_h(i,j,k,l)
c     Double Higgs contribution
      implicit double precision (a-h,o-z)
      dl_slr_h = (0.d0,0.d0)
      return
      end

      double complex function dl_slr_hg(i,j,k,l)
c     W-Higgs contribution
      implicit double precision (a-h,o-z)
      double complex yh_eff_r
      double complex ckm
      common/km_mat/ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_slr_hg = (0.d0,0.d0)
      if (k.ne.l) return
      do m=1,3
         do n=1,2
            dl_slr_hg = dl_slr_hg + dconjg(ckm(i,m)*yh_eff_r(j,m,n))
     $           * zh(1,n)*dp1(um(m),cm(n),wm,0.d0)
         end do
      end do
      dl_slr_hg = e2/st2/2*yl(k)*dl_slr_hg
      return
      end

      double complex function dl_slr_c(i,j,k,l)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd,zv,zl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_slr_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do nn=1,3
               do mm=1,6
                 dl_slr_c = dl_slr_c + vl_duc(i,mm,m)*zv(l,nn)
     $              *dconjg(vr_duc(j,mm,n)*zpos(1,m)*zneg(2,n)*zv(k,nn))
     $              * dp1(fcm(m),fcm(n),sum(mm),vm(nn))
               end do
            end do
         end do
      end do
      dl_slr_c = - e/st*yl(k)/2*dl_slr_c
      return
      end

      double complex function dl_slr_n(i,j,k,l)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,vl_lln,vr_lln,zn,zu,zd,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_slr_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do mm=1,6
               do nn=1,6
                  dl_slr_n = dl_slr_n
     $                 - (vr_lln(k,nn,n)*dconjg(vl_lln(l,nn,m))
     $                 + vr_lln(k,nn,m)*dconjg(vl_lln(l,nn,n)))
     $                 * vl_ddn(i,mm,m)*dconjg(vr_ddn(j,mm,n))/2
     $                 * dp1(fnm(m),fnm(n),sdm(mm),slm(nn))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dl_slr_box(i,j,k,l)
c     Full box A^S_LR formfactor
      implicit double precision (a-h,o-z)
      double complex dl_slr_g,dl_slr_h,dl_slr_hg,dl_slr_c,dl_slr_n
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      dl_slr_box = (ih*dl_slr_g(i,j,k,l) + ih*dl_slr_h(i,j,k,l)
     $     + ic*dl_slr_c(i,j,k,l) + in*dl_slr_n(i,j,k,l)
     $     + ih*dl_slr_hg(i,j,k,l))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^S_RL                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_srl_g(i,j,k,l)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      dl_srl_g = (0.d0,0.d0)
      return
      end

      double complex function dl_srl_h(i,j,k,l)
c     Double Higgs contribution
      implicit double precision (a-h,o-z)
      dl_srl_h = (0.d0,0.d0)
      return
      end

      double complex function dl_srl_hg(i,j,k,l)
c     W-Higgs contribution
      implicit double precision (a-h,o-z)
      double complex yh_eff_r
      double complex ckm
      common/km_mat/ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_srl_hg = (0.d0,0.d0)
      if (k.ne.l) return
      do m=1,3
         do n=1,2
            dl_srl_hg = dl_srl_hg + ckm(j,m)*yh_eff_r(i,m,n)*zh(1,n)
     $           * dp1(um(m),cm(n),wm,0.d0)
         end do
      end do
      dl_srl_hg = e2/st2/2*yl(l)*dl_srl_hg
      return
      end

      double complex function dl_srl_c(i,j,k,l)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd,zv,zl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_srl_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do nn=1,3
               do mm=1,6
                 dl_srl_c = dl_srl_c + vr_duc(i,mm,m)*zpos(1,n)*zv(l,nn)
     $              * dconjg(vl_duc(j,mm,n)*zv(k,nn))*zneg(2,m)
     $              * dp1(fcm(m),fcm(n),sum(mm),vm(nn))
               end do
            end do
         end do
      end do
      dl_srl_c = - e/st*yl(l)/2*dl_srl_c
      return
      end

      double complex function dl_srl_n(i,j,k,l)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,vl_lln,vr_lln,zn,zu,zd,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_srl_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do mm=1,6
               do nn=1,6
                  dl_srl_n = dl_srl_n
     $                 - (vl_lln(k,nn,n)*dconjg(vr_lln(l,nn,m))
     $                 + vl_lln(k,nn,m)*dconjg(vr_lln(l,nn,n)))
     $                 * vr_ddn(i,mm,m)*dconjg(vl_ddn(j,mm,n))/2
     $                 * dp1(fnm(m),fnm(n),sdm(mm),slm(nn))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dl_srl_box(i,j,k,l)
c     Full box A^S_RL formfactor
      implicit double precision (a-h,o-z)
      double complex dl_srl_g,dl_srl_h,dl_srl_hg,dl_srl_c,dl_srl_n
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      dl_srl_box = (ih*dl_srl_g(i,j,k,l) + ih*dl_srl_h(i,j,k,l)
     $     + ic*dl_srl_c(i,j,k,l) + in*dl_srl_n(i,j,k,l)
     $     + ih*dl_srl_hg(i,j,k,l))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^T_L                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_tl_g(i,j,k,l)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      dl_tl_g = (0.d0,0.d0)
      return
      end

      double complex function dl_tl_h(i,j,k,l)
c     Double Higgs contribution
      implicit double precision (a-h,o-z)
      dl_tl_h = (0.d0,0.d0)
      return
      end

      double complex function dl_tl_hg(i,j,k,l)
c     W-Higgs contribution
      implicit double precision (a-h,o-z)
      dl_tl_hg = (0.d0,0.d0)
      return
      end

      double complex function dl_tl_c(i,j,k,l)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd,zv,zl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_tl_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do nn=1,3
               do mm=1,6
                 dl_tl_c = dl_tl_c + vl_duc(i,mm,m)*zpos(1,n)*zv(l,nn)
     $                 * zneg(2,m)*dconjg(vr_duc(j,mm,n)*zv(k,nn))
     $                 * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(mm),vm(nn))
               end do
            end do
         end do
      end do
      dl_tl_c = - e/st*yl(l)/8*dl_tl_c
      return
      end

      double complex function dl_tl_n(i,j,k,l)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,vl_lln,vr_lln,zn,zu,zd,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_tl_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do mm=1,6
               do nn=1,6
                  dl_tl_n = dl_tl_n
     $                 - (vl_lln(k,nn,n)*dconjg(vr_lln(l,nn,m))
     $                 - vl_lln(k,nn,m)*dconjg(vr_lln(l,nn,n)))
     $                 * vl_ddn(i,mm,m)*dconjg(vr_ddn(j,mm,n))/8
     $                 *fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(mm),slm(nn))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dl_tl_box(i,j,k,l)
c     Full box A^T_L formfactor
      implicit double precision (a-h,o-z)
      double complex dl_tl_g,dl_tl_h,dl_tl_hg,dl_tl_c,dl_tl_n
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      dl_tl_box = (ih*dl_tl_g(i,j,k,l) + ih*dl_tl_h(i,j,k,l)
     $     + ic*dl_tl_c(i,j,k,l) + in*dl_tl_n(i,j,k,l)
     $     + ih*dl_tl_hg(i,j,k,l))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^T_R                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_tr_g(i,j,k,l)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      dl_tr_g = (0.d0,0.d0)
      return
      end

      double complex function dl_tr_h(i,j,k,l)
c     Double Higgs contribution
      implicit double precision (a-h,o-z)
      dl_tr_h = (0.d0,0.d0)
      return
      end

      double complex function dl_tr_hg(i,j,k,l)
c     W-Higgs contribution
      implicit double precision (a-h,o-z)
      dl_tr_hg = (0.d0,0.d0)
      return
      end

      double complex function dl_tr_c(i,j,k,l)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd,zv,zl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/yukawa/yl(3),yu(3),yd(3)
      dl_tr_c = (0.d0,0.d0)
      do m=1,2
         do n=1,2
            do nn=1,3
               do mm=1,6
                 dl_tr_c = dl_tr_c + vr_duc(i,mm,m)*zv(l,nn)
     $              *dconjg(vl_duc(j,mm,n)*zpos(1,m)*zneg(2,n)*zv(k,nn))
     $              * fcm(m)*fcm(n)*dp0(fcm(m),fcm(n),sum(mm),vm(nn))
               end do
            end do
         end do
      end do
      dl_tr_c = - e/st*yl(k)/8*dl_tr_c
      return
      end

      double complex function dl_tr_n(i,j,k,l)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,vl_lln,vr_lln,zn,zu,zd,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      dl_tr_n = (0.d0,0.d0)
      do m=1,4
         do n=1,4
            do mm=1,6
               do nn=1,6
                  dl_tr_n = dl_tr_n
     $                 - (vr_lln(k,nn,n)*dconjg(vl_lln(l,nn,m))
     $                 - vr_lln(k,nn,m)*dconjg(vl_lln(l,nn,n)))
     $                 * vr_ddn(i,mm,m)*dconjg(vl_ddn(j,mm,n))/8
     $                 *fnm(m)*fnm(n)*dp0(fnm(m),fnm(n),sdm(mm),slm(nn))
               end do
            end do
         end do
      end do
      return
      end

      double complex function dl_tr_box(i,j,k,l)
c     Full box A^T_R formfactor
      implicit double precision (a-h,o-z)
      double complex dl_tr_g,dl_tr_h,dl_tr_hg,dl_tr_c,dl_tr_n
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      dl_tr_box = (ih*dl_tr_g(i,j,k,l) + ih*dl_tr_h(i,j,k,l)
     $     + ic*dl_tr_c(i,j,k,l) + in*dl_tr_n(i,j,k,l)
     $     + ih*dl_tr_hg(i,j,k,l))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Sum of the box, Z- and Higgs penguin contributions           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function dl_vll(i,j,k,l)
c     Full A^V_LL formfactor
      implicit double precision (a-h,o-z)
      double complex dl_vll_box,zdd_vl,dv_sig0,da_sig0
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      if (ibox.eq.1) then 
         dl_vll = dl_vll_box(i,j,k,l)
      else
         dl_vll = (0.d0,0.d0)
      end if
      if ((k.eq.l).and.(izpeng.eq.1)) then
         dl_vll = dl_vll - e*(1 - 2*st2)/2/sct/zm2*(zdd_vl(i,j)
     $        - e/2/sct*(1 - 2.d0/3*st2)*(dv_sig0(i,j) - da_sig0(i,j)))
      end if
      return
      end

      double complex function dl_vrr(i,j,k,l)
c     Full A^V_RR formfactor
      implicit double precision (a-h,o-z)
      double complex dl_vrr_box,zdd_vr,dv_sig0,da_sig0
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      if (ibox.eq.1) then 
         dl_vrr = dl_vrr_box(i,j,k,l)
      else
         dl_vrr = (0.d0,0.d0)
      end if
      if ((k.eq.l).and.(izpeng.eq.1)) then
         dl_vrr = dl_vrr + e*st/ct/zm2*(zdd_vr(i,j)
     $        + e*st/3/ct*(dv_sig0(i,j) + da_sig0(i,j)))
      end if
      return
      end

      double complex function dl_vlr(i,j,k,l)
c     Full A^V_LR formfactor
      implicit double precision (a-h,o-z)
      double complex dl_vlr_box,zdd_vl,dv_sig0,da_sig0
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      if (ibox.eq.1) then 
         dl_vlr = dl_vlr_box(i,j,k,l)
      else
         dl_vlr = (0.d0,0.d0)
      end if
      if ((k.eq.l).and.(izpeng.eq.1)) then
         dl_vlr = dl_vlr + e*st/ct/zm2*(zdd_vl(i,j)
     $        - e/2/sct*(1 - 2.d0/3*st2)*(dv_sig0(i,j) - da_sig0(i,j)))
      end if
      return
      end

      double complex function dl_vrl(i,j,k,l)
c     Full A^V_RL formfactor
      implicit double precision (a-h,o-z)
      double complex dl_vrl_box,zdd_vr,dv_sig0,da_sig0
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      if (ibox.eq.1) then 
         dl_vrl = dl_vrl_box(i,j,k,l)
      else
         dl_vrl = (0.d0,0.d0)
      end if
      if ((k.eq.l).and.(izpeng.eq.1)) then
         dl_vrl = dl_vrl - e*(1 - 2*st2)/2/sct/zm2*(zdd_vr(i,j)
     $        + e*st/3/ct*(dv_sig0(i,j) + da_sig0(i,j)))
      end if
      return
      end

      double complex function dl_sll(i,j,k,l)
c     Full A^S_LL formfactor
      implicit double precision (a-h,o-z)
      double complex dl_sll_box,form(2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      if (ibox.eq.1) then 
         dl_sll = dl_sll_box(i,j,k,l)
      else
         dl_sll = (0.d0,0.d0)
      end if
      if ((k.eq.l).and.(ihpeng.eq.1)) then
         do m=1,2
            call sdd_vert_eff(i,j,m,form)
            dl_sll = dl_sll + yl(k)*zr(1,m)/sq2/rm(m)/rm(m)*form(1)
            call pdd_vert_eff(i,j,m,form)
            dl_sll = dl_sll + yl(k)*zh(1,m)/sq2/pm(m)/pm(m)*form(1)
         end do
      end if
      return
      end

      double complex function dl_srr(i,j,k,l)
c     Full A^S_RR formfactor
      implicit double precision (a-h,o-z)
      double complex dl_srr_box,form(2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      if (ibox.eq.1) then 
         dl_srr = dl_srr_box(i,j,k,l)
      else
         dl_srr = (0.d0,0.d0)
      end if
      if ((k.eq.l).and.(ihpeng.eq.1)) then
         do m=1,2
            call sdd_vert_eff(i,j,m,form)
            dl_srr = dl_srr + yl(k)*zr(1,m)/sq2/rm(m)/rm(m)*form(2)
            call pdd_vert_eff(i,j,m,form)
            dl_srr = dl_srr - yl(k)*zh(1,m)/sq2/pm(m)/pm(m)*form(2)
         end do
      end if
      return
      end

      double complex function dl_slr(i,j,k,l)
c     Full A^S_LR formfactor
      implicit double precision (a-h,o-z)
      double complex dl_slr_box,form(2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      if (ibox.eq.1) then 
         dl_slr = dl_slr_box(i,j,k,l)
      else
         dl_slr = (0.d0,0.d0)
      end if
      if ((k.eq.l).and.(ihpeng.eq.1)) then
         do m=1,2
            call sdd_vert_eff(i,j,m,form)
            dl_slr = dl_slr + yl(k)*zr(1,m)/sq2/rm(m)/rm(m)*form(1)
            call pdd_vert_eff(i,j,m,form)
            dl_slr = dl_slr - yl(k)*zh(1,m)/sq2/pm(m)/pm(m)*form(1)
         end do
      end if
      return
      end

      double complex function dl_srl(i,j,k,l)
c     Full A^S_RL formfactor
      implicit double precision (a-h,o-z)
      double complex dl_srl_box,form(2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      if (ibox.eq.1) then 
         dl_srl = dl_srl_box(i,j,k,l)
      else
         dl_srl = (0.d0,0.d0)
      end if
      if ((k.eq.l).and.(ihpeng.eq.1)) then
         do m=1,2
            call sdd_vert_eff(i,j,m,form)
            dl_srl = dl_srl + yl(k)*zr(1,m)/sq2/rm(m)/rm(m)*form(2)
            call pdd_vert_eff(i,j,m,form)
            dl_srl = dl_srl + yl(k)*zh(1,m)/sq2/pm(m)/pm(m)*form(2)
         end do
      end if
      return
      end

      double complex function dl_tl(i,j,k,l)
c     Full A^T_L formfactor
      implicit double precision (a-h,o-z)
      double complex dl_tl_box
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      if (ibox.eq.1) then 
         dl_tl = dl_tl_box(i,j,k,l)
      else
         dl_tl = (0.d0,0.d0)
      end if
      return
      end

      double complex function dl_tr(i,j,k,l)
c     Full A^T_R formfactor
      implicit double precision (a-h,o-z)
      double complex dl_tr_box
      common/diag_type/ihpeng,izpeng,ifpeng,ibox
      if (ibox.eq.1) then 
         dl_tr = dl_tr_box(i,j,k,l)
      else
         dl_tr = (0.d0,0.d0)
      end if
      return
      end



