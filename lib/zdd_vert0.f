c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: ZDD_VERT0.F
c     Released: 07:08:2007 (J.R.)
c     Corrected: 22.08.2007(A.Dedes, P.Tanedo.,J.R.)
c     Gauge + Higgs contributions rewritten
c     Corrected: 02:11:2007
c     sign error in zdd_vl_hg corrected

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expression for the coefficient of the        c
c     Hamiltonian for the B->pi+v+vbar and the K->pi+v+vbar           c
c     processes (at the quark level).                                 c
c     General form of the vertex is:                                  c
c     i\bar d^J (A^V_L \gamma_mu P_L + A^V_R \gamma_mu P_R) d^I Z^mu  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_L                                            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function zdd_vl_g(i,j)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      double complex ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      zdd_vl_g = (0.d0,0.d0)
      do m=1,3
         zdd_vl_g = zdd_vl_g + dconjg(ckm(i,m))*ckm(j,m)
c     Wuu contribution
     $        * ((1 - 4.d0/3*st2)*(cp1(wm,um(m),um(m)) - 0.5d0)
     $        + 8.d0/3*st2*um(m)*um(m)*cp0(wm,um(m),um(m))
c     WWu contribution
     $        + ct2*(6*cp1(wm,wm,um(m)) - 1))
      end do
      zdd_vl_g = e*e2/4/sct/st2*zdd_vl_g
      return
      end
      
      double complex function zdd_vl_hg(i,j)
c     Higgs-gauge contribution
      implicit double precision (a-h,o-z)
      double complex ckm,yh_eff_l
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      zdd_vl_hg = (0.d0,0.d0)
      do m=1,3
         zdd_vl_hg = zdd_vl_hg + (dconjg(ckm(i,m)*yh_eff_l(j,m,2)) 
     $        + ckm(j,m)*yh_eff_l(i,m,2))*um(m)*cp0(wm,wm,um(m))
      end do
      zdd_vl_hg = e2*wm/sq2/ct*zdd_vl_hg
      return
      end
      
      double complex function zdd_vl_h(i,j)
c     Higgs contribution
      implicit double precision (a-h,o-z)
      double complex yh_eff_l
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      zdd_vl_h = (0.d0,0.d0)
      do m=1,3
         do l=1,2
           zdd_vl_h = zdd_vl_h - yh_eff_l(i,m,l)*dconjg(yh_eff_l(j,m,l))
c     Huu contribution
     $           * (2*st2/3.d0*(cp1(cm(l),um(m),um(m)) - 0.5d0)
     $           + (1 - 4*st2/3.d0)*um(m)*um(m)*cp0(cm(l),um(m),um(m))
c     HHu contribution
     $           + (ct2 - st2)/2*(cp1(cm(l),cm(l),um(m)) + 0.5d0))
        end do
      end do
      zdd_vl_h = e/2/sct*zdd_vl_h
      return
      end
      
      double complex function zdd_vl_c(i,j)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vl_duc,vl_ccz,vr_ccz,v_uuz,zpos,zneg,zu,zd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      zdd_vl_c = (0.d0,0.d0)
c     CCU contribution
      do m=1,2
        do n=1,2
          do mm=1,6
             zdd_vl_c = zdd_vl_c - vl_duc(i,mm,m)*dconjg(vl_duc(j,mm,n))
     $         * (vl_ccz(m,n)*(cp1(sum(mm),fcm(m),fcm(n)) - 0.5d0)
     $         - 2*fcm(m)*fcm(n)*vr_ccz(m,n)*cp0(sum(mm),fcm(m),fcm(n)))
            end do
         end do
      end do
c     CUU contribution
      do mm=1,6
         do nn=1,6
            do m=1,2  
               zdd_vl_c = zdd_vl_c + v_uuz(nn,mm)*vl_duc(i,mm,m)
     $              * dconjg(vl_duc(j,nn,m))
     $              * (cp1(fcm(m),sum(mm),sum(nn)) + 0.5d0)
            end do
         end do
      end do
      zdd_vl_c = e/4/sct*zdd_vl_c
      return
      end

      double complex function zdd_vl_n(i,j)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vl_nnz,v_ddz,zn,zu,zd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      zdd_vl_n = (0.d0,0.d0)
c     NND contribution
      do m=1,4
        do n=1,4
          do mm=1,6
            zdd_vl_n = zdd_vl_n + vl_ddn(i,mm,m)*dconjg(vl_ddn(j,mm,n))
     $            * (dconjg(vl_nnz(n,m))
     $            * (cp1(sdm(mm),fnm(m),fnm(n)) - 0.5d0)
     $            + 2*vl_nnz(n,m)*fnm(m)*fnm(n)
     $            * cp0(sdm(mm),fnm(m),fnm(n)))
            end do
         end do
      end do
c     NDD contribution
      do mm=1,6
         do nn=1,6
            do m=1,4
               zdd_vl_n = zdd_vl_n - v_ddz(mm,nn)*vl_ddn(i,mm,m)
     $              * dconjg(vl_ddn(j,nn,m))
     $              * (cp1(fnm(m),sdm(mm),sdm(nn)) + 0.5d0)
            end do
         end do
      end do
      zdd_vl_n = e/4/sct*zdd_vl_n
      return
      end

      double complex function zdd_vl_gl(i,j)
c     Gluino contribution
      implicit double precision (a-h,o-z)
      double complex zu,zd,gm2,gm3,v_ddz
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      zdd_vl_gl = (0.d0,0.d0)
      g3 = 4*pi*alfas(gm1)
c     DDG contribution only
      do mm=1,6
         do nn=1,6
            zdd_vl_gl = zdd_vl_gl - v_ddz(mm,nn)*zd(i,mm)
     $           * dconjg(zd(j,nn))*(cp1(gm1,sdm(mm),sdm(nn)) + 0.5d0)
         end do
      end do
      zdd_vl_gl = 2*e*g3/3.d0/sct*zdd_vl_gl
      return
      end

      double complex function zdd_vl(i,j)
c     Full A^V_L formfactor
      implicit double precision (a-h,o-z)
      double complex zdd_vl_g,zdd_vl_hg,zdd_vl_h,zdd_vl_c,zdd_vl_n,
     $     zdd_vl_gl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      zdd_vl = (ih*(zdd_vl_g(i,j) + zdd_vl_hg(i,j) + zdd_vl_h(i,j))
     $     + ic*zdd_vl_c(i,j) + in*zdd_vl_n(i,j) 
     $     + ig*zdd_vl_gl(i,j))/16/pi/pi
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Formfactor A^V_R                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      double complex function zdd_vr_g(i,j)
c     Gauge contribution
      implicit double precision (a-h,o-z)
      zdd_vr_g = (0.d0,0.d0)
      return
      end
      
      double complex function zdd_vr_hg(i,j)
c     Higgs-gauge contribution
      implicit double precision (a-h,o-z)
      zdd_vr_hg = (0.d0,0.d0)
      return
      end
      
      double complex function zdd_vr_h(i,j)
c     Higgs contribution
      implicit double precision (a-h,o-z)
      double complex yh_eff_r
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      zdd_vr_h = (0.d0,0.d0)
      do m=1,3
         do l=1,2
           zdd_vr_h = zdd_vr_h + yh_eff_r(i,m,l)*dconjg(yh_eff_r(j,m,l))
c     Huu contribution
     $           * (((1 - 4*st2/3.d0)*(cp1(cm(l),um(m),um(m)) - 0.5d0)
     $           + 2*st2/3.d0*um(m)*um(m)*cp0(cm(l),um(m),um(m)))
c     HHu contribution
     $           - (ct2 - st2)*(cp1(cm(l),cm(l),um(m)) + 0.5d0))
        end do
      end do
      zdd_vr_h = e/4/sct*zdd_vr_h
      return
      end
      
      double complex function zdd_vr_c(i,j)
c     Chargino contribution
      implicit double precision (a-h,o-z)
      double complex vr_duc,vl_ccz,vr_ccz,v_uuz,zpos,zneg,zu,zd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      zdd_vr_c = (0.d0,0.d0)
c     CCU contribution
      do m=1,2
        do n=1,2
          do mm=1,6
             zdd_vr_c = zdd_vr_c - vr_duc(i,mm,m)*dconjg(vr_duc(j,mm,n))
     $         * (vr_ccz(m,n)*(cp1(sum(mm),fcm(m),fcm(n)) - 0.5d0)
     $         - 2*fcm(m)*fcm(n)*vl_ccz(m,n)*cp0(sum(mm),fcm(m),fcm(n)))
            end do
         end do
      end do
c     CUU contribution
      do mm=1,6
         do nn=1,6
            do m=1,2  
               zdd_vr_c = zdd_vr_c + v_uuz(nn,mm)*vr_duc(i,mm,m)
     $              * dconjg(vr_duc(j,nn,m))
     $              * (cp1(fcm(m),sum(mm),sum(nn)) + 0.5d0)
            end do
         end do
      end do
      zdd_vr_c = e/4/sct*zdd_vr_c
      return
      end

      double complex function zdd_vr_n(i,j)
c     Neutralino contribution
      implicit double precision (a-h,o-z)
      double complex vr_ddn,vl_nnz,v_ddz,zn,zu,zd
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      zdd_vr_n = (0.d0,0.d0)
c     NND contribution
      do m=1,4
        do n=1,4
           do mm=1,6
             zdd_vr_n = zdd_vr_n - vr_ddn(i,mm,m)*dconjg(vr_ddn(j,mm,n))
     $             * (vl_nnz(n,m)*(cp1(sdm(mm),fnm(m),fnm(n)) - 0.5d0)
     $             + 2*dconjg(vl_nnz(n,m))*fnm(m)*fnm(n)
     $             * cp0(sdm(mm),fnm(m),fnm(n)))
            end do
         end do
      end do
c     NDD contribution
      do mm=1,6
         do nn=1,6
            do m=1,4
               zdd_vr_n = zdd_vr_n - v_ddz(mm,nn)*vr_ddn(i,mm,m)
     $              * dconjg(vr_ddn(j,nn,m))*cp1(fnm(m),sdm(mm),sdm(nn))
            end do
         end do
      end do
      zdd_vr_n = e/4/sct*zdd_vr_n
      return
      end

      double complex function zdd_vr_gl(i,j)
c     Gluino contribution
      implicit double precision (a-h,o-z)
      double complex zu,zd,v_ddz,gm2,gm3
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      zdd_vr_gl = (0.d0,0.d0)
      g3 = 4*pi*alfas(gm1)
c     DDG contribution only
      do mm=1,6
         do nn=1,6
            zdd_vr_gl = zdd_vr_gl - v_ddz(mm,nn)*zd(i+3,mm)
     $           * dconjg(zd(j+3,nn))*(cp1(gm1,sdm(mm),sdm(nn)) + 0.5d0)
         end do
      end do
      zdd_vr_gl = 2*e*g3/3.d0/sct*zdd_vr_gl
      return
      end
      
      double complex function zdd_vr(i,j)
c     Full A^V_R formfactor
      implicit double precision (a-h,o-z)
      double complex zdd_vr_g,zdd_vr_hg,zdd_vr_h,zdd_vr_c,zdd_vr_n,
     $     zdd_vr_gl
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      zdd_vr = (ih*(zdd_vr_g(i,j) + zdd_vr_hg(i,j) + zdd_vr_h(i,j)) 
     $     + ic*zdd_vr_c(i,j) + in*zdd_vr_n(i,j) 
     $     + ig*zdd_vr_gl(i,j))/16/pi/pi
      return
      end


      
