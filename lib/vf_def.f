c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM} 
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: VF_DEF.F
c     Released: 5: 1:1993(J.R.)
c     Last revised: 28: 1:1993 (P.Ch.)
c     Vertices vl(r)_lsnc (lepton-sneutrino-chargino) added
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for fermion vertices              c
c     Convention is such that the Feynman rule for the vertex is       c
c     (-i)*(expression calculated below)*(projector P_l or P_R)        c
c     Compare with the paper: J.Rosiek@Phys.Rev.D41(1990)p.3464;       c
c     erratum, hep-ph/9511250                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine init_yukawa()
c     Initialization of Yukawa coupling
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/vev/v1,v2
      do i=1,3
         yu(i) = sq2*um(i)/v2
         yd(i) = - sq2*dm(i)/v1
         yl(i) = - sq2*em(i)/v1
      end do
      return
      end

      double complex function yh_eff_l(i,j,k)
      implicit double precision (a-h,o-z)
c     temporary simplified version, no large tan(beta) resummation
      double complex ckm
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      yh_eff_l = yu(j)*zh(2,k)*dconjg(ckm(i,j)) 
      return
      end

      double complex function yh_eff_r(i,j,k)
      implicit double precision (a-h,o-z)
c     temporary simplified version, no large tan(beta) resummation
      double complex ckm
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      yh_eff_r = - yd(i)*zh(1,k)*dconjg(ckm(i,j))
      return
      end


cccccccccccccccccccccccccccccccccc
c     SUSY vertices              c
cccccccccccccccccccccccccccccccccc

      double complex function v_nnn(i,j,k)
c     Neutrino-sneutrino-neutralino vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/neut/fnm(4),zn(4,4)
      v_nnn = -e/sct/sq2*dconjg(zv(i,j))*(zn(1,k)*st - zn(2,k)*ct)
      return
      end
 
      double complex function vl_lln(i,j,k)
c     Lepton-slepton-neutralino left vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vl_lln = -e/sct/sq2*zl(i,j)*(zn(1,k)*st + zn(2,k)*ct)
     $     - yl(i)*zl(i+3,j)*zn(3,k)
      return
      end
 
      double complex function vr_lln(i,j,k)
c     Lepton-slepton-neutralino right vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vr_lln = sq2*e/ct*zl(i+3,j)*dconjg(zn(1,k))
     $     - yl(i)*zl(i,j)*dconjg(zn(3,k))
      return
      end
 
      double complex function vl_lsnc(i,j,k)
c     Lepton-sneutrino-chargino left vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zpos,zneg
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      vl_lsnc = e/st*dconjg(zv(i,j))*zpos(1,k)
      return
      end
 
      double complex function vr_lsnc(i,j,k)
c     Lepton-sneutrino-chargino right vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zpos,zneg
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      vr_lsnc = yl(i)*dconjg(zv(i,j)*zneg(2,k))
      return
      end
 
      double complex function v_nlc(i,j,k)
c     Neutrino-slepton-chargino vertex
      implicit double precision (a-h,o-z)
      double complex zv,zl,zpos,zneg
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      v_nlc = e/st*zl(i,j)*zneg(1,k) + yl(i)*zl(i+3,j)*zneg(2,k)
      return
      end
 
      double complex function vl_uun(i,j,k)
c     Up quark-up squark-neutralino left vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vl_uun =  e/sct/sq2*dconjg(zu(i,j))*(zn(1,k)*st/3 + zn(2,k)*ct)
     $     + yu(i)*dconjg(zu(i+3,j))*zn(4,k)
      return
      end
 
      double complex function vr_uun(i,j,k)
c     Up quark-up squark-neutralino right vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vr_uun = -2*sq2*e/3.d0/ct*dconjg(zu(i+3,j)*zn(1,k))
     $     + yu(i)*dconjg(zu(i,j)*zn(4,k))
      return
      end
 
      double complex function vl_ddn(i,j,k)
c     Down quark-down squark-neutralino left vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vl_ddn =  e/sct/sq2*zd(i,j)*(zn(1,k)*st/3 - zn(2,k)*ct)
     1       - yd(i)*zd(i+3,j)*zn(3,k)
      return
      end
 
      double complex function vr_ddn(i,j,k)
c     Down quark-down squark-neutralino right vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zn
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      common/yukawa/yl(3),yu(3),yd(3)
      vr_ddn = sq2*e/3/ct*zd(i+3,j)*dconjg(zn(1,k))
     $     - yd(i)*zd(i,j)*dconjg(zn(3,k))
      return
      end
 
      double complex function vl_duc(i,j,k)
c     Down quark-up squark-chargino left vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zpos,zneg
      double complex ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      vl_duc = (0.d0,0.d0)
      do l=1,3
        vl_duc = vl_duc 
     $        + dconjg(ckm(i,l))*(e/st*dconjg(zu(l,j))*zpos(1,k)
     $        - yu(l)*dconjg(zu(l+3,j))*zpos(2,k))
      end do
      return
      end
 
      double complex function vr_duc(i,j,k)
c     Down quark-up squark-chargino right vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zpos,zneg
      double complex ckm
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      vr_duc = (0.d0,0.d0)
      do l=1,3
         vr_duc = vr_duc + yd(i)*dconjg(zu(l,j)*zneg(2,k)*ckm(i,l))
      end do
      return
      end
 
      double complex function vl_udc(i,j,k)
c     Up quark-down squark-chargino left vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zpos,zneg
      double complex ckm
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      vl_udc = (0.d0,0.d0)
      do l=1,3
         vl_udc = vl_udc + (e/st*zd(l,j)*zneg(1,k)
     $        + yd(l)*zd(l+3,j)*zneg(2,k))*ckm(l,i)
      end do
      return
      end
 
      double complex function vr_udc(i,j,k)
c     Up quark-down squark-chargino right vertex
      implicit double precision (a-h,o-z)
      double complex zu,zd,zpos,zneg
      double complex ckm
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      vr_udc = (0.d0,0.d0)
      do l=1,3
         vr_udc = vr_udc - yu(i)*zd(l,j)*dconjg(zpos(2,k))*ckm(l,i)
      end do
      return
      end



