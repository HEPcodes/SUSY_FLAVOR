c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: DD_GAMMA.F
c     Released: 25:03:1996 (J.R.)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for the coefficients of the   c
c     Hamiltonian for the d^J -> d^I + gamma decay, e.g.           c
c     b -> s gamma at the MZ energy scale                          c
c     General form of the Hamiltonian is:                          c
c     -iH = SUM_{i=1}^5 (A^i_L H^i_L + A^i_R H^i_R)                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dd_gam_w(i,j,cfl,cfr)
c     W and u quark contributions (SM like)
      implicit double precision (a-h,o-z)
      double complex ckm
      double complex cfl(5),cfr(5)
      double complex ai,aj
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      qf = - 2.d0/3
      qv = 1
      do k=1,3
        ai = e/sq2/st*ckm(i,k)
        aj = e/sq2/st*dconjg(ckm(j,k))
        call dd_ffv(ai,aj,qf,um(k),wm,cfl,cfr)
        call dd_vvf(ai,aj,qv,um(k),wm,cfl,cfr)
      end do
      return
      end

      subroutine dd_gam_h(i,j,cfl,cfr)
c     Higgs and Goldstone contributions
      implicit double precision (a-h,o-z)
      double complex ckm
      double complex cfl(5),cfr(5)
      double complex ai,aj,bi,bj
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/hangle/ca,sa,cb,sb
      qf = - 2.d0/3
      qs = - 1
      do k=1,3
        do l=1,2
          ai = e/sq2/st*um(k)/wm*zh(2,l)/sb*ckm(i,k)
          bi = e/sq2/st*dm(i)/wm*zh(1,l)/cb*ckm(i,k)
          aj = e/sq2/st*um(k)/wm*zh(2,l)/sb*dconjg(ckm(j,k))
          bj = e/sq2/st*dm(j)/wm*zh(1,l)/cb*dconjg(ckm(j,k))
          call dd_ffs(ai,aj,bi,bj,qf,um(k),cm(l),cfl,cfr)
          call dd_ssf(ai,aj,bi,bj,qs,um(k),cm(l),cfl,cfr)
        end do
      end do
      return
      end

      subroutine dd_gam_c(i,j,cfl,cfr)
c     chargino/up squark contributions
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ai,aj,bi,bj
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      qf = 1
      qs = 2.d0/3
      do k=1,6
        do l=1,2
          ai = dconjg(vl_duc(i,k,l))
          bi = dconjg(vr_duc(i,k,l))
          aj = vl_duc(j,k,l)
          bj = vr_duc(j,k,l)
          call dd_ffs(ai,aj,bi,bj,qf,fcm(l),sum(k),cfl,cfr)
          call dd_ssf(ai,aj,bi,bj,qs,fcm(l),sum(k),cfl,cfr)
        end do
      end do
      return
      end

      subroutine dd_gam_n(i,j,cfl,cfr)
c     neutralino-down squark contributions
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ai,aj,bi,bj
      double complex vl_ddn,vr_ddn,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      qs = - 1.d0/3
      do k=1,6
        do l=1,4
          ai = dconjg(vl_ddn(i,k,l))
          bi = dconjg(vr_ddn(i,k,l))
          aj = vl_ddn(j,k,l)
          bj = vr_ddn(j,k,l)
          call dd_ssf(ai,aj,bi,bj,qs,fnm(l),sdm(k),cfl,cfr)
        end do
      end do
      return
      end

      subroutine dd_gam_g(i,j,cfl,cfr)
c     gluino-down squark contributions
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      double complex ai,aj,bi,bj
      double complex zu,zd,gm2,gm3
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      do k=1,6
         qs = 32/9.d0*pi*alfas((gm1 + sdm(k))/2)
         ai = - dconjg(zd(i,k))
         bi =   dconjg(zd(i+3,k))
         aj = - zd(j,k)
         bj =   zd(j+3,k)
         call dd_ssf(ai,aj,bi,bj,qs,gm1,sdm(k),cfl,cfr)
      end do
      return
      end

      subroutine dd_gam(i,j,cfl,cfr)
c     Full coefficients
      implicit double precision (a-h,o-z)
      double complex cfl(5),cfr(5)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do k=1,5
        cfl(k) = (0.d0,0.d0)
        cfr(k) = (0.d0,0.d0)
      end do
      call dd_gam_w(i,j,cfl,cfr)
      call dd_gam_h(i,j,cfl,cfr)
      call dd_gam_c(i,j,cfl,cfr)
      call dd_gam_n(i,j,cfl,cfr)
      call dd_gam_g(i,j,cfl,cfr)
      do k=1,5
        cfl(k) = cfl(k)/16/pi/pi
        cfr(k) = cfr(k)/16/pi/pi
      end do
      return
      end


