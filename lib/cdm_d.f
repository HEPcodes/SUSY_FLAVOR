c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
      
c     FILENAME: CDM_D.F
c     Released: 25.03.1998 (J.R.)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for CDM of d quarks            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      double precision function cdm_d_c(i)
c     chargino-up squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_duc,vr_duc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      cdm_d_c = 0
      do k=1,6
         do l=1,2
            cdm_d_c  = cdm_d_c + fcm(l)*cp12(sum(k),fcm(l))
     $           * dimag(vl_duc(i,k,l)*dconjg(vr_duc(i,k,l)))
         end do
      end do
      cdm_d_c = cdm_d_c/8/pi 
      return
      end
      
      double precision function cdm_d_n(i)
c     neutralino-down squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_ddn,vr_ddn,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      cdm_d_n = 0
      do k=1,6
        do l=1,4
          cdm_d_n = cdm_d_n + fnm(l)*cp12(sdm(k),fnm(l))
     $        * dimag(vl_ddn(i,k,l)*dconjg(vr_ddn(i,k,l)))
        end do
      end do
      cdm_d_n = cdm_d_n/8/pi
      return
      end
      
      double precision function cdm_d_g(i)
c     gluino-up squark contributions
      implicit double precision (a-h,o-z)
      double complex zu,zd,gm2,gm3
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      cdm_d_g = 0
      do k=1,6
         cdm_d_g = cdm_d_g + dimag(zd(i,k)*dconjg(zd(i+3,k)))
     $        * (3*cp11(gm1,sdm(k)) + cp12(sdm(k),gm1)/6)
      end do
      cdm_d_g = alfas(zm)*gm1*cdm_d_g 
      return
      end
      
      double precision function cdm_d(i)
c     Full down quark CDM
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      cdm_d = sqrt(alfas(zm)/pi)*(cdm_d_n(i) + cdm_d_c(i) + cdm_d_g(i))
      return
      end




