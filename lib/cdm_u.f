c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
      
c     FILENAME: CDM_U.F
c     Released: 25.03.1998 (J.R.)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for CDM of u quarks            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      double precision function cdm_u_c(i)
c     chargino-down squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_udc,vr_udc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      cdm_u_c = 0
      do k=1,6
         do l=1,2
            cdm_u_c  = cdm_u_c + fcm(l)*cp12(sdm(k),fcm(l))
     $           * dimag(vl_udc(i,k,l)*dconjg(vr_udc(i,k,l))) 
         end do
      end do
      cdm_u_c = cdm_u_c/8/pi 
      return
      end
      
      double precision function cdm_u_n(i)
c     neutralino-up squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      cdm_u_n = 0
      do k=1,6
        do l=1,4
          cdm_u_n = cdm_u_n + fnm(l)*cp12(sum(k),fnm(l))
     $        * dimag(vl_uun(i,k,l)*dconjg(vr_uun(i,k,l)))
        end do
      end do
      cdm_u_n = cdm_u_n/8/pi
      return
      end
      
      double precision function cdm_u_g(i)
c     gluino-up squark contributions
      implicit double precision (a-h,o-z)
      double complex zu,zd,gm2,gm3
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      cdm_u_g = 0
      do k=1,6
         cdm_u_g = cdm_u_g + dimag(zu(i,k)*dconjg(zu(i+3,k)))
     $        * (3*cp11(gm1,sum(k)) + cp12(sum(k),gm1)/6)
      end do
      cdm_u_g = - alfas(zm)*gm1*cdm_u_g 
      return
      end
      
      double precision function cdm_u(i)
c     Full up quark CDM
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      cdm_u = sqrt(alfas(zm)/pi)*(cdm_u_n(i) + cdm_u_c(i) + cdm_u_g(i))
      return
      end


