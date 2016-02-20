c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
      
c     FILENAME: EDM_U.F
c     Released: 25.03.1998 (J.R.)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for EDM of u quarks            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      double precision function edm_u_c(i)
c     chargino-down squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_udc,vr_udc,zpos,zneg,zu,zd
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      edm_u_c = 0
      do k=1,6
        do l=1,2
          edm_u_c  = edm_u_c 
     $        + fcm(l)*dimag(vl_udc(i,k,l)*dconjg(vr_udc(i,k,l)))
     $        * (cp11(fcm(l),sdm(k)) + cp12(sdm(k),fcm(l))/6)
        end do
      end do
      edm_u_c = - e/8/pi/pi*edm_u_c 
      return
      end
      
      double precision function edm_u_n(i)
c     neutralino-up squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_uun,vr_uun,zn,zu,zd
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      edm_u_n = 0
      do k=1,6
        do l=1,4
          edm_u_n = edm_u_n + fnm(l)*cp12(sum(k),fnm(l))
     $        * dimag(vl_uun(i,k,l)*dconjg(vr_uun(i,k,l)))
        end do
      end do
      edm_u_n = e/24/pi/pi*edm_u_n 
      return
      end
      
      double precision function edm_u_g(i)
c     gluino-up squark contributions
      implicit double precision (a-h,o-z)
      double complex zu,zd,gm2,gm3
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/gmass/gm1,gm2,gm3
      edm_u_g = 0
      do k=1,6
        edm_u_g = edm_u_g 
     $      + dimag(zu(i,k)*dconjg(zu(i+3,k)))*cp12(sum(k),gm1)
      end do
      edm_u_g = 4/9.d0/pi*e*alfas(zm)*gm1*edm_u_g 
      return
      end
      
      double precision function edm_u(i)
c     Full up quark EDM
      implicit double precision (a-h,o-z)
      edm_u = edm_u_n(i) + edm_u_c(i) + edm_u_g(i)
      return
      end


