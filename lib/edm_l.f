c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
      
c     FILENAME: EDM_L.F
c     Released: 25.03.1998 (J.R.)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for EDM of leptons             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      double precision function edm_l_c(i)
c     chargino-sneutrino contributions
      implicit double precision (a-h,o-z)
      double complex zpos,zneg,zv,zl
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
       common/yukawa/yl(3),yu(3),yd(3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      edm_l_c = 0
      do k=1,3
        do l=1,2
          edm_l_c  = edm_l_c + fcm(l)*abs(zv(i,k))**2
     $          *dimag(zpos(1,l)*zneg(2,l))*cp11(fcm(l),vm(k))
        end do
      end do
      edm_l_c = alpha/2/pi/st*yl(i)*edm_l_c
      return
      end
      
      double precision function edm_l_n(i)
c     neutralino-slepton squark contributions
      implicit double precision (a-h,o-z)
      double complex vl_lln,vr_lln,zn,zv,zl
      common/neut/fnm(4),zn(4,4)
      common/slmass/vm(3),slm(6),zv(3,3),zl(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      edm_l_n = 0
      do k=1,6
        do l=1,4
          edm_l_n = edm_l_n + fnm(l)*cp12(slm(k),fnm(l))
     $        *dimag(vl_lln(i,k,l)*dconjg(vr_lln(i,k,l)))
        end do
      end do
      edm_l_n =  - e/16/pi/pi*edm_l_n
      return
      end
      
      double precision function edm_l(i)
c     Full lepton EDM
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/ph_units/hbar,gev
      external init_units
      edm_l = (edm_l_c(i) + edm_l_n(i))/e/gev
      return
      end
      
