c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
      
c     FILENAME: EDM_N.F
c     Released: 25.03.1998 (J.R.)
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for EDM of neutron             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      double precision function edm_n()
c     complete neutron EDM. Included:
c       - electric dipole moment of quarks
c       - chromoelectric dipole moment of quarks
c       - chromoelectric dipole moment of gluon
c       - approximate QCD renormalization factors
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/edm_qcd/eta_e,eta_c,eta_g,alamx
      common/ph_units/hbar,gev
      external init_units
      edm_n = (eta_e/3.d0*(4*edm_d(1) - edm_u(1))
     $     + e/12.d0/pi*eta_c*(4*cdm_d(1) - cdm_u(1))
     $     + e/4.d0/pi*eta_g*alamx*cdm_g())/e/gev
      return
      end

      double precision function edm_d_phys(i)
c     d quark effective EDM
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/edm_qcd/eta_e,eta_c,eta_g,alamx
      common/ph_units/hbar,gev
      edm_d_phys = (eta_e*edm_d(i) + e/4/pi*eta_c*cdm_d(i))/e/gev
      return
      end 

      double precision function edm_u_phys(i)
c     u quark effective EDM
      implicit double precision (a-h,o-z)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/edm_qcd/eta_e,eta_c,eta_g,alamx
      common/ph_units/hbar,gev
      edm_u_phys = (eta_e*edm_u(i) + e/4/pi*eta_c*cdm_u(i))/e/gev
      return
      end 
      
      block data init_edm
      implicit double precision (a-h,o-z)
c     QCD correction factors and chiral symmetry breaking scale
      common/edm_qcd/eta_e,eta_c,eta_g,alamx
      data eta_e,eta_c,eta_g/1.53d0,3.4d0,3.4d0/
      data alamx/1.18d0/
      end


