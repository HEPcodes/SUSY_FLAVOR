c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: SDD_VERT0.F
c     Released: 12:09:2008(J.R.)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains neutral scalar Higgs - down quark-down quark c
c     vertex formfactors at vanishing external momenta                c
c     (1PI irreducible contributions only)                            c
c     Currently it calculates correctly only flavour violating        c
c     couplings, i.e. I<>J. No renormalization required in this case  c
c                                                                     c
c                                                                     c
c                              V3                                     c
c                               ------<------ d^J                     c
c                             /|     -q                               c
c                        L1 /  |                                      c
c                         /    |                                      c
c              S_k      /      | L2                                   c
c               ~~~~~~~~\V1    |                                      c
c       p+q (incoming)    \    |                                      c
c                        L3 \  |                                      c
c                             \|              d^I                     c
c                              ------->------                         c
c                             V2     p                                c
c                                                                     c
c       General form of the vertex                                    c
c                                                                     c
c       V = V_tree + i(F_L P_R + F_R P_R)                             c
c                                                                     c
c      Momentum arguments in formfactors: p^2=q^2=(p+q)^2=0           c
c      Other arguments:                                               c
c       form:      complex output array containing formfactor values  c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sdd_vert0_uwh(i,j,k,form)
c     Up quark-W-charged higgs in loop
      implicit double precision (a-h,o-z)
      double complex form(2),ckm,tmp
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,3
         do m=1,2
            tmp = e2/st2/2/sq2*am(k,m)*zh(1,m)*dconjg(ckm(j,l))*ckm(i,l)
     $           * cp1(um(l),wm,cm(m))
            form(1) =  form(1) - yd(i)*tmp
            form(2) =  form(2) + yd(j)*tmp
         end do
      end do
      return
      end

      subroutine sdd_vert0_uuw(i,j,k,form)
c     Up quark-up-quark-W in loop
      implicit double precision (a-h,o-z)
      double complex form(2),ckm,tmp
      common/km_mat/ckm(3,3)
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,3
         tmp = e2/st2/sq2*zr(2,k)*um(l)*yu(l)*dconjg(ckm(j,l))*ckm(i,l)
     $        *cp12(um(l),wm)
         form(1) =  form(1) + dm(i)*tmp
         form(2) =  form(2) + dm(j)*tmp
      end do
      return
      end

      subroutine sdd_vert0_wwu(i,j,k,form)
c     Up quark-W-W in loop
      implicit double precision (a-h,o-z)
      double complex form(2),ckm,tmp
      common/km_mat/ckm(3,3)
      common/fmass/em(3),um(3),dm(3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,3
         tmp = e2*e2/st2/st2/2*cr(k)*dconjg(ckm(j,l))*ckm(i,l)
     $        * cp11(wm,um(l))
         form(1) =  form(1) + dm(i)*tmp
         form(2) =  form(2) + dm(j)*tmp
      end do
      return
      end

      subroutine sdd_vert0_uuh(i,j,k,form)
c     Up quark-up-quark-charged higgs in loop
      implicit double precision (a-h,o-z)
      double complex form(2),ckm,tmp
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,3
         do m=1,2
            tmp = zr(2,k)*zh(1,m)*zh(2,m)*yu(l)*yu(l)*dconjg(ckm(j,l))
     $           * ckm(i,l)/sq2*(cp1(um(l),um(l),cm(m)) + um(l)*um(l)
     $           * cp0(um(l),um(l),cm(m)))
            form(1) =  form(1) - yd(i)*tmp
            form(2) =  form(2) - yd(j)*tmp
         end do
      end do
      return
      end

      subroutine sdd_vert0_hhu(i,j,k,form)
c     Charged higgs-charged higgs-up quark in loop
      implicit double precision (a-h,o-z)
      double complex form(2),ckm,tmp
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/yukawa/yl(3),yu(3),yd(3)
      common/km_mat/ckm(3,3)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,3
         do m=1,2
            do n=1,2
               tmp = e/2/st*v_hhs(k,n,m)*yu(l)*um(l)*dconjg(ckm(j,l))
     $              * ckm(i,l)*cp0(um(l),cm(m),cm(n))
               form(1) = form(1) - yd(i)*zh(1,n)*zh(2,m)*tmp
               form(2) = form(2) - yd(j)*zh(1,m)*zh(2,n)*tmp
            end do
         end do
      end do
      return
      end

      subroutine sdd_vert0_ccu(i,j,k,form)
c     Up squark-chargino-chargino in loop
      implicit double precision (a-h,o-z)
      double complex form(2),a1,a2
      double complex v_ccs,vl_duc,vr_duc
      double complex zu,zd,zpos,zneg
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,2
         do m=1,2
            a1 = dconjg(v_ccs(m,l,k))
            a2 = v_ccs(l,m,k)
            do n=1,6
               form(1) = form(1) + e/st/sq2
     $              * dconjg(vr_duc(i,n,l))*vl_duc(j,n,m)
     $              * (a1*cp1(fcm(l),fcm(m),sum(n)) 
     $              + a2*fcm(l)*fcm(m)*cp0(fcm(l),fcm(m),sum(n)))
               form(2) = form(2) + e/st/sq2
     $              * dconjg(vl_duc(i,n,l))*vr_duc(j,n,m)
     $              * (a2*cp1(fcm(l),fcm(m),sum(n)) 
     $              + a1*fcm(l)*fcm(m)*cp0(fcm(l),fcm(m),sum(n)))
            end do
         end do
      end do
      return
      end

      subroutine sdd_vert0_uuc(i,j,k,form)
c     Chargino-up squark-up squark in loop
      implicit double precision (a-h,o-z)
      double complex form(2),cv
      double complex vl_duc,vr_duc,v_uus
      double complex zu,zd,zpos,zneg
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      do l=1,6
         do m=1,6
            cv = v_uus(l,m,k)
            do n=1,2
               form(1) = form(1) - cv*dconjg(vr_duc(i,m,n))
     $              * vl_duc(j,l,n)*fcm(n)*cp0(sum(l),sum(m),fcm(n))
               form(2) = form(2) - cv*dconjg(vl_duc(i,m,n))
     $              * vr_duc(j,l,n)*fcm(n)*cp0(sum(l),sum(m),fcm(n))
          end do
        end do
      end do
      return
      end

      subroutine sdd_vert0_nnd(i,j,k,form)
c     Down squark-neutralino-neutralino in loop
      implicit double precision (a-h,o-z)
      double complex form(2),vv
      double complex v_nns,vl_ddn,vr_ddn
      double complex zu,zd,zn
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,4
         do m=1,4
            vv = v_nns(l,m,k)
            do n=1,6
               form(1) = form(1) - e/2/sct
     $              * dconjg(vr_ddn(i,n,l))*vl_ddn(j,n,m)
     $              * (dconjg(vv)*cp1(fnm(l),fnm(m),sdm(n)) 
     $              + vv*fnm(l)*fnm(m)*cp0(fnm(l),fnm(m),sdm(n)))
               form(2) = form(2) - e/2/sct
     $              * dconjg(vl_ddn(i,n,l))*vr_ddn(j,n,m)
     $              * (vv*cp1(fnm(l),fnm(m),sdm(n)) 
     $              +dconjg(vv)*fnm(l)*fnm(m)*cp0(fnm(l),fnm(m),sdm(n)))
            end do
         end do
      end do
      return
      end

      subroutine sdd_vert0_ddn(i,j,k,form)
c     Neutralino-down squark- down squark in loop
      implicit double precision (a-h,o-z)
      double complex form(2),cv
      double complex vl_ddn,vr_ddn,v_dds
      double complex zu,zd,zn
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      do l=1,6
         do m=1,6
            cv = v_dds(l,m,k)
            do n=1,4
               form(1) = form(1) - cv*dconjg(vr_ddn(i,m,n))
     $              * vl_ddn(j,l,n)*fnm(n)*cp0(sdm(l),sdm(m),fnm(n))
               form(2) = form(2) - cv*dconjg(vl_ddn(i,m,n))
     $              * vr_ddn(j,l,n)*fnm(n)*cp0(sdm(l),sdm(m),fnm(n))
            end do
         end do
      end do
      return
      end

      subroutine sdd_vert0_ddg(i,j,k,form)
c     Gluino-down squark- down squark in loop
      implicit double precision (a-h,o-z)
      logical init_alpha_susy
      double complex form(2)
      double complex v_dds,tmp
      double complex zu,zd
      double complex gm2,gm3
      common/gmass/gm1,gm2,gm3
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      if (init_alpha_susy) call init_alpha_s_susy
      al = 8/3.d0*g3d*g3d
      do l=1,6
         do m=1,6
            tmp = al*v_dds(l,m,k)*gm1*cp0(sdm(l),sdm(m),gm1)
            form(1) = form(1) + tmp*dconjg(zd(i+3,m))*zd(j,l)
            form(2) = form(2) + tmp*dconjg(zd(i,m))*zd(j+3,l)
         end do
      end do
      return
      end

      subroutine sdd_vert0(i,j,k,form)
c     Full bare Sdd formfactors at vanishing external momenta
      implicit double precision (a-h,o-z)
      double complex form(2)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      do l=1,2
         form(l) = (0.d0,0.d0)
      end do
c     Higgs and gauge boson contributions
      if (ih.eq.1) then
         call sdd_vert0_uwh(i,j,k,form)
         call sdd_vert0_uuw(i,j,k,form)
         call sdd_vert0_wwu(i,j,k,form)
         call sdd_vert0_hhu(i,j,k,form)
         call sdd_vert0_uuh(i,j,k,form)
      end if
c     chargino contributions
      if (ic.eq.1) then
         call sdd_vert0_ccu(i,j,k,form)
         call sdd_vert0_uuc(i,j,k,form)
      end if
c     neutralino contributions
      if (in.eq.1) then
         call sdd_vert0_nnd(i,j,k,form)
         call sdd_vert0_ddn(i,j,k,form)
      end if
c     gluino contributions
      if (ig.eq.1) call sdd_vert0_ddg(i,j,k,form)
      do l=1,2
         form(l) = form(l)/16/pi/pi
      end do
      return
      end

      subroutine sdd_vert_eff(i,j,k,form)
c     Full flavour violating Sdd formfactors at vanishing external momenta
c     including reducible contributions
      implicit double precision (a-h,o-z)
      double complex form(2)
      double complex dsl_sig0,dsr_sig0
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/vev/v1,v2
      if (i.eq.j) stop 'sdd_vert_eff defined only for i<>j'
      call sdd_vert0(i,j,k,form)
      form(1) = form(1) - zr(1,k)/v1*dsl_sig0(i,j)
      form(2) = form(2) - zr(1,k)/v1*dsr_sig0(i,j)
      end

c      subroutine sdd_ren0(i,j,k,form)
c     Full renormalized Sdd formfactors at vanishing external momenta
c      implicit double precision (a-h,o-z)
c      double complex form(2)
c      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
c      common/frconst/dzll(3),dzre(3),dzlq(3),dzru(3),dzrd(3)
c      common/grconst/dza,dzb,dz2,dx
c      call sdd_vert0(i,j,k,form)
c      if (i.ne.j) return
c     Renormalization has to be done ...
c      form(1) = form(1) + ...
c      form(2) = form(2) + ...
c      stop 'Renormalization of flavour-diagonal Hdd vertex not included'
c      end

