c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: SUU_VERT0.F
c     Released: 26: 6:1999(P.Ch.)
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains scalar-up q-up q vertex formfactors          c
c                                                                     c
c     Incoming scalar:      momentum (p+q)                            c
c     Outgoing antiquark :  momentum  p                               c
c     Outgoing quark :      momentum  q                               c
c                                                                     c
c                              V3            ___                      c
c                               ____________ u^J                      c
c                             /|      p (outgoing)                    c
c                        L1 /  |                                      c
c                         /    |                                      c
c              S_i      /      | L2                                   c
c               ~~~~~~~~\V1    |                                      c
c       p+q (incoming)    \    |                                      c
c                        L3 \  |                                      c
c                             \|_____________ u^K                     c
c                              V2     q (outgoing)                    c
c                                                                     c
c       General form of the vertex (G(mu) = gamma(mu),                c
c       G(5)= gamma(5)):                                              c
c                                                                     c
c       V = V_tree + i(F_1 - i F_2 G(5))                              c
c                                                                     c
c      Momentum arguments in formfactors:                             c
c      p=p^2=0    q=q^2=0    pq = 1/2 ((p + q)^2 - p^2 - q^2)=0       c
c      Other arguments:                                               c
c       form:      complex output array containing formfactor values  c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine suu_vert0_1(i,j,k,form)
c     Down squark-chargino-chargino in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),a1,a2,a3,b1,b2,b3
      double complex v_ccs,vl_udc,vr_udc
      double complex zu,zd,zpos,zneg
      double complex cz,co,ci
      common/num/cz,co,ci,zero,one
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,2
        do m=1,2
          a1 = v_ccs(m,l,k) + dconjg(v_ccs(l,m,k))
          b1 = ci*(- v_ccs(m,l,k) + dconjg(v_ccs(l,m,k)))
          do n=1,6
            a2 = dconjg(vr_udc(j,n,l) + vl_udc(j,n,l))
            b2 = ci*dconjg(vl_udc(j,n,l) - vr_udc(j,n,l))
            a3 = vl_udc(i,n,m) + vr_udc(i,n,m)
            b3 = ci*(vr_udc(i,n,m) - vl_udc(i,n,m))
            call sff_svert0(fcm(m),sdm(n),fcm(l),a1,b1,a2,b2,a3,b3,tmp)
            do kk=1,2
              form(kk) = form(kk) - e/sq2/8/st*tmp(kk)
            end do
          end do
        end do
      end do
      return
      end
 
      subroutine suu_vert0_2(i,j,k,form)
c     Chargino-down squark-down squark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),a2,a3,b2,b3,cv
      double complex vl_udc,vr_udc,v_dds
      double complex zu,zd,zpos,zneg
      double complex cz,co,ci
      common/num/cz,co,ci,zero,one
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      do l=1,6
        do m=1,6
          cv = v_dds(m,l,k)
          do n=1,2
            a3 = vl_udc(i,m,n) + vr_udc(i,m,n)
            b3 = ci*(- vl_udc(i,m,n) + vr_udc(i,m,n))
            a2 = dconjg(vr_udc(j,l,n) + vl_udc(j,l,n))
            b2 = ci*dconjg(vl_udc(j,l,n) - vr_udc(j,l,n))
            call fss_svert0(sdm(m),fcm(n),sdm(l),a2,b2,a3,b3,tmp)
            do kk=1,2
              form(kk) = form(kk) + cv/4*tmp(kk)
            end do
          end do
        end do
      end do
      return
      end
 
      subroutine suu_vert0_3(i,j,k,form)
c     Up squark-neutralino-neutralino in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),a1,a2,a3,b1,b2,b3
      double complex v_nns,vl_uun,vr_uun
      double complex zu,zd,zn
      double complex cz,co,ci
      common/num/cz,co,ci,zero,one
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,4
        do m=1,4
          a1 = v_nns(l,m,k) + dconjg(v_nns(m,l,k))
          b1 = ci*(- v_nns(l,m,k) + dconjg(v_nns(m,l,k)))
          do n=1,6
            a2 = dconjg(vl_uun(j,n,l) + vr_uun(j,n,l))
            b2 = ci*dconjg(vl_uun(j,n,l) - vr_uun(j,n,l))
            a3 = vl_uun(i,n,m) + vr_uun(i,n,m)
            b3 = ci*(vr_uun(i,n,m) - vl_uun(i,n,m))
            call sff_svert0(fnm(m),sum(n),fnm(l),a1,b1,a2,b2,a3,b3,tmp)
            do kk=1,2
              form(kk) = form(kk) + e/16/sct*tmp(kk)
            end do
          end do
        end do
      end do
      return
      end
 
      subroutine suu_vert0_4(i,j,k,form)
c     Neutralino-up squark- up squark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),a2,a3,b2,b3
      double complex vl_uun,vr_uun,v_uus
      double complex zu,zd,zn
      double complex cz,co,ci,cv
      common/num/cz,co,ci,zero,one
      common/neut/fnm(4),zn(4,4)
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      do l=1,6
        do m=1,6
          cv = v_uus(l,m,k)
          do n=1,4
            a3 = vl_uun(i,m,n) + vr_uun(i,m,n)
            b3 = ci*(vr_uun(i,m,n) - vl_uun(i,m,n))
            a2 = dconjg(vl_uun(j,l,n) + vr_uun(j,l,n))
            b2 = ci*dconjg(vl_uun(j,l,n) - vr_uun(j,l,n))
            call fss_svert0(sum(m),fnm(n),sum(l),a2,b2,a3,b3,tmp)
            do kk=1,2
              form(kk) = form(kk) + cv/4*tmp(kk)
            end do
          end do
        end do
      end do
      return
      end
 
      subroutine suu_vert0_5(i,j,k,form)
c     Gluino-up squark- up squark in loop
      implicit double precision (a-h,o-z)
      double complex tmp(2),form(2),a2,a3,b2,b3
      double complex v_uus
      double complex zu,zd
      double complex cz,co,ci
      double complex gm2,gm3
      common/num/cz,co,ci,zero,one
      common/gmass/gm1,gm2,gm3
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      al = alfas(zm)  
      do l=1,6
        do m=1,6
          a2 = - zu(j,l) + zu(j+3,l)
          b2 = - ci*(zu(j,l) + zu(j+3,l))
          a3 = dconjg(- zu(i,m) + zu(i+3,m))
          b3 = ci*dconjg(zu(i,m) + zu(i+3,m))
          call fss_svert0(sum(m),gm1,sum(l),a2,b2,a3,b3,tmp)
          do n=1,2
            form(n) = form(n) + 8*pi*al/3.d0*v_uus(l,m,k)*tmp(n)
          end do
        end do
      end do
      return
      end
 
 
      subroutine suu_svert0(i,j,k,form)
c     Full bare Suu formfactors at vanishing external momenta
      implicit double precision (a-h,o-z)
      double complex form(2)
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      do l=1,2
         form(l) = (0.d0,0.d0)
      end do
      call suu_vert0_1(i,j,k,form)
      call suu_vert0_2(i,j,k,form)
      call suu_vert0_3(i,j,k,form)
      call suu_vert0_4(i,j,k,form)
      call suu_vert0_5(i,j,k,form)
c     No charged Higgs here at present
c      call suu_vert0_6(i,j,k,form)
c      call suu_vert0_7(i,j,k,form)
      do l=1,2
         form(l) = form(l)/16/pi/pi
      end do
      return
      end
 
c      subroutine suu_ren0(i,j,k,form)
c     Full renormalized Suu formfactors
c      implicit double precision (a-h,o-z)
c      double complex form(2),tmp(2)
c      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
c      common/frconst/dzll(3),dzre(3),dzlq(3),dzru(3),dzrd(3)
c      common/grconst/dza,dzb,dz2,dx
c      do l=1,2
c         form(l) = (0.d0,0.d0)
c      end do
c      call suu_svert0(i,j,k,tmp)
c      do l=1,2
c         form(l) = form(l) + tmp(l)
c      end do
c     Renormalization has to be done ...
c      form(1) = form(1) + ...
c      form(2) = form(2) + ...
c      return
c      end
