c     PROGRAMS FOR ONE-LOOP ON-SHELL CALCULATIONS IN THE MSSM
c     Authors: P.H.Chankowski, S.Pokorski, J.Rosiek
c     e-mail: rosiek@fuw.edu.pl
c             chank@padova.infn.it

c     FILENAME: U_SELF0.FOR
c     Released: 1: 4:1994 (P.Ch.)
c     Revised:  4:02:2011 (J.R.)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains expressions for up quarks self-energy          c
c     function at s=0 and its renormalization.                          c
c                                                                       c
c     The definition of the self energy as                              c
c     follows (arguments are i,j):                                      c
c                                                                       c
c       k      ____    k                                                c
c             |    |          =    -i (G(k)uv_sig + G(k)G(5)ua_sig      c
c       ~~~~~~|____|~~~~~~               + us_sig +     G(5)up_sig)     c
c      l_i             l_j                                              c
c                                                                       c
c                                                                       c
c     i and j are the flavors of incoming and outgoing quark            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Vector self-energy (proportional to k(mu)gamma(mu) = G(k)))       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function uv_sig0_4(i,j)
c     Down quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b1,ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      uv_sig0_4 = 0
      do l=1,3
         do k=1,2
            uv_sig0_4 = uv_sig0_4 + (yu(j)*yu(i)*zh(2,k)**2
     $           + (zh(1,k)*yd(l))**2)/2
     $           * ckm(l,i)*dconjg(ckm(l,j))*b1(0.d0,dm(l),cm(k))
         end do
      end do
      return
      end

      double complex function uv_sig0_7(i,j)
c     Neutralino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zn
      double complex vl_uun,vr_uun
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      uv_sig0_7 = 0
      do k=1,6
        do l=1,4
           uv_sig0_7 = uv_sig0_7 
     $          + (vl_uun(i,k,l)*dconjg(vl_uun(j,k,l))
     $          + vr_uun(i,k,l)*dconjg(vr_uun(j,k,l)))/2
     $          * b1(0.d0,fnm(l),sum(k))
        end do
      end do
      return
      end

      double complex function uv_sig0_8(i,j)
c     Chargino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zpos,zneg
      double complex vl_udc,vr_udc
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      uv_sig0_8 = 0
      do k=1,6
        do l=1,2
           uv_sig0_8 = uv_sig0_8 
     $          + (vl_udc(i,k,l)*dconjg(vl_udc(j,k,l))
     $          + vr_udc(i,k,l)*dconjg(vr_udc(j,k,l)))/2
     $          * b1(0.d0,fcm(l),sdm(k))
        end do
      end do
      return
      end

      double complex function uv_sig0_9(i,j)
c     gluino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu0,zd0
      double complex zu
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      al = 4/3.d0*g3u*g3u
      uv_sig0_9 = 0
      do k=1,6
         uv_sig0_9 = uv_sig0_9 + al*(zu(j,k)*dconjg(zu(i,k))
     $        + zu(j+3,k)*dconjg(zu(i+3,k)))*b1(0.d0,gm1,sum(k))
      end do
      return
      end

      double complex function uv_sig0(i,j)
c     Full bare up quark self-energy, vector part
      implicit double precision (a-h,o-z)
      double complex uv_sig0_4,uv_sig0_7,uv_sig0_8,uv_sig0_9
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      common/u_self_stat/iv,ia,is,ip
      external init_u_self
      uv_sig0 = (0.d0,0.d0)
      if (ih.eq.1) uv_sig0 = uv_sig0 + uv_sig0_4(i,j) 
      if (in.eq.1) uv_sig0 = uv_sig0 + uv_sig0_7(i,j) 
      if (ic.eq.1) uv_sig0 = uv_sig0 + uv_sig0_8(i,j) 
      if (ig.eq.1) uv_sig0 = uv_sig0 + uv_sig0_9(i,j) 
      uv_sig0 = iv*uv_sig0/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Axial self-energy (proportional to k(mu)gamma(mu)gamma(5))        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function ua_sig0_4(i,j)
c     Down quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b1,ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      ua_sig0_4 = 0
      do l=1,3
        do k=1,2
           ua_sig0_4 = ua_sig0_4 + (yu(i)*yu(j)*zh(2,k)**2
     $          - (zh(1,k)*yd(l))**2)/2
     $          * ckm(l,i)*dconjg(ckm(l,j))*b1(0.d0,dm(l),cm(k))
        end do
      end do
      return
      end

      double complex function ua_sig0_7(i,j)
c     Neutralino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zn
      double complex vl_uun,vr_uun
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      ua_sig0_7 = 0
      do k=1,6
        do l=1,4
           ua_sig0_7 = ua_sig0_7 
     $          - (vl_uun(i,k,l)*dconjg(vl_uun(j,k,l))
     $          - vr_uun(i,k,l)*dconjg(vr_uun(j,k,l)))/2
     $          * b1(0.d0,fnm(l),sum(k))
        end do
      end do
      return
      end

      double complex function ua_sig0_8(i,j)
c     Chargino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu,zd,zpos,zneg
      double complex vl_udc,vr_udc
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      ua_sig0_8 = 0
      do k=1,6
        do l=1,2
           ua_sig0_8 = ua_sig0_8 
     $          - (vl_udc(i,k,l)*dconjg(vl_udc(j,k,l))
     $          - vr_udc(i,k,l)*dconjg(vr_udc(j,k,l)))/2
     $          * b1(0.d0,fcm(l),sdm(k))
        end do
      end do
      return
      end

      double complex function ua_sig0_9(i,j)
c     gluino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b1
      double complex zu0,zd0
      double complex zu
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      al = 4/3.d0*g3u*g3u
      ua_sig0_9 = 0
      do k=1,6
         ua_sig0_9 = ua_sig0_9 - al*(zu(j,k)*dconjg(zu(i,k))
     $        - zu(j+3,k)*dconjg(zu(i+3,k)))*b1(0.d0,gm1,sum(k))
      end do
      return
      end

      double complex function ua_sig0(i,j)
c     Full bare up quark self-energy, axial part
      implicit double precision (a-h,o-z)
      double complex ua_sig0_4,ua_sig0_7,ua_sig0_8,ua_sig0_9
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      common/u_self_stat/iv,ia,is,ip
      external init_u_self
      ua_sig0 = (0.d0,0.d0)
      if (ih.eq.1) ua_sig0 = ua_sig0 + ua_sig0_4(i,j) 
      if (in.eq.1) ua_sig0 = ua_sig0 + ua_sig0_7(i,j) 
      if (ic.eq.1) ua_sig0 = ua_sig0 + ua_sig0_8(i,j) 
      if (ig.eq.1) ua_sig0 = ua_sig0 + ua_sig0_9(i,j) 
      ua_sig0 = ia*ua_sig0/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Scalar self-energy                                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function us_sig0_4(i,j)
c     Down quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b0,ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      us_sig0_4 = 0
      do l=1,3
        do k=1,2
           us_sig0_4 = us_sig0_4 + zh(2,k)*zh(1,k)*dm(l)*yd(l)
     $          *(yu(i) + yu(j))*ckm(l,i)*dconjg(ckm(l,j))/2
     $          * b0(0.d0,dm(l),cm(k))
        end do
      end do
      return
      end

      double complex function us_sig0_7(i,j)
c     Neutralino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zn
      double complex vl_uun,vr_uun
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      us_sig0_7 = 0
      do k=1,6
        do l=1,4
           us_sig0_7 = us_sig0_7 
     $          - (vl_uun(i,k,l)*dconjg(vr_uun(j,k,l))
     $          + vr_uun(i,k,l)*dconjg(vl_uun(j,k,l)))*fnm(l)/2
     $          * b0(0.d0,fnm(l),sum(k))
        end do
      end do
      return
      end

      double complex function us_sig0_8(i,j)
c     Chargino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zpos,zneg
      double complex vl_udc,vr_udc
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      us_sig0_8 = 0
      do k=1,6
        do l=1,2
           us_sig0_8 = us_sig0_8 
     $          - (vl_udc(i,k,l)*dconjg(vr_udc(j,k,l))
     $          + vr_udc(i,k,l)*dconjg(vl_udc(j,k,l)))*fcm(l)/2
     $          * b0(0.d0,fcm(l),sdm(k))
        end do
      end do
      return
      end

      double complex function us_sig0_9(i,j)
c     gluino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu0,zd0
      double complex zu
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      al = 8/3.d0*g3u*g3u
      us_sig0_9 = 0
      do k=1,6
         us_sig0_9 = us_sig0_9 + al*(zu(j,k)*dconjg(zu(i+3,k))
     $        + zu(j+3,k)*dconjg(zu(i,k)))*gm1/2*b0(0.d0,gm1,sum(k))
      end do
      return
      end

      double complex function us_sig0(i,j)
c     Full bare up quark self-energy, scalar part
      implicit double precision (a-h,o-z)
      double complex us_sig0_4,us_sig0_7,us_sig0_8,us_sig0_9
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      common/u_self_stat/iv,ia,is,ip
      external init_u_self
      us_sig0 = (0.d0,0.d0)
      if (ih.eq.1) us_sig0 = us_sig0 + us_sig0_4(i,j) 
      if (in.eq.1) us_sig0 = us_sig0 + us_sig0_7(i,j) 
      if (ic.eq.1) us_sig0 = us_sig0 + us_sig0_8(i,j) 
      if (ig.eq.1) us_sig0 = us_sig0 + us_sig0_9(i,j) 
      us_sig0 = is*us_sig0/16/pi/pi
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Pseudoscalar self-energy (proportional to G(5)                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function up_sig0_4(i,j)
c     Down quark + charged Higgs in loop
      implicit double precision (a-h,o-z)
      double complex b0,ckm
      double complex yl,yu,yd
      common/hmass/cm(2),rm(2),pm(2),zr(2,2),zh(2,2)
      common/fmass/em(3),um(3),dm(3)
      common/km_mat/ckm(3,3)
      common/yukawa/yl(3),yu(3),yd(3)
      up_sig0_4 = 0
      do l=1,3
        do k=1,2
           up_sig0_4 = up_sig0_4 - zh(2,k)*zh(1,k)*dm(l)*yd(l)
     $          * (yu(i) - yu(j))*ckm(l,i)*dconjg(ckm(l,j))/2
     $          * b0(0.d0,dm(l),cm(k))
        end do
      end do
      return
      end

      double complex function up_sig0_7(i,j)
c     Neutralino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zn
      double complex vl_uun,vr_uun
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/neut/fnm(4),zn(4,4)
      up_sig0_7 = 0
      do k=1,6
        do l=1,4
           up_sig0_7 = up_sig0_7 
     $          - (vl_uun(i,k,l)*dconjg(vr_uun(j,k,l))
     $          - vr_uun(i,k,l)*dconjg(vl_uun(j,k,l)))*fnm(l)/2
     $          * b0(0.d0,fnm(l),sum(k))
        end do
      end do
      return
      end

      double complex function up_sig0_8(i,j)
c     Chargino and down squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu,zd,zpos,zneg
      double complex vl_udc,vr_udc
      common/sqmass/sum(6),sdm(6),zu(6,6),zd(6,6)
      common/charg/fcm(2),zpos(2,2),zneg(2,2)
      up_sig0_8 = 0
      do k=1,6
         do l=1,2
            up_sig0_8 = up_sig0_8 
     $           - (vl_udc(i,k,l)*dconjg(vr_udc(j,k,l))
     $           - vr_udc(i,k,l)*dconjg(vl_udc(j,k,l)))*fcm(l)/2
     $           * b0(0.d0,fcm(l),sdm(k))
         end do
      end do
      return
      end

      double complex function up_sig0_9(i,j)
c     gluino and up squark in loop
      implicit double precision (a-h,o-z)
      double complex b0
      double complex zu0,zd0
      double complex zu
      double complex gm2,gm3
      logical init_alpha_susy
      common/alpha_s_susy/g3u,g3d,init_alpha_susy
      common/sqmass/sum(6),sdm(6),zu0(6,6),zd0(6,6)
      common/gmass/gm1,gm2,gm3
      if (init_alpha_susy) call init_alpha_s_susy
      al = 8/3.d0*g3u*g3u
      up_sig0_9 = 0
      do k=1,6
         up_sig0_9 = up_sig0_9 - al*(zu(j,k)*dconjg(zu(i+3,k))
     $        - zu(j+3,k)*dconjg(zu(i,k)))*gm1/2*b0(0.d0,gm1,sum(k))
      end do
      return
      end

      double complex function up_sig0(i,j)
c     Full bare up quark self-energy, pseudoscalar part
      implicit double precision (a-h,o-z)
      double complex up_sig0_4,up_sig0_7,up_sig0_8,up_sig0_9
      common/vpar/st,ct,st2,ct2,sct,sct2,e,e2,alpha,wm,wm2,zm,zm2,pi,sq2
      common/debug_4q/ih,ic,in,ing,ig
      common/u_self_stat/iv,ia,is,ip
      external init_u_self
      up_sig0 = (0.d0,0.d0)
      if (ih.eq.1) up_sig0 = up_sig0 + up_sig0_4(i,j) 
      if (in.eq.1) up_sig0 = up_sig0 + up_sig0_7(i,j) 
      if (ic.eq.1) up_sig0 = up_sig0 + up_sig0_8(i,j) 
      if (ig.eq.1) up_sig0 = up_sig0 + up_sig0_9(i,j) 
      up_sig0 = ip*up_sig0/16/pi/pi
      return
      end

      block data init_u_self
      implicit double precision (a-h,o-z)
      common/u_self_stat/iv,ia,is,ip
      data iv,ia,is,ip/4*1/
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Basis change: 1,gamma(5) -> P_L,P_R                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double complex function uvl_sig0(i,j)
      implicit double precision (a-h,o-z)
      double complex uv_sig0,ua_sig0
      uvl_sig0 = uv_sig0(i,j) - ua_sig0(i,j)
      return
      end

      double complex function uvr_sig0(i,j)
      implicit double precision (a-h,o-z)
      double complex uv_sig0,ua_sig0
      uvr_sig0 = uv_sig0(i,j) + ua_sig0(i,j)
      return
      end

      double complex function usl_sig0(i,j)
      implicit double precision (a-h,o-z)
      double complex us_sig0,up_sig0
      usl_sig0 = us_sig0(i,j) - up_sig0(i,j)
      return
      end

      double complex function usr_sig0(i,j)
      implicit double precision (a-h,o-z)
      double complex us_sig0,up_sig0
      usr_sig0 = us_sig0(i,j) + up_sig0(i,j)
      return
      end


