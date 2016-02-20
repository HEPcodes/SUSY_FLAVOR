c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor
 
c     FILENAME: SFF_FUN0.F
c     Released: 26: 6: 1999(J.R.)
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This file contains standard functions for the                   c 
c     scalar(pseudoscalar)-fermion-fermion vertex formfactor          c
c     in the approximation of vanishing external masses and momenta   c
c                                                                     c
c                              V3            _                        c
c                               ____________ f                        c
c                             /|      p (outgoing)                    c
c                        L1 /  |                                      c
c                         /    |                                      c
c              S        /      | L2                                   c
c               ~~~~~~~~\V1    |                                      c
c       p+q (incoming)    \    |                                      c
c                        L3 \  |                                      c
c                             \|_____________ f'                      c
c                              V2     q (outgoing)                    c
c                                                                     c
c     General form of the vertex                                      c
c                                                                     c
c     V = V_tree + i(F_1 - i F_2 gamma(5))                            c
c                                                                     c
c     Momentum arguments in formfactors: p^2=q^2=(p+q)^2=0            c
c                                                                     c
c     Other arguments:                                                c
c     d1,d2,d3:        masses of particles circulating in loop        c 
c                      on lines L1,L2,L3 respectively                 c
c     v_i,a_i,s_i,p_i: vector and axial parts of couplings            c
c                      in the vertices (complex in general)           c
c     form:      complex output array containing formfactor values    c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine sff_svert0(d1,d2,d3,s1,p1,s2,p2,s3,p3,form)
c     Scalar - fermion - fermion vertex
c     Fermion-fermion-scalar loop
c     Vertices in loop:
c       V1: i (s1 - i p1 G(5))
c       V2: i (s2 - i p2 G(5))
c       V3: i (s3 - i p3 G(5))
      implicit double precision (a-h,o-z)
      double complex s1,p1,s2,p2,s3,p3
      double complex form(2)
c     Scalar formfactor
      form(1) = - (s1*(s2*s3 - p2*p3) + p1*(p2*s3 + s2*p3))
     $    * cp1(d1,d2,d3)
     $    - (s1*(s2*s3 - p2*p3) - p1*(p2*s3 + s2*p3))
     $    * d1*d3*cp0(d1,d2,d3)
c     Pseudoscalar formfactor
      form(2) = - (s1*(s2*p3 + p2*s3) - p1*(s2*s3 - p2*p3))
     $     * cp1(d1,d2,d3)
     $     - (s1*(s2*p3 + p2*s3) + p1*(s2*s3 - p2*p3))
     $     * d1*d3*cp0(d1,d2,d3)
      return  
      end

      subroutine fss_svert0(d1,d2,d3,s2,p2,s3,p3,form)
c     Scalar - fermion - fermion vertex
c     Scalar-scalar-fermion loop
c     Vertices in loop:
c       V1: i a, a  - factorized constant
c       V2: i (s2 - i p2 G(5))
c       V3: i (s3 - i p3 G(5)) 
      implicit double precision (a-h,o-z)
      double complex s2,p2,s3,p3
      double complex form(2)
c     Scalar formfactor
      form(1) = - (s2*s3 - p2*p3)*d2*cp0(d1,d2,d3)
c     Pseudoscalar formfactor
      form(2) = - (s2*p3 + p2*s3)*d2*cp0(d1,d2,d3)
      return
      end

      subroutine vff_svert0(d1,d2,d3,s1,p1,v2,a2,form)
c     Scalar - fermion - fermion vertex
c     Fermion-fermion-vector boson loop
c     Vertices in loop:
c       V1: i (s1 - i p1 G(5))
c       V2: i A G(al)(v2 - a2 G(5))  
c       V3: i B G(be)(v2 - a2 G(5))
c     A, B - factorized constants
      implicit double precision (a-h,o-z)
      double complex s1,p1,v2,a2
      double complex form(2)
      common/dimreg/idflag
c     Scalar formfactor
      form(1) =  4*s1*(v2*v2 - a2*a2)
     $     * (cp1(d1,d2,d3) + d1*d3*cp0(d1,d2,d3) 
c     Add this line to get the DIMREG result:
     $     - idflag/2.d0)
c     Pseudoscalar formfactor
      form(2) = 4*p1*(v2*v2 - a2*a2)
     $     * (cp1(d1,d2,d3) - d1*d3*cp0(d1,d2,d3)  
c     Add this line to get the DIMREG result:
     $     - idflag/2.d0)
      return
      end

      subroutine fvv_svert0(d1,d2,d3,v2,a2,v3,a3,form)
c     Scalar - fermion - fermion vertex
c     Vector boson-vector boson-fermion loop
c     Vertices in loop:
c       V1: i S; S factorized 
c       V2: i G(al)(v2 - a2 G(5))
c       V3: i G(be)(v3 - a3 G(5))
      implicit double precision (a-h,o-z)
      double complex v2,a2,v3,a3
      double complex form(2)
c     Scalar formfactor
      form(1) = 4*(v2*v3 - a2*a3)*d2*cp0(d1,d2,d3)
c     Pseudoscalar formfactor
      form(2) = 0
      return
      end

      subroutine fsv_svert0(d1,d2,d3,v2,a2,s3,p3,form)
c     Scalar - fermion - fermion vertex
c     Scalar-vector boson-fermion loop
c     Vertices in loop:
c       V1: i c^v (p+q+k);  c^v - factorized constant
c       V2: i G(al)(v2 - a2 G(5))
c       V3: i (s3 - i p3 G(5))
      implicit double precision (a-h,o-z)
      double complex v2,a2,s3,p3
      double complex form(2)
      double complex cz,co,ci
      common/num/cz,co,ci,zero,one
c     Scalar formfactor
      form(1) = (v2*s3 - ci*a2*p3)*cp1(d1,d2,d3) 
c     Pseudoscalar formfactor
      form(2) = (v2*p3 + ci*a2*s3)*cp1(d1,d2,d3) 
      return
      end

      subroutine fvs_svert0(d1,d2,d3,s2,p2,v3,a3,form)
c     Scalar - fermion - fermion vertex
c     Vector boson-scalar-fermion loop
c     Vertices in loop:
c       V1: i c^v (p+q+k)  c^v - factorized constant
c       V2: i (s2 - i p2 G(5))
c       V3: i G(al)(v3 - a3 G(5))
      implicit double precision (a-h,o-z)
      double complex s2,p2,v3,a3
      double complex form(2)
      double complex cz,co,ci
      common/num/cz,co,ci,zero,one
c     Scalar formfactor
      form(1) = - (s2*v3 + ci*p2*a3)*cp1(d1,d2,d3) 
c     Pseudoscalar formfactor
      form(2) = - (p2*v3 - ci*s2*a3)*cp1(d1,d2,d3) 
      return
      end

