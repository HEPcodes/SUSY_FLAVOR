c     SUSY_FLAVOR: program for calculating rare decays in the general MSSM
c     Maintainer: Janusz Rosiek (janusz.rosiek@fuw.edu.pl)
c     Program web page: http://www.fuw.edu.pl/susy_flavor

c     FILENAME: SUSY_FLAVOR_FILE.F
c     Released: 25:10:2013(J.R.)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Driver file for SUSY_FLAVOR library                             c
c     Example of MSSM parameter initialization: input and output read c
c     from file/wrote to file                                         c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program susy_flavor_header
      implicit double precision (a-h,o-z)

      call sflav_input(ilev,ierr) ! read parameters from susy_flavor.in
      if (ierr.ne.0) write(*,*) 'Error in parameter initialization!' 

      call set_resummation_level(ilev,ierr) ! resummation of chiral corrections
      if (ierr.ne.0) write(*,*)ierr,
     $     'Error in chiral corrections resummation!'

      call susy_flavor          ! main routine calculating physical observables

      call sflav_output(ilev,ierr) ! output written to susy_flavor.out

      end

