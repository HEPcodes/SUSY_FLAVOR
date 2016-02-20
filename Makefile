all: sflav

# FOPT = -O0 -g -fno-automatic -fbounds-check -Wall
FOPT = -O -fno-automatic -Wall
F77 = gfortran
COMP = $(F77) $(FOPT)

obj/rombint.o: lib/rombint.f
	$(COMP) -c lib/rombint.f
	mv rombint.o obj
obj/bsg_nl.o: lib/bsg_nl.f
	$(COMP) -c lib/bsg_nl.f
	mv bsg_nl.o obj
obj/ddg_fun.o: lib/ddg_fun.f
	$(COMP) -c lib/ddg_fun.f
	mv ddg_fun.o obj
obj/dd_gamma.o: lib/dd_gamma.f
	$(COMP) -c lib/dd_gamma.f
	mv dd_gamma.o obj
obj/dd_gluon.o: lib/dd_gluon.f
	$(COMP) -c lib/dd_gluon.f
	mv dd_gluon.o obj
obj/dd_vv.o: lib/dd_vv.f
	$(COMP) -c lib/dd_vv.f
	mv dd_vv.o obj
obj/dd_ll.o: lib/dd_ll.f
	$(COMP) -c lib/dd_ll.f
	mv dd_ll.o obj
obj/dd_mix.o: lib/dd_mix.f
	$(COMP) -c lib/dd_mix.f
	mv dd_mix.o obj
obj/uu_mix.o: lib/uu_mix.f
	$(COMP) -c lib/uu_mix.f
	mv uu_mix.o obj
obj/sff_fun0.o: lib/sff_fun0.f
	$(COMP) -c lib/sff_fun0.f
	mv sff_fun0.o obj
obj/sdd_vert0.o: lib/sdd_vert0.f
	$(COMP) -c lib/sdd_vert0.f
	mv sdd_vert0.o obj
obj/pdd_vert0.o: lib/pdd_vert0.f
	$(COMP) -c lib/pdd_vert0.f
	mv pdd_vert0.o obj
obj/zdd_vert0.o: lib/zdd_vert0.f
	$(COMP) -c lib/zdd_vert0.f
	mv zdd_vert0.o obj
obj/d_self0.o: lib/d_self0.f
	$(COMP) -c lib/d_self0.f
	mv d_self0.o obj
obj/suu_vert0.o: lib/suu_vert0.f
	$(COMP) -c lib/suu_vert0.f
	mv suu_vert0.o obj
obj/puu_vert0.o: lib/puu_vert0.f
	$(COMP) -c lib/puu_vert0.f
	mv puu_vert0.o obj
obj/u_self0.o: lib/u_self0.f
	$(COMP) -c lib/u_self0.f
	mv u_self0.o obj
obj/phen_4q.o: lib/phen_4q.f
	$(COMP) -c lib/phen_4q.f
	mv phen_4q.o obj
obj/phen_2q.o: lib/phen_2q.f
	$(COMP) -c lib/phen_2q.f
	mv phen_2q.o obj
obj/mh_diag.o: lib/mh_diag.f
	$(COMP) -c lib/mh_diag.f
	mv mh_diag.o obj
#
obj/cd_fun.o: lib/cd_fun.f
	$(COMP) -c lib/cd_fun.f
	mv cd_fun.o obj
obj/b_fun.o: lib/b_fun.f
	$(COMP) -c lib/b_fun.f
	mv b_fun.o obj
obj/db_fun.o: lib/db_fun.f
	$(COMP) -c lib/db_fun.f
	mv db_fun.o obj
obj/c_fun.o: lib/c_fun.f
	$(COMP) -c lib/c_fun.f
	mv c_fun.o obj
obj/qcd_fun.o: lib/qcd_fun.f
	$(COMP) -c lib/qcd_fun.f
	mv qcd_fun.o obj
obj/eisch1.o: lib/eisch1.f
	$(COMP) -c lib/eisch1.f
	mv eisch1.o obj
obj/vg_def.o: lib/vg_def.f
	$(COMP) -c lib/vg_def.f
	mv vg_def.o obj
obj/vf_def.o: lib/vf_def.f
	$(COMP) -c lib/vf_def.f
	mv vf_def.o obj
obj/vh_def.o: lib/vh_def.f
	$(COMP) -c lib/vh_def.f
	mv vh_def.o obj
obj/mh_init.o: lib/mh_init.f
	$(COMP) -c lib/mh_init.f
	mv mh_init.o obj
#
obj/cdm_d.o: lib/cdm_d.f
	$(COMP) -c lib/cdm_d.f
	mv cdm_d.o obj
obj/cdm_u.o: lib/cdm_u.f
	$(COMP) -c lib/cdm_u.f
	mv cdm_u.o obj
obj/edm_d.o: lib/edm_d.f 
	$(COMP) -c lib/edm_d.f
	mv edm_d.o obj
obj/edm_u.o: lib/edm_u.f 
	$(COMP) -c lib/edm_u.f
	mv edm_u.o obj
obj/edm_l.o: lib/edm_l.f
	$(COMP) -c lib/edm_l.f
	mv edm_l.o obj
obj/vegas.o: lib/vegas.f
	$(COMP) -c lib/vegas.f
	mv vegas.o obj
obj/cdm_g.o: lib/cdm_g.f
	$(COMP) -c lib/cdm_g.f
	mv cdm_g.o obj
obj/edm_n.o: lib/edm_n.f
	$(COMP) -c lib/edm_n.f
	mv edm_n.o obj
#
obj/sflav_io.o: lib/sflav_io.f
	$(COMP) -c lib/sflav_io.f
	mv sflav_io.o obj


libfcnc.a: obj/bsg_nl.o obj/ddg_fun.o obj/dd_gamma.o obj/dd_gluon.o	\
	obj/dd_vv.o obj/dd_ll.o obj/dd_mix.o obj/uu_mix.o		\
	obj/sff_fun0.o obj/zdd_vert0.o obj/sdd_vert0.o obj/pdd_vert0.o	\
	obj/d_self0.o obj/suu_vert0.o obj/puu_vert0.o obj/phen_4q.o	\
	obj/phen_2q.o obj/mh_diag.o obj/cd_fun.o obj/b_fun.o		\
	obj/db_fun.o obj/c_fun.o obj/qcd_fun.o obj/vg_def.o		\
	obj/vf_def.o obj/vh_def.o obj/mh_init.o obj/eisch1.o		\
	obj/rombint.o obj/cdm_d.o obj/cdm_u.o obj/cdm_g.o obj/vegas.o	\
	obj/edm_d.o obj/edm_l.o obj/edm_n.o obj/edm_u.o obj/sflav_io.o

	ar crs libfcnc.a obj/bsg_nl.o obj/ddg_fun.o obj/dd_gamma.o	\
	obj/dd_gluon.o obj/dd_vv.o obj/dd_ll.o obj/dd_mix.o		\
	obj/uu_mix.o obj/sff_fun0.o obj/zdd_vert0.o obj/sdd_vert0.o	\
	obj/pdd_vert0.o obj/d_self0.o obj/suu_vert0.o obj/puu_vert0.o	\
	obj/phen_4q.o obj/phen_2q.o obj/mh_diag.o obj/cd_fun.o		\
	obj/b_fun.o obj/db_fun.o obj/c_fun.o obj/qcd_fun.o		\
	obj/vg_def.o obj/vf_def.o obj/vh_def.o obj/mh_init.o		\
	obj/eisch1.o obj/rombint.o obj/cdm_d.o obj/cdm_u.o obj/cdm_g.o	\
	obj/vegas.o obj/edm_d.o obj/edm_l.o obj/edm_n.o obj/edm_u.o	\
	obj/sflav_io.o
#
#	Useful but cause problems on BSD systems:
#	ar ts libfcnc.a
#

sflav:	susy_flavor.f libfcnc.a
	$(COMP) -o sflav susy_flavor.f -L. -lfcnc
	@echo "Executing SUSY_FLAVOR..."
	./sflav

clean:
	rm -rf *.o obj/*.o *.a *~ lib/*~ sflav

tar:	
	rm -rf obj/*
	tar zcvf susy_flavor.tgz susy_flavor* *sample* obj lib/*.f Makefile
	mv susy_flavor.tgz $(HOME)