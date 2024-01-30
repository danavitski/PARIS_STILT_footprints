#
#  GNU-Makefile for Stochastic Time-Inverted Lagrangian Transport (STILT) particle dispersion model
#  Produces executables called 'hymodelc', 'getsurf', 'getfield' or 'chk_file'
#
#  The behaviour is affected by two definitions during invocation,
#     1) SITE=Jena|... selects the compiler
#     2) STYLE=opt|... selects the set of compiler options ("opt" for optimizing, other or none for
#        debug)
#  Invocation example: make SITE=Jena STYLE=opt
#
#  $Id: Makefile,v 1.34 2009-02-16 12:58:12 gerbig Exp $
#

###################################
ifeq ($(SITE),Jena)
#---- 1) Intel Fortran compiler with hardware as at the MPG Jena Linux cluster
#FC     = ifort
#ifeq ($(STYLE),opt)
#FFLAGS = -assume byterecl -fpe0 -ipo -O3 -no-prec-div -static -xW
#else
#FFLAGS = -assume byterecl -fpe0 -warn all,nodec,interfaces,unused -ftrapuv -check all -traceback \
#   -O0 -fp-stack-check -g -xW
#endif
#
#---- 2) Portland Group Fortran compiler 7.1 or higher (MPG Jena Linux cluster)
FC     = pgf95
ifeq ($(STYLE),opt)
FFLAGS = -fastsse -O4 -Mipa=fast -Ktrap=fp
else
#FFLAGS = -g -Mx,0,0x1000 -Mbounds -Mextend #use this on AER Linux cluster beehive
FFLAGS = -g -Mbounds -Ktrap=fp
#  Floating point exception handling is done by core dump (allow for sufficient
#  core dump size via "ulimit -c" or "limit coredumpsize") and calling the PGI
#  debugger: "pgdbg -core core <exefile>".
#
endif
else
#---- 3) Portland Group Fortran compiler 5.2 on Harvard DEAS Grid
FC     = pgf90
FFLAGS = -O -L/deas/vendor/depot/pgi-5.2/liblf -Mx,0,0x1000 -Mbounds -Mextend -Mlfs
#
#---- 4) Intel Fortran compiler 4.2 Itanium (compiles and runs, but results not verified yet)
#FC     = ifort
#FFLAGS = -align dcommons -align sequence -assume byterecl -fast -fno-alias -fno-fnalias -fpe0 -ftz -traceback
#
#---- 5) NEC SX-6 cross-compiler (compiles and runs, but results not verified yet)
#FC     = sxf90
#FFLAGS = -sx6 -Wf'-pvctl nomsg'
#
# IBM xlf compiler on GNU Linux, PowerPC architecture (AER humboldt cluster)
#FC=xlf90
#FFLAGS=-O -qsuffix=f=f90 -qmaxmem=-1
endif
###################################


#
#  Files needed to compile 'hymodelc'
#

OBJECTS1 = adviec.o advint.o advint2D.o advintWIND.o advintZLOC.o advintZML.o  \
           advmet.o advmetGRELL.o advpnt.o advrng.o advsrt.o basg2m.o basm2g.o \
           cgrell.o cgrnll.o cgszll.o cgszxy.o cll2xy.o condsk.o conini.o      \
           conset.o consum.o conzro.o cpolll.o cxy2ll.o depdry.o depelm.o      \
           deprad.o depset.o depsus.o emsdat.o emsgrd.o emspnt.o emsset.o      \
           gasdev.o gbl2ll.o gbl2xy.o gbldxy.o gblset.o geo_ll.o hymodelc.o    \
           limmath.o ll_geo.o map_xe.o map_xy.o metgrd.o metini.o metinp.o     \
           metold.o metpos.o metset.o metslp.o metwnd.o mpstrt.o               \
           module_defconc.o module_defgrid.o module_defmeto.o module_defsize.o \
           module_defspot.o module_map_utils.o constants_module.o              \
           misc_definitions_module.o mprintf.o output.o pakinp.o pardsp.o      \
           prfcom.o prfecm.o prfprs.o prfsig.o prfter.o prfwrf.o proj_3d.o     \
           pufdel.o pufdsp.o pufpar.o ran3.o ext_zsg.o rsat.o rtsafe.o         \
           runset.o sfcinp.o sobstr.o stbanl.o stbhor.o stbsnd.o stcm1p.o      \
           stlmbr.o sunang.o sunave.o sunclr.o sunflx.o thes.o tm2day.o        \
           tm2min.o x_prod.o xe_xy.o xy_map.o

#
#  Files needed to compile 'getsurf'
#

OBJECTS2 = basg2m.o basm2g.o cg2cxy.o cgrnll.o cgrnxy.o cgszll.o cll2xy.o   \
           cpolll.o cpolxy.o cxy2ll.o gbl2ll.o gbldxy.o gblset.o geo_ll.o   \
           getsurf.o limmath.o ll_geo.o map_xe.o map_xy.o metini.o metold.o \
           metset.o module_defsize.o module_defgrid.o mpstrt.o pakinp.o     \
           module_map_utils.o constants_module.o misc_definitions_module.o  \
           mprintf.o proj_3d.o sobstr.o sound.o stcm1p.o stlmbr.o tm2day.o  \
           tm2min.o x_prod.o xe_xy.o xy_map.o


#
#  Files needed to compile 'getfield'
#

OBJECTS3 = basg2m.o basm2g.o cg2cxy.o cgrnll.o cgrnxy.o cgszll.o cll2xy.o    \
           cpolll.o cpolxy.o cxy2ll.o gbl2ll.o gbldxy.o gblset.o geo_ll.o    \
           getfield.o limmath.o ll_geo.o map_xe.o map_xy.o metini.o metold.o \
           metset.o module_defsize.o module_defgrid.o mpstrt.o pakinp.o      \
           proj_3d.o sobstr.o sound.o stcm1p.o stlmbr.o tm2day.o tm2min.o    \
           x_prod.o xe_xy.o xy_map.o \
           module_map_utils.o constants_module.o misc_definitions_module.o mprintf.o


#
#  Files needed to compile 'chk_file'
#

OBJECTS4 = basg2m.o basm2g.o cgrnll.o cgszll.o chk_file.o cll2xy.o cpolll.o \
           cxy2ll.o gbl2ll.o gbl2xy.o gbldxy.o gblset.o geo_ll.o limmath.o  \
           ll_geo.o map_xe.o map_xy.o metini.o metold.o metset.o            \
           module_defsize.o module_defgrid.o mpstrt.o proj_3d.o sobstr.o    \
           stcm1p.o stlmbr.o tm2day.o tm2min.o x_prod.o xe_xy.o xy_map.o    \
           module_map_utils.o constants_module.o misc_definitions_module.o mprintf.o


#
#  define implicit rule
#

.SUFFIXES:
.SUFFIXES: .f90 .o
VPATH = .:diagnostics
%.o: %.f90
	$(FC) $(FFLAGS) -c $<


#
#  define targets
#

.PHONY: clean

hymodelc: Makefile $(OBJECTS1)
	$(FC) $(FFLAGS) $(OBJECTS1) -o $@

getsurf: Makefile $(OBJECTS2)
	$(FC) $(FFLAGS) $(OBJECTS2) -o $@

getfield: Makefile $(OBJECTS3)
	$(FC) $(FFLAGS) $(OBJECTS3) -o $@

chk_file: Makefile $(OBJECTS4)
	$(FC) $(FFLAGS) $(OBJECTS4) -o $@

chk_data: Makefile chk_data.o
	$(FC) $(FFLAGS) chk_data.o -o $@
clean:
	rm -f *.o *.oo *.mod *.ipo hymodelc getsurf getfield chk_file chk_data


#
#  Dependencies
#

advintZLOC.o: module_defsize.o
advintZML.o: module_defsize.o
advmet.o: module_defmeto.o
advmetGRELL.o: module_defmeto.o
advpnt.o: module_defconc.o module_defgrid.o module_defmeto.o
advrng.o: module_defgrid.o module_defmeto.o
advsrt.o: module_defsize.o
cgszll.o: module_map_utils.o module_defgrid.o
cgszxy.o: module_map_utils.o module_defgrid.o
chk_file.o: module_defgrid.o module_defsize.o
cll2xy.o: module_map_utils.o module_defgrid.o
condsk.o: module_defconc.o
conini.o: module_defconc.o module_defgrid.o module_defspot.o
conset.o: module_defconc.o
consum.o: module_defconc.o
conzro.o: module_defconc.o
cxy2ll.o: module_map_utils.o module_defgrid.o
depdry.o: module_defconc.o module_defspot.o
depelm.o: module_defconc.o
deprad.o: module_defconc.o
depset.o: module_defconc.o
depsus.o: module_defconc.o module_defgrid.o module_defmeto.o
emsdat.o: module_defgrid.o module_defspot.o
emsgrd.o: module_defconc.o module_defgrid.o module_defspot.o
emspnt.o: module_defconc.o module_defspot.o
emsset.o: module_defconc.o
gbl2ll.o: module_defgrid.o
gbl2xy.o: module_defgrid.o
gbldxy.o: module_defgrid.o
gblset.o: module_defgrid.o
getfield.o: module_defgrid.o module_defsize.o
getsurf.o: module_defgrid.o module_defsize.o
hymodelc.o: module_defconc.o module_defgrid.o module_defmeto.o module_defspot.o
metgrd.o: module_defgrid.o
metini.o: module_map_utils.o module_defgrid.o
metinp.o: module_defgrid.o
metold.o: module_defgrid.o
metpos.o: module_defgrid.o
metset.o: module_defgrid.o
module_defconc.o: module_defsize.o
module_defgrid.o: module_map_utils.o module_defsize.o
module_defmeto.o: module_defsize.o
module_defspot.o: module_defsize.o
module_map_utils.o: constants_module.o misc_definitions_module.o mprintf.o
pardsp.o: module_defsize.o
prfcom.o: module_defgrid.o
pufdel.o: module_defsize.o
pufdsp.o: module_defsize.o
pufpar.o: module_defgrid.o
runset.o: module_defgrid.o module_defmeto.o module_defspot.o
sound.o: module_defgrid.o module_defsize.o
stbanl.o: module_defsize.o
stbhor.o: module_defsize.o
stcm1p.o: module_defgrid.o
