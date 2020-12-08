# Makefile created by mkmf $Id: mkmf,v 18.0.18.4.2.1.2.3 2013/08/21 15:04:06 Seth.Underwood Exp $ 

SRCROOT = ../src.acc

include ../nvhpc.mk


.DEFAULT:
	-echo $@ does not exist.
all: c_sw.x
Main.o: ${SRCROOT}/Main.F90 fv_arrays.o sw_core.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	${SRCROOT}/Main.F90
sw_core.o: ${SRCROOT}/sw_core.F90 fv_arrays.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	${SRCROOT}/sw_core.F90
fv_arrays.o: ${SRCROOT}/fv_arrays.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	${SRCROOT}/fv_arrays.F90
./sw_core.F90: ${SRCROOT}/sw_core.F90
	cp ${SRCROOT}/sw_core.F90 .
./fv_arrays.F90: ${SRCROOT}/fv_arrays.F90
	cp ${SRCROOT}/fv_arrays.F90 .
./Main.F90: ${SRCROOT}/Main.F90
	cp /lustre/f1/unswept/Zhi.Liang/tests/gfdl_fv_mapz/src/./Main.F90 .
SRC = ${SRCROOT}/fv_arrays.F90 ${SRCROOT}/sw_core.F90 ${SRCROOT}/Main.F90
OBJ = fv_arrays.o sw_core.o Main.o
OFF = ${SRCROOT}/fv_arrays.F90 ${SRCROOT}/sw_core.F90 ${SRCROOT}/Main.F90
clean: neat
	-rm -f .c_sw.x.cppdefs $(OBJ) c_sw.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
c_sw.x: $(OBJ) 
	$(LD) $(OBJ) -o c_sw.x  $(LDFLAGS)
