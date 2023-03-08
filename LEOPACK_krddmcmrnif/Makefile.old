#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        krddmcmrnif
#
krddmcmrnif:	krddmcmrnif.f
	$(FORTRAN) -o ./krddmcmrnif krddmcmrnif.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
