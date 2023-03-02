#
# Makefile for PROGRAMS directory.
# Comment out one of the two following lines:
OPTIM=   -O3
#  OPTIM=        -O3 -fallow-argument-mismatch
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
