#
# Makefile for PROGRAMS directory.
# Comment out one of the two following lines
OPTIM=   -O3
# OPTIM=        -O3 -fallow-argument-mismatch
FORTRAN= gfortran
#
all:	\
        blscnlsic_evecs
#
blscnlsic_evecs:	blscnlsic_evecs.f
	$(FORTRAN) -o ./blscnlsic_evecs blscnlsic_evecs.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
