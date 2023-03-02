#
# Makefile for PROGRAMS directory.
# Comment out one of the two following lines:
OPTIM=   -O3
# OPTIM=        -O3 -fallow-argument-mismatch
#
FORTRAN= gfortran
#
all:	\
        blscnlsic
#
blscnlsic:	blscnlsic.f
	$(FORTRAN) -o ./blscnlsic blscnlsic.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
