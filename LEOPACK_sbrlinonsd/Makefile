#
# Makefile for PROGRAMS directory.
# Comment out one of the following two lines:
OPTIM=   -O3
# OPTIM=        -O3 -fallow-argument-mismatch
FORTRAN= gfortran
#
all:	\
        sbrlinonsd
#
sbrlinonsd:	sbrlinonsd.f
	$(FORTRAN) -o ./sbrlinonsd sbrlinonsd.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
