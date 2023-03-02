#
# Makefile for PROGRAMS directory.
# Comment out one of the following two lines:
OPTIM=   -O3
# OPTIM=        -O3 -fallow-argument-mismatch
FORTRAN= gfortran
#
all:	\
        sbrlinons1
#
sbrlinons1:	sbrlinons1.f
	$(FORTRAN) -o ./sbrlinons1 sbrlinons1.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
