#
# Makefile for PROGRAMS directory.
# Comment out one of the two following lines:
OPTIM=   -O3
#  OPTIM=        -O3 -fallow-argument-mismatch
FORTRAN= gfortran
#
all:	\
        krssgeps
#
krssgeps:	krssgeps.f
	$(FORTRAN) -o ./krssgeps krssgeps.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
