#
# Makefile for PROGRAMS directory.
# Comment out one of the following two lines
OPTIM=   -O3
# OPTIM=   -O3 -fallow-argument-mismatch
FORTRAN= gfortran
#
all:	\
        linons2
#
linons2:	linons2.f
	$(FORTRAN) -o ./linons2 linons2.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
