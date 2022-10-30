#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
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
