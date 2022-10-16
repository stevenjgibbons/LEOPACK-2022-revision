#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        linons1
#
linons1:	linons1.f
	$(FORTRAN) -o ./linons1 linons1.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
