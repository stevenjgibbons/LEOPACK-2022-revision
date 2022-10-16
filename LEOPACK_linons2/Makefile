#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
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
