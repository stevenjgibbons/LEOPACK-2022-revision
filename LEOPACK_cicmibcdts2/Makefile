#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        cicmibcdts2
#
cicmibcdts2:	cicmibcdts2.f
	$(FORTRAN) -o ./cicmibcdts2 cicmibcdts2.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
