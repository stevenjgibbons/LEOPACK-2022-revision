#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        cicibcdts2
#
cicibcdts2:	cicibcdts2.f
	$(FORTRAN) -o ./cicibcdts2 cicibcdts2.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
