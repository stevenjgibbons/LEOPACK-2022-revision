#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        cicmubcdts2
#
cicmubcdts2:	cicmubcdts2.f
	$(FORTRAN) -o ./cicmubcdts2 cicmubcdts2.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
