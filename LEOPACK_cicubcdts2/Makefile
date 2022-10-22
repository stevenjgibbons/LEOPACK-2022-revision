#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        cicubcdts2
#
cicubcdts2:	cicubcdts2.f
	$(FORTRAN) -o ./cicubcdts2 cicubcdts2.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
