#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        cicmsvpnsmap
#
cicmsvpnsmap:	cicmsvpnsmap.f
	$(FORTRAN) -o ./cicmsvpnsmap cicmsvpnsmap.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
