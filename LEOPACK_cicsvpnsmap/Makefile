#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        cicsvpnsmap
#
cicsvpnsmap:	cicsvpnsmap.f
	$(FORTRAN) -o ./cicsvpnsmap cicsvpnsmap.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
