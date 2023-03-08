#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        djiepgrf
#
djiepgrf:	djiepgrf.f
	$(FORTRAN) -o ./djiepgrf djiepgrf.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
