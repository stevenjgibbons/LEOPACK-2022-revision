#
# Makefile for PROGRAMS directory.
# Comment out one of the two lines below
OPTIM=   -O3
# OPTIM=   -O3 -fallow-argument-mismatch
FORTRAN= gfortran
#
all:	\
        kriepgrf
#
kriepgrf:	kriepgrf.f
	$(FORTRAN) -o ./kriepgrf kriepgrf.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
