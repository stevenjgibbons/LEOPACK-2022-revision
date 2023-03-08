#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
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
