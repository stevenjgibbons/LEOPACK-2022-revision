#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        o2ubcdts2
#
o2ubcdts2:	o2ubcdts2.f
	$(FORTRAN) -o ./o2ubcdts2 o2ubcdts2.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
