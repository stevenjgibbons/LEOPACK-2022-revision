#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        o2ubtctsc2
#
o2ubtctsc2:	o2ubtctsc2.f
	$(FORTRAN) -o ./o2ubtctsc2 o2ubtctsc2.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
