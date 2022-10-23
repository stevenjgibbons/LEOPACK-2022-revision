#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        o2ibtctsc2
#
o2ibtctsc2:	o2ibtctsc2.f
	$(FORTRAN) -o ./o2ibtctsc2 o2ibtctsc2.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
