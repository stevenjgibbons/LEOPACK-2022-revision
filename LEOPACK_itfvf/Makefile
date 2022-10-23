#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        itfvf
#
itfvf:	itfvf.f
	$(FORTRAN) -o ./itfvf itfvf.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
