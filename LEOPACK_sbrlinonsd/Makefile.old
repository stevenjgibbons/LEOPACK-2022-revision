#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        sbrlinonsd
#
sbrlinonsd:	sbrlinonsd.f
	$(FORTRAN) -o ./sbrlinonsd sbrlinonsd.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
