#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        blscnlsc
#
blscnlsc:	blscnlsc.f
	$(FORTRAN) -o ./blscnlsc blscnlsc.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
