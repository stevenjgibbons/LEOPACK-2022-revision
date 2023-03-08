#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        blscnlsic
#
blscnlsic:	blscnlsic.f
	$(FORTRAN) -o ./blscnlsic blscnlsic.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
