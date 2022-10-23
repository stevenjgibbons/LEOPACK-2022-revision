#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        msvip
#
msvip:	msvip.f
	$(FORTRAN) -o ./msvip msvip.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
