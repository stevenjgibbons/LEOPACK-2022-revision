#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        svpnsmap
#
svpnsmap:	svpnsmap.f
	$(FORTRAN) -o ./svpnsmap svpnsmap.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
