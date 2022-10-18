#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        cicm2ocdisplay
#
cicm2ocdisplay:	cicm2ocdisplay.f
	$(FORTRAN) -o ./cicm2ocdisplay cicm2ocdisplay.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
