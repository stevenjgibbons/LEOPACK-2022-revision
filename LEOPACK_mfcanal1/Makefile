#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        mfcanal1
#
mfcanal1:	mfcanal1.f
	$(FORTRAN) -o ./mfcanal1 mfcanal1.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
