#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        iic2cicsc
#
iic2cicsc:	iic2cicsc.f
	$(FORTRAN) -o ./iic2cicsc iic2cicsc.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
