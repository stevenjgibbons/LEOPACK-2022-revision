#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        svenspec
#
svenspec:	svenspec.f
	$(FORTRAN) -o ./svenspec svenspec.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
