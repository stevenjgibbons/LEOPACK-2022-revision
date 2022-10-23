#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
#
all:	\
        rsvfg
#
rsvfg:	rsvfg.f
	$(FORTRAN) -o ./rsvfg rsvfg.f $(OPTIM)
	\rm -f *.o
#
clean:
	sh clean.sh
#
