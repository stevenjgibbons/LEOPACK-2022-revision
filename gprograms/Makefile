#
# Makefile for directory GPROGRAMS
# This Makefile must be modified to give the
# location of the PGPLOT libraries in your
# file system.
#
# Set FCOMPL to be the FORTRAN 77 compiler used.
#
FCOMPL=gfortran
#
# Change the following flags if the current settings
# are unsatisfactory
#
FFLAGD=  -O3
#
# Set CPG to location of file libgcc_path.sh
# This should be within the PGPLOT directory tree.
#
## CPG=./libgcc_path.sh
## CPG=$(HOME)/ext_programs/PGPLOT/pgplot/cpg/libgcc_path.sh
#
# Next line is probably somewhat system dependent.
# The following works for Sun SPARC/Solaris
#
# LIBS=-L/usr/openwin/lib -lX11 `$(CPG)` -lgcc -lm
#
# Set PGPLOT_LIB to the location of libpgplot.a
# This should be within the PGPLOT directory tree.
#
PGPLOT_LIB= $(HOME)/ext_programs/PGPLOT/pgplot/libpgplot.a
LEOPACK_BIN= ../bin
LEOPACK_DIR= ..
#
# Now we need to set variables for the various libraries
# of subroutines compiled within the LEOPACK distribution.
# First set SGLIB to be the location of the standard
# subprograms from the SUBS directory.
#
SGLIB=  $(LEOPACK_DIR)/lib/subslib.a
#
# Enter for LALIB the location of the LAPACK library
# (this may be a machine specific library)
#
LALIB=  $(LEOPACK_DIR)/lib/lalib.a
#
# Enter for BLLIB the location of the BLAS library
# (this may be a machine specific library)
#
BLLIB=  $(LEOPACK_DIR)/lib/bllib.a
#
# Finally, set PGSLIB to the location of the library generated
# from the subprograms in the directory GSUBS.
#
PGSLIB= $(LEOPACK_DIR)/lib/pgslib.a
#
NEWLIBS= $(PGSLIB) $(SGLIB) $(LALIB) $(BLLIB)
LIBRARIES= $(NEWLIBS) $(PGPLOT_LIB) $(LIBS)
#
# Now include instructions for the various programs
#
arrows_const_r3:	arrows_const_r3.f $(NEWLIBS)
	$(FCOMPL) $(FFLAGD) -o $(LEOPACK_BIN)/arrows_const_r3 arrows_const_r3.f $(LIBRARIES)
#
arrows_z_eq_merid4:	arrows_z_eq_merid4.f $(NEWLIBS)
	$(FCOMPL) $(FFLAGD) -o $(LEOPACK_BIN)/arrows_z_eq_merid4 arrows_z_eq_merid4.f $(LIBRARIES)
#
cicm_arrows_z_eq_merid4:	cicm_arrows_z_eq_merid4.f $(NEWLIBS)
	$(FCOMPL) $(FFLAGD) -o $(LEOPACK_BIN)/cicm_arrows_z_eq_merid4 cicm_arrows_z_eq_merid4.f $(LIBRARIES)
#
cic_arrows_z_eq_merid4:	cic_arrows_z_eq_merid4.f $(NEWLIBS)
	$(FCOMPL) $(FFLAGD) -o $(LEOPACK_BIN)/cic_arrows_z_eq_merid4 cic_arrows_z_eq_merid4.f $(LIBRARIES)
#
continent_arrows_const_r3:	continent_arrows_const_r3.f $(NEWLIBS)
	$(FCOMPL) $(FFLAGD) -o $(LEOPACK_BIN)/continent_arrows_const_r3 continent_arrows_const_r3.f $(LIBRARIES)
#
continent_full_sphere_plot:	continent_full_sphere_plot.f $(NEWLIBS)
	$(FCOMPL) $(FFLAGD) -o $(LEOPACK_BIN)/continent_full_sphere_plot continent_full_sphere_plot.f $(LIBRARIES)
#
cutout_sphere_plot:	cutout_sphere_plot.f $(NEWLIBS)
	$(FCOMPL) $(FFLAGD) -o $(LEOPACK_BIN)/cutout_sphere_plot cutout_sphere_plot.f $(LIBRARIES)
#
full_sphere_plot:	full_sphere_plot.f $(NEWLIBS)
	$(FCOMPL) $(FFLAGD) -o $(LEOPACK_BIN)/full_sphere_plot full_sphere_plot.f $(LIBRARIES)
#
ps_2plot_z_eq_merid2:	ps_2plot_z_eq_merid2.f $(NEWLIBS)
	$(FCOMPL) $(FFLAGD) -o $(LEOPACK_BIN)/ps_2plot_z_eq_merid2 ps_2plot_z_eq_merid2.f $(LIBRARIES)
#
shc_sphere_plot:	shc_sphere_plot.f $(NEWLIBS)
	$(FCOMPL) $(FFLAGD) -o $(LEOPACK_BIN)/shc_sphere_plot shc_sphere_plot.f $(LIBRARIES)
#
