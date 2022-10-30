# LEOPACK-2022-revision
A set of FORTRAN77 codes for convection and magnetic field generation in rotating spherical geometry

Steven J. Gibbons (NGI) 2022-10-16

Between November 2001 and February 2002, I was working at the Institutt for geologi at the
University of Oslo, Norway, completing a postdoc position funded through the University of
Leeds in the UK. In this time, I collected and documented a set of codes for thermal convection and 
magnetic field generation in rotating spherical geometry that I had written and/or adapted to
solve diverse problems in Magnetohydrodynamics. Some of the codes had been started during
an earlier postdoc position at the University of Exeter. The product was a tar file containing
source code, input files, sample outputs, and a 600 page user guide.
The codes were subsequently adapted and modified by co-workers - David Gubbins, Chris Davies, Ashley Willis -
for subsequent research projects.

I named the software LEOPACK - Leeds, Exeter, Oslo - for want of a better name.

In this repository, I will add the codes - largely as they stood in February 2002 - first checking
that they work and making small modifications as necessary. On the one hand, there may be limited 
interest in the codes given that they are written in Fortran 77 and given that major
advances have been made in Magnetohydrodynamics research since February 2002.
On the other hand, the codes still work and compile with minimal or no modifications
(this can rarely be taken for granted after 20 years!) and there are many programs
for calculating solutions to fundamental problems for which I am not aware of any
alternative open codes.

The codes are largely self-contained. 
They use the LAPACK and BLAS linear algebra libraries.
The only additional requirements are from those codes which make 
graphical displays as these rely on the PGPLOT library:
https://sites.astro.caltech.edu/~tjp/pgplot/
As of today's date, I can obtain and compile the PGPLOT software from this source
and it still works (tested on Ubuntu 20 Linux). (In the original programs you could choose
gif, png, and postscript as output. I currently only get the postscript output to work.
I have not yet understood why the other formats are not working.)

To compile the graphics programs, make sure that the line
PGPLOT_LIB= [location of libpgplot.a]
in the file gprograms/Makefile
is set correctly prior to running the script
compile_graphics_programs.sh
and make sure that the environment variable PGPLOT_FONT is set to the location
of the file grfont.dat when you run any of the graphics programs.

I will add directories as and when I can test and document them.
I want each code to have its own pdf document, rather than the 600 page document
that came with the original package. I found this rather difficult to navigate
and, in the years that followed, I soon lost track of which programs were present.
This README file will be adapted with every directory I add.

Directories:
============

linalg (added 2022-10-18: contains all the BLAS, LAPACK, and ARPACK routines needed)

subs (added 2022-10-18: contains all the FORTRAN77 subroutines needed for the codes)

gsubs (added 2022-10-18: contains all the graphics subroutines that call PGPLOT routines)

lib, bin (added 2022-10-18: made to contain the libraries and graphics code binaries on compilation)

(Programs for the linear onset of thermal convection)
=====================================================

LEOPACK_linons1 (added 2022-10-16: Linear Onset of Thermal Convection 1)

LEOPACK_linons2 (added 2022-10-16: Linear Onset of Thermal Convection 2)

LEOPACK_sbrlinons1 (added 2022-10-24: Solid Body Rotation Linear Onset of Thermal Convection 1)

LEOPACK_sbrlinonsd (added 2022-10-24: Solid Body Rotation Linear Onset of Thermal Convection Drift Rate Solve)

(Programs for boundary locked convection)
=========================================

LEOPACK_blscnlsc (added 2022-10-16: Boundary Locked Steady Convection Non-Linear Solution Calculate)

LEOPACK_blscnlsic (added 2022-10-16: Boundary Locked Steady Convection Non-Linear Solution and Instability Calculate)

LEOPACK_blscnlsic_evecs (added 2022-10-17: Boundary Locked Steady Convection Non-Linear Solution and Instability Calculate with Eigenvectors)

(Programs for the Kinematic Dynamo Problem)
===========================================

LEOPACK_djiepgrf (added 2022-10-30: Dudley James velocity Instability Eigenvalue Problem Growth Rate Find)

LEOPACK_kriepgrf (added 2022-10-30: Kumar Roberts velocity Instability Eigenvalue Problem Growth Rate Find)

(Programs for full dynamo problem)
==================================

LEOPACK_o2ubcdts2 (added 2022-10-23: Uniform Boundary Convective Dynamo Time Step Code 2)

LEOPACK_cicibcdts2 (added 2022-10-17: Conducting Inner Core Inhomogeneous Boundary Convective Dynamo Time-Step code)

LEOPACK_cicmubcdts2 (added 2022-10-18: Conducting Inner Core and Mantle Uniform Boundary Convective Dynamo Time-Step code)

LEOPACK_cicubcdts2 (added 2022-10-22: Conducting Inner Core Uniform Boundary Convective Dynamo Time-Step code)

LEOPACK_cicmibcdts2 (added 2022-10-18: Conducting Inner Core and Mantle Inhomogeneous Boundary Convective Dynamo Time-Step code)

(Programs for non-magnetic thermal convection)
==============================================

LEOPACK_o2ibtctsc2 (added 2022-10-23: Inhomog. Boundary Thermal Conv. Time Step Code 2 )

LEOPACK_o2ubtctsc2 (added 2022-10-23: Uniform Boundary Thermal Conv. Time Step Code 2 )

(Programs for manipulation of solution vectors)
===============================================

LEOPACK_itfvf (added 2022-10-23: Inhomogeneous Temperature Function Vector Form)

LEOPACK_msvip (added 2022-10-23: Multiple Solution Vector Ingestion Program)

LEOPACK_rsvfg (added 2022-10-23: Random Solution Vector File Generate)

LEOPACK_iic2cicsc (added 2022-10-22: Insulating Inner Core 2 Conducting Inner Core Solution Convert)

LEOPACK_svenspec (added 2022-10-24: Solution Vector Energy SPECtrum)

LEOPACK_cicm2ocdisplay (added 2022-10-18: Conducting Inner Core and Mantle 2 Outer Core Display)

LEOPACK_svpnsmap (added 2022-10-30: Solution Vector Perturbation and New Spatial Mesh Adaptation Program)

LEOPACK_cicsvpnsmap (added 2022-10-30: Conducting Inner Core Solution Vector Perturbation and New Spatial Mesh Adaptation Program)

LEOPACK_cicmsvpnsmap (added 2022-10-30: Conducting Inner Core and Mantle Solution Vector Perturbation and New Spatial Mesh Adaptation Program)

LEOPACK_mfcanal1 (added 2022-10-30: Magnetic Field Component Analysis version 1)


(Programs for displaying solutions)
===================================

GRAPHICS_arrows_const_r3 (added 2022-10-29: Constant radius plot with arrows)

GRAPHICS_arrows_z_eq_merid4 (added 2022-10-29: Constant z or meridian plot with arrows)

GRAPHICS_cic_arrows_z_eq_merid4 (added 2022-10-29: Conducting Inner Core Constant z or meridian plot with arrows)

GRAPHICS_cicm_arrows_z_eq_merid4 (added 2022-10-29: Conducting Inner Core and Mantle Constant z or meridian plot with arrows)

GRAPHICS_continent_arrows_const_r3 (added 2022-10-29: Constant radius plot with arrows with Continents)

GRAPHICS_continent_full_sphere_plot (added 2022-10-29: Full Sphere Plot with Continents)

GRAPHICS_cutout_sphere_plot (added 2022-10-29: Cut Out Sphere Plot)

GRAPHICS_full_sphere_plot (added 2022-10-29: Full Sphere Plot)

GRAPHICS_ps_2plot_z_eq_merid2 (added 2022-10-29: Side-by-side Constant z or meridian plots)

GRAPHICS_shc_sphere_plot (added 2022-10-29: Spherical Harmonic Coefficient Sphere Plot)
