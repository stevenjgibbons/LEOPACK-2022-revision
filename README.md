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
and it still works (tested on Ubuntu 20 Linux).

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

LEOPACK_linons1 (added 2022-10-16: Linear Onset of Thermal Convection 1)

LEOPACK_linons2 (added 2022-10-16: Linear Onset of Thermal Convection 2)

LEOPACK_blscnlsc (added 2022-10-16: Boundary Locked Steady Convection Non-Linear Solution Calculate)

LEOPACK_blscnlsic (added 2022-10-16: Boundary Locked Steady Convection Non-Linear Solution and Instability Calculate)

LEOPACK_blscnlsic_evecs (added 2022-10-17: Boundary Locked Steady Convection Non-Linear Solution and Instability Calculate with Eigenvectors)

LEOPACK_cicibcdts2 (added 2022-10-17: Conducting Inner Core Inhomogeneous Boundary Convective Dynamo Time-Step code

LEOPACK_cicmubcdts2 (added 2022-10-18: Conducting Inner Core and Mantle Uniform Boundary Convective Dynamo Time-Step code
