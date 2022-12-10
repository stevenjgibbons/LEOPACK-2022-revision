Steve Gibbons
2022-12-10 (NGI)

We want to reproduce numbers and figures from the study of Gibbons et al. (2007)
https://doi.org/10.1080/03091920701472550

Table 1 in that paper gives marginal Rayleigh numbers for the onset of thermal convection.

We have two input files for the program sbrlinonsd which generate
the numbers in the two columns to the right of Table 1 in this paper.

Type

sbrlinonsd < Table1_column1_sbrlinonsd.input

to generate the left results column (i.e. column 3 in the table)
and type

sbrlinonsd < Table1_column2_sbrlinonsd.input

to generate the right results column (i.e. column 4 in the table).

You can see that the values of NH, NR, M in the two input files
correspond to the numbers in the tables.
The outputs are in the files
Table1_column1.res and Table1_column2.res respectively.

The files Figure2_sbrlinonsd.input and Figure4_sbrlinonsd.input
are input files for the sbrlinonsd program to calculate the solutions
displayed in Figures 2 and 4 respectively.

Running

sbrlinonsd < Figure2_sbrlinonsd.input

will generate files

Figure2_sbrlinonsd_OUTPUT.run001.ints  
Figure2_sbrlinonsd_OUTPUT.run001.veci  
Figure2_sbrlinonsd_OUTPUT.run001.vecr  
Figure2_sbrlinonsd_OUTPUT.run001.xarr  
Figure2_sbrlinonsd_OUTPUT.run002.ints  
Figure2_sbrlinonsd_OUTPUT.run002.veci  
Figure2_sbrlinonsd_OUTPUT.run002.vecr  
Figure2_sbrlinonsd_OUTPUT.run002.xarr  
Figure2_sbrlinonsd_OUTPUT.run003.ints  
Figure2_sbrlinonsd_OUTPUT.run003.veci   
Figure2_sbrlinonsd_OUTPUT.run003.vecr  
Figure2_sbrlinonsd_OUTPUT.run003.xarr  
Figure2_sbrlinonsd_OUTPUT.run004.ints  
Figure2_sbrlinonsd_OUTPUT.run004.veci  
Figure2_sbrlinonsd_OUTPUT.run004.vecr  
Figure2_sbrlinonsd_OUTPUT.run004.xarr  

An input file for plotting the contours and arrows for the first panel of Figure 2 is

Figure2topleft_arrows_z_eq_merid4.input

The file Y22_coeffs.txt contains only a single
line specifying the spherical harmonic to have as the outer boundary heat flux.

We generate the vectors to contain this thermal boundary condition using

itfvf < itfvf.input

and this generates the files

Y22_flux_function.ints   
Y22_flux_function.vecs   
Y22_flux_function.xarr   



