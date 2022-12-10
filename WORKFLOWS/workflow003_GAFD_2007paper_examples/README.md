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
