Steve Gibbons
2022-11-19 (NGI)

We want to generate initial conditions for the Benchmark dynamo
study of Christensen et al. (2001)

https://doi.org/10.1016/S0031-9201(01)00275-8

Codes in this release of LEOPACK can perform simulations
for Case 0 (non-magnetic convection) and 
Case 1 (full dynamo with insulating inner core).

Step 1: compile the code dbpisvecf

Type 

make dbpisvecf

This should generate an executable ./dbpisvecf ...

Step2: Generate the input for the Case 0 convection calculation.

./dbpisvecf < case0_init.dbpisvecf.input

where the file case0_init.dbpisvecf.input specifies a solution vector with 40 grid nodes with
Chebyshev spacing and an aspect ratio of 0.35. Maximum spherical harmonic l of 32 for the
temperature and the vorticity and no magnetic field.
We impose a 4-fold symmetry and consider wavenumbers up to 20.

This should generate the files
 case0_init.xarr
 case0_init.ints
 case0_init.vecs

The Case0 calculation could then be carried out by typing

o2ubtctsc2 < case0_o2ubtctsc2.input

which will read in the newly generated files and perform the non-magnetic time-stepping.
Adjust the parameters of this input file as required.

Step3: Generate the input file for the Case 1 dynamo calculation.

./dbpisvecf < case1_init.dbpisvecf.input

This should generate the three files
case1_init.xarr case1_init.ints case1_init.vecs

The dynamo simulation can be executed with the program o2ubcdts2 using
o2ubcdts2 < case1_o2ubcdts2.input
where the parameters in the o2ubcdts2 input file are adjusted as necessary.
