Steve Gibbons
2022-11-19 (NGI)

The file lower_mantle_coefs.txt was provided by Guy Masters:
a seismic anomaly map for the lowermost mantle which we want
to use as a thermal boundary condition.
To prepare it for files to use with 
the various codes, we have two small and self-contained FORTRAN
codes to make the set of files needed.

Step1: Convert Guy's file to the LEOPACK format spherical harmonic coefficients.

make gmshcc  
./gmshcc

This should generate a file called "output_coeffs.txt"

Step2: Convert the coefficients to the format needed for the program
itfvf which writes out the temperature solution vector which enforces
this boundary condition.

First, compile the program shc2itfvf

Type

make shc2itfvf

This should result in an executable file shc2itfvf.

Then, execute the script run_step2.sh
This should generate the files swave_temp.ints swave_temp.xarr and swave_temp.vecs
for the geometry we specify in the run_step2.sh script.

Step3: (optional) Display our heat-flux pattern at the outer boundary
by running the graphics program continent_arrows_const_r3 with
the input file continent_arrows_const_r3.input
