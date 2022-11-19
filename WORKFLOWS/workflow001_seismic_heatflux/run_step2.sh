#!/bin/sh
if test ! -x ./shc2itfvf
then
  echo No executable ./shc2itfvf found ...
  exit 1
fi
if test ! -r output_coeffs.txt
then
  echo No file output_coeffs.txt found ...
  exit 1
fi
cat << EOF > shc2itfvf.input
output_coeffs.txt
swave_itfvf_coeffs.txt
EOF
./shc2itfvf < shc2itfvf.input
#
# This should have generated a file called swave_itfvf_coeffs.txt
# with 120 lines. Abort if this file has not been made
#
if test ! -r swave_itfvf_coeffs.txt
then
  echo No file swave_itfvf_coeffs.txt found ...
  exit 1
fi
numlines=`wc swave_itfvf_coeffs.txt | awk '{print $1}'`
if [ $numlines -ne 120 ]
then
  echo We were expecting 120 lines for the file swave_itfvf_coeffs.txt
  exit 1
fi
#
# Now generate an input file for itfvf
#
cat << EOF > itfvf.input
* input file for itfvf
swave_temp                               : Output filename stem
swave_itfvf_coeffs.txt                   : File containing coeff.s
 40  2   4   0.53846154   1.5384615  1   : nr isp ifrmf ri ro ithebc
 0.0     1.0                             : epsi  epso
EOF
itfvf < itfvf.input
