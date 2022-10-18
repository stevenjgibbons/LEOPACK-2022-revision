#!/bin/sh
PGPLIB=`grep PGPLOT_LIB= gprograms/Makefile | awk '{print $2}'`
echo "Make sure that $PGPLIB exists"
cd linalg
make 
cd ../subs
make
cd ../gsubs
make
cd ../gprograms
sh ./compile_all.sh
