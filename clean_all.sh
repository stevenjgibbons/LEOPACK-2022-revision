#!/bin/sh
#
for dir in \
	LEOPACK_linons1 \
	LEOPACK_linons2
do
    cd $dir
    make clean
    cd ..
done
#
