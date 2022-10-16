#!/bin/sh
#
for dir in \
	LEOPACK_linons1 \
	LEOPACK_linons2 \
	LEOPACK_blscnlsc \
	LEOPACK_blscnlsic
do
    cd $dir
    pwd
    make clean
    cd ..
done
#
