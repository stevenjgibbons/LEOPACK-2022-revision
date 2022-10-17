#!/bin/sh
#
for dir in \
	LEOPACK_linons1 \
	LEOPACK_linons2 \
	LEOPACK_blscnlsc \
	LEOPACK_blscnlsic  \
	LEOPACK_blscnlsic_evecs \
	LEOPACK_cicibcdts2
do
    cd $dir
    pwd
    make clean
    cd ..
done
#
