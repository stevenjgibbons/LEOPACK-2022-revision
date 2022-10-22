#!/bin/sh
#
for dir in \
	bin \
	lib \
	LEOPACK_iic2cicsc \
	LEOPACK_linons1 \
	LEOPACK_linons2 \
	LEOPACK_blscnlsc \
	LEOPACK_blscnlsic  \
	LEOPACK_blscnlsic_evecs \
	LEOPACK_cicibcdts2 \
	LEOPACK_cicubcdts2 \
	LEOPACK_cicmibcdts2 \
	LEOPACK_cicmubcdts2
do
    cd $dir
    pwd
    make clean
    cd ..
done
#
