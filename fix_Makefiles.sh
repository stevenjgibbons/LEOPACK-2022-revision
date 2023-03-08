#!/bin/bash
#
# Steve Gibbons 2023/03/02
# Script obtained from page
# https://forum.mmm.ucar.edu/threads/how-to-fix-error-rank-mismatch-between-actual-argument-at-1-and-actual-argument-at-2.11998/#post-27521
#

#IF statement for GNU compiler issue
export GCC_VERSION=$(/usr/bin/gcc -dumpfullversion | awk '{print$1}')
export GFORTRAN_VERSION=$(/usr/bin/gfortran -dumpfullversion | awk '{print$1}')
export GPLUSPLUS_VERSION=$(/usr/bin/g++ -dumpfullversion | awk '{print$1}')

export GCC_VERSION_MAJOR_VERSION=$(echo $GCC_VERSION | awk -F. '{print $1}')
export GFORTRAN_VERSION_MAJOR_VERSION=$(echo $GFORTRAN_VERSION | awk -F. '{print $1}')
export GPLUSPLUS_VERSION_MAJOR_VERSION=$(echo $GPLUSPLUS_VERSION | awk -F. '{print $1}')

export version_10="10"

echo "##########################################"
if [ $GCC_VERSION_MAJOR_VERSION -ge $version_10 ] || [ $GFORTRAN_VERSION_MAJOR_VERSION -ge $version_10 ] || [ $GPLUSPLUS_VERSION_MAJOR_VERSION -ge $version_10 ]
then
  export fallow_argument=-fallow-argument-mismatch
  export boz_argument=-fallow-invalid-boz
  cp LEOPACK_blscnlsic_evecs/Makefile.new LEOPACK_blscnlsic_evecs/Makefile
  cp LEOPACK_blscnlsic/Makefile.new LEOPACK_blscnlsic/Makefile
  cp LEOPACK_djiepgrf/Makefile.new LEOPACK_djiepgrf/Makefile
  cp LEOPACK_krcmrnif/Makefile.new LEOPACK_krcmrnif/Makefile
  cp LEOPACK_krddmcmrnif/Makefile.new LEOPACK_krddmcmrnif/Makefile
  cp LEOPACK_kriepgrf/Makefile.new LEOPACK_kriepgrf/Makefile
  cp LEOPACK_krssgeps/Makefile.new LEOPACK_krssgeps/Makefile
  cp LEOPACK_linons1/Makefile.new LEOPACK_linons1/Makefile
  cp LEOPACK_linons2/Makefile.new LEOPACK_linons2/Makefile
  cp LEOPACK_sbrlinons1/Makefile.new LEOPACK_sbrlinons1/Makefile
  cp LEOPACK_sbrlinonsd/Makefile.new LEOPACK_sbrlinonsd/Makefile
  cp linalg/Makefile.new linalg/Makefile
else
  export fallow_argument=
  export boz_argument=
  cp LEOPACK_blscnlsic_evecs/Makefile.old LEOPACK_blscnlsic_evecs/Makefile
  cp LEOPACK_blscnlsic/Makefile.old LEOPACK_blscnlsic/Makefile
  cp LEOPACK_djiepgrf/Makefile.old LEOPACK_djiepgrf/Makefile
  cp LEOPACK_krcmrnif/Makefile.old LEOPACK_krcmrnif/Makefile
  cp LEOPACK_krddmcmrnif/Makefile.old LEOPACK_krddmcmrnif/Makefile
  cp LEOPACK_kriepgrf/Makefile.old LEOPACK_kriepgrf/Makefile
  cp LEOPACK_krssgeps/Makefile.old LEOPACK_krssgeps/Makefile
  cp LEOPACK_linons1/Makefile.old LEOPACK_linons1/Makefile
  cp LEOPACK_linons2/Makefile.old LEOPACK_linons2/Makefile
  cp LEOPACK_sbrlinons1/Makefile.old LEOPACK_sbrlinons1/Makefile
  cp LEOPACK_sbrlinonsd/Makefile.old LEOPACK_sbrlinonsd/Makefile
  cp linalg/Makefile.old linalg/Makefile
fi


export FFLAGS=$fallow_argument
export FCFLAGS=$fallow_argument


echo "FFLAGS = $FFLAGS"
echo "FCFLAGS = $FCFLAGS"
echo "##########################################"
