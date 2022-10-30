#!/bin/sh
for program in \
	blscnlsc     \
	blscnlsic     \
	blscnlsic_evecs     \
	cicibcdts2     \
	cicm2ocdisplay     \
	cicmibcdts2     \
	cicmsvpnsmap     \
	cicmubcdts2     \
	cicsvpnsmap     \
	cicubcdts2     \
	djiepgrf     \
	iic2cicsc     \
	itfvf     \
	krcmrnif     \
	krddmcmrnif     \
	kriepgrf     \
	krssgeps     \
	linons1     \
	linons2     \
	mfcanal1     \
	msvip     \
	o2ibtctsc2     \
	o2ubcdts2     \
	o2ubtctsc2     \
	rsvfg     \
	sbrlinons1     \
	sbrlinonsd     \
	svenspec     \
	svpnsmap
do
	make $program
done
