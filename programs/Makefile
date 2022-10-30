#
# Makefile for PROGRAMS directory.
#
OPTIM=   -O3
FORTRAN= gfortran
LEOPACK_BIN= ../bin
LEOPACK_DIR= ..
SGLIB=  $(LEOPACK_DIR)/lib/subslib.a
LALIB=  $(LEOPACK_DIR)/lib/lalib.a
BLLIB=  $(LEOPACK_DIR)/lib/bllib.a
ARLIB=  $(LEOPACK_DIR)/lib/arlib.a
#
all:	blscnlsc
#
blscnlsc:	blscnlsc.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/blscnlsc blscnlsc.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
blscnlsic:	blscnlsic.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/blscnlsic blscnlsic.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
blscnlsic_evecs:	blscnlsic_evecs.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/blscnlsic_evecs blscnlsic_evecs.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
cicibcdts2:	cicibcdts2.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/cicibcdts2 cicibcdts2.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
cicm2ocdisplay:	cicm2ocdisplay.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/cicm2ocdisplay cicm2ocdisplay.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
cicmibcdts2:	cicmibcdts2.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/cicmibcdts2 cicmibcdts2.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
cicmsvpnsmap:	cicmsvpnsmap.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/cicmsvpnsmap cicmsvpnsmap.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
cicmubcdts2:	cicmubcdts2.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/cicmubcdts2 cicmubcdts2.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
cicsvpnsmap:	cicsvpnsmap.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/cicsvpnsmap cicsvpnsmap.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
cicubcdts2:	cicubcdts2.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/cicubcdts2 cicubcdts2.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
djiepgrf:	djiepgrf.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/djiepgrf djiepgrf.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
iic2cicsc:	iic2cicsc.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/iic2cicsc iic2cicsc.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
itfvf:	itfvf.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/itfvf itfvf.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
krcmrnif:	krcmrnif.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/krcmrnif krcmrnif.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
krddmcmrnif:	krddmcmrnif.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/krddmcmrnif krddmcmrnif.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
kriepgrf:	kriepgrf.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/kriepgrf kriepgrf.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
krssgeps:	krssgeps.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/krssgeps krssgeps.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
linons1:	linons1.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/linons1 linons1.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
linons2:	linons2.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/linons2 linons2.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
mfcanal1:	mfcanal1.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/mfcanal1 mfcanal1.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
msvip:	msvip.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/msvip msvip.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
o2ibtctsc2:	o2ibtctsc2.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/o2ibtctsc2 o2ibtctsc2.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
o2ubcdts2:	o2ubcdts2.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/o2ubcdts2 o2ubcdts2.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
o2ubtctsc2:	o2ubtctsc2.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/o2ubtctsc2 o2ubtctsc2.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
rsvfg:	rsvfg.f $(SGLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/rsvfg rsvfg.f $(SGLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
sbrlinonsd:	sbrlinonsd.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/sbrlinonsd sbrlinonsd.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
sbrlinons1:	sbrlinons1.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/sbrlinons1 sbrlinons1.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
svenspec:	svenspec.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/svenspec svenspec.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
svpnsmap:	svpnsmap.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB)
	$(FORTRAN) -o $(LEOPACK_BIN)/svpnsmap svpnsmap.f $(SGLIB) $(ARLIB) $(LALIB) $(BLLIB) $(OPTIM)
#
