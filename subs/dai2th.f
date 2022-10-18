C*********************************************************************
C Double Precision Function DAI2TH                            ********
C Steve Gibbons Wed Mar 21 13:17:43 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Transform inclination to distance to VGP.                          C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DINCL     : Inclination                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DAI2TH( DINCL )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION DINCL, DAI2TH
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION PI, DA, DFAC, DZERO
      PARAMETER ( PI=3.14159265358979312D0, DZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( DINCL.EQ.DZERO .OR. DINCL.EQ.PI ) THEN     
        DA = 0.5d0*PI
      ELSE
        DFAC = 2.0d0/DTAN( DINCL )
        DA   = DATAN( DFAC )
        IF ( DA.LT.DZERO ) DA = PI + DA
      ENDIF
      DAI2TH = DA
C     .
      RETURN
      END
C*********************************************************************
