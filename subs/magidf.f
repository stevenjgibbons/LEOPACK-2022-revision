C*********************************************************************
C subroutine MAGnetic Inclination and Declination Find ***************
C            ---      -               -           -    ***************
C Steve Gibbons Wed Mar 21 12:43:22 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Takes in the double precision values BRAD, BTHE and PHI and        C
C calculates inclination DINCL and declination DDECL.                C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     BRAD      : Radial component of magnetic field.                C
C     BTHE      : Theta component of magnetic field.                 C
C     BPHI      : Phi component of magnetic field.                   C
C                                                                    C
C     DINCL     : Inclination                                        C
C     DDECL     : Declination                                        C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MAGIDF( BRAD, BTHE, BPHI, DINCL, DDECL )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION BRAD, BTHE, BPHI, DINCL, DDECL
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION DX, DY, DZ, DF, DFAC2
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DX = BTHE*(-1.0d0)
      DY = BPHI
      DZ = BRAD
C     .
      DF = DSQRT( DX*DX + DY*DY + DZ*DZ )
      DFAC2 = DZ/DF
C     .
      DDECL = DATAN2( DY, DX )
      DINCL = DASIN( DFAC2 )
C     .
      RETURN
      END
C*********************************************************************

