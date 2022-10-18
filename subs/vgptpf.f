C*********************************************************************
C subroutine Virtual Geomagnetic Pole Theta and Phi Find *************
C            -       -           -    -         -   -    *************
C Steve Gibbons Wed Mar 21 13:08:10 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Calculates THEP and PHIP - the site of a virtual geomagnetic pole  C
C given inclination and declination DINCL and DDECL, and a site with C
C theta and phi ( THES and PHIS ). Note that theta is colatitude     C
C and not latitude. All angles in radians.                           C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     THEP      : Colatitude of virtual geomagnetic pole.            C
C     PHIP      : Longitude  of virtual geomagnetic pole.            C
C                                                                    C
C     THES      : Colatitude of observation site.                    C
C     PHIS      : Longitude  of observation site.                    C
C                                                                    C
C     DINCL     : Inclination                                        C
C     DDECL     : Declination                                        C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VGPTPF( THEP, PHIP, THES, PHIS, DINCL, DDECL )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION THEP, PHIP, THES, PHIS, DINCL, DDECL
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION DA, PI, DLOW, DCTP, SBETA, DONE, BETA,
     1                 DFAC1, DFAC2, DAI2TH
      EXTERNAL         DAI2TH
      PARAMETER      ( PI=3.14159265358979312D0, DLOW = 1.0d-5,
     1                 DONE = 1.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Transform inclination to distance
C
      DA   = DAI2TH( DINCL )
      DCTP = DCOS( THES )*DCOS( DA ) +
     1       DSIN( THES )*DSIN( DA )*DCOS( DDECL )
      THEP = DACOS( DCTP )
C
      IF ( THEP.GT.DLOW .AND. (PI-THEP).GT.DLOW ) THEN
       SBETA = DSIN( DDECL )*DSIN( DA )/DSIN( THEP )
       IF ( DABS( SBETA ).GT.DONE ) SBETA = SBETA/DABS( SBETA )
       BETA  = DASIN( SBETA )
       DFAC1 = DCOS( DA )
       DFAC2 = DCOS( THES )*DCOS( THEP )
       IF ( DFAC1 .GE. DFAC2 ) THEN
         PHIP = PHIS + BETA
       ELSE
         PHIP = PHIS + PI - BETA
       ENDIF
      ELSE
       PHIP = 0.0d0
      ENDIF
      DFAC2 = 2.0d0*PI
      IF ( PHIP.GT.DFAC2 ) PHIP = PHIP - DFAC2
C     .
      RETURN
      END
C*********************************************************************

