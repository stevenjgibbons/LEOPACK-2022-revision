C*********************************************************************
C subroutine CONTour FILL in RECTangular coords **********************
C            ----    ----    ----               **********************
C Steve Gibbons Tue Jan 23 14:40:58 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Requires COMMON block / PARAMC / THETA1, THETA2, PHI1, PHI2        C
C                                                                    C
C Fills (colours) in between the contour lines in polar geometry.    C
C The region between contour levels (ILEV-1) and ILEV is filled      C
C with the colour with code ILEV. (See for example my routine        C
C HLS_CONT_SET for definition of colours from contour levels).       C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NPHI      : Number of phi grid nodes.                          C
C     NTHE      : Number of theta grid nodes.                        C
C     NLEV      : Number of contour levels.                          C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Dim (NPHI,NTHE). Function to contour.              C
C     CONTL     : Dim (NLEV). Resulting contour levels.              C
C                                                                    C
C*********************************************************************
      SUBROUTINE CONT_FILL_RECT( NPHI, NTHE, F, NLEV, CONTL )
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      REAL    THETA1, THETA2, PHI1, PHI2
      COMMON  / PARAMC / THETA1, THETA2, PHI1, PHI2
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHI, NTHE, NLEV
      REAL    F( NPHI, NTHE ), CONTL( NLEV )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ILEV, IX1, IX2, IY1, IY2
      REAL    CLEV, CLEV1, TR( 6 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      TR( 1 ) = PHI1
      TR( 2 ) = (PHI2 - PHI1)/REAL( NPHI - 1 )
      TR( 3 ) = 0.0
      TR( 4 ) = THETA1
      TR( 5 ) = 0.0
      TR( 6 ) = (THETA2 - THETA1)/REAL( NTHE - 1 )
C
      IF ( NLEV.LT.2 ) THEN
        PRINT *,' Subroutine CONT_FILL_RECT.'
        PRINT *,' NLEV = ', NLEV
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IX1    = 1
      IX2    = NPHI
      IY1    = 1
      IY2    = NTHE
C
      DO ILEV = 2, NLEV
C
C Sets the colour for the interval between contours
C ILEV - 1 and ILEV.
C
        CALL PGSCI( ILEV )
C
        CLEV  = CONTL( ILEV-1 )
        CLEV1 = CONTL( ILEV )
        CALL PGCONF( F, NPHI, NTHE, IX1, IX2, IY1, IY2, CLEV,
     1                     CLEV1, TR )
      ENDDO
C
      RETURN
      END
C*********************************************************************
