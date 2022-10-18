C*********************************************************************
C subroutine CONTour DRAW in RECT coords *****************************
C            ----    ----    ----        *****************************
C Steve Gibbons Tue Jan 23 15:40:44 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Draws contour lines in Cart. coordinates. The PGPLOT default of    C
C solid for positive values and dashed for negative values is used.  C
C                                                                    C
C Requires COMMON block / PARAMC / THETA1, THETA2, PHI1, PHI2        C
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
C     IW        : Width of lines.                                    C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Dim (NPHI,NTHE). Function to contour.              C
C     CONTL     : Dim (NLEV). Resulting contour levels.              C
C                                                                    C
C*********************************************************************
      SUBROUTINE CONT_DRAW_RECT( NPHI, NTHE, F, NLEV, CONTL, IW )
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      REAL    THETA1, THETA2, PHI1, PHI2
      COMMON  / PARAMC / THETA1, THETA2, PHI1, PHI2
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHI, NTHE, NLEV, IW
      REAL    F( NTHE, NTHE ), CONTL( NLEV )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ILEV, NC, IX1, IX2, IY1, IY2
      REAL    CLEV, TR( 6 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NLEV.LT.2 ) THEN
        PRINT *,' Subroutine CONT_DRAW_RECT.'
        PRINT *,' NLEV = ', NLEV
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      TR( 1 ) = PHI1
      TR( 2 ) = (PHI2 - PHI1)/REAL( NPHI - 1 )
      TR( 3 ) = 0.0
      TR( 4 ) = THETA1
      TR( 5 ) = 0.0
      TR( 6 ) = (THETA2 - THETA1)/REAL( NTHE - 1 )
C
      NC     = 1
      IX1    = 1
      IX2    = NPHI
      IY1    = 1
      IY2    = NTHE
      CALL PGSLW( IW )
      DO ILEV = 1, NLEV
        CLEV  = CONTL( ILEV )
C       .
        CALL PGCONT( F, NPHI, NTHE, IX1, IX2, IY1, IY2, CLEV, NC,
     1               TR )
      ENDDO
      CALL PGSLW( 1 )
C
      RETURN
      END
C*********************************************************************
