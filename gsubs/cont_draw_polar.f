C*********************************************************************
C subroutine CONTour DRAW in POLAR coords ****************************
C            ----    ----    -----        ****************************
C Steve Gibbons Fri Jan 19 14:28:44 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Draws contour lines in polar coordinates. The PGPLOT default of    C
C solid for positive values and dashed for negative values is used.  C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NRAD      : Number of radial grid nodes.                       C
C     NTHE      : Number of theta grid nodes.                        C
C     NLEV      : Number of contour levels.                          C
C     IW        : Width of lines.                                    C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Dim (NRAD,NTHE). Function to contour.              C
C     CONTL     : Dim (NLEV). Resulting contour levels.              C
C                                                                    C
C*********************************************************************
      SUBROUTINE CONT_DRAW_POLAR( NRAD, NTHE, F, NLEV, CONTL, IW )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NRAD, NTHE, NLEV, IW
      REAL    F( NRAD, NTHE ), CONTL( NLEV )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ILEV, NC, IX1, IX2, IY1, IY2
      REAL    CLEV
      EXTERNAL POLAR_SUB
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NLEV.LT.2 ) THEN
        PRINT *,' Subroutine CONT_DRAW_POLAR.'
        PRINT *,' NLEV = ', NLEV
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NC     = 1
      IX1    = 1
      IX2    = NRAD
      IY1    = 1
      IY2    = NTHE
      CALL PGSLW( IW )
      DO ILEV = 1, NLEV
        CLEV  = CONTL( ILEV )
C       .
        CALL PGCONX( F, NRAD, NTHE, IX1, IX2, IY1, IY2, CLEV, NC,
     1               POLAR_SUB )
      ENDDO
      CALL PGSLW( 1 )
C
      RETURN
      END
C*********************************************************************
