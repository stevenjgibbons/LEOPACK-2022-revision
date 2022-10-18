C*********************************************************************
C subroutine ARRAY POLAR DRAW ****************************************
C            ----- ----- ---- ****************************************
C Steve Gibbons Mon Jan 22 14:27:52 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Draws contour lines for an array in polar coords.                  C
C It automatically calculates appropriate contours etc.              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NLEV      : Number of contour levels.                          C
C                 (NLEV.EQ.-1 is alias for 17 for colouring)         C
C     NRAD      : Number of radial grid nodes.                       C
C     NTHE      : Number of theta grid nodes.                        C
C     IW        : Width of lines.                                    C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Data array with dimensions (NRAD,NTHE)             C
C     CONTL     : Dim (NLEV). Work array (to store contour values)   C
C     DELTA     : Difference between highest contour level and       C
C                 maximum array value, or minimum array value and    C
C                 lowest contour level.                              C
C                                                                    C
C*********************************************************************
      SUBROUTINE ARRAY_POLAR_DRAW( F, NRAD, NTHE, NLEV, CONTL, DELTA,
     1                             IW )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NRAD, NTHE, NLEV, IW
      REAL    F( NRAD, NTHE ), CONTL( NLEV ), DELTA
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRAD, ITHE, NLEV2
      REAL    FMIN, FMAX, ZERO, FLIMIT, FEXTRM
      PARAMETER ( ZERO = 0.0, FLIMIT = 1.0e8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NLEV2 = NLEV
      IF ( NLEV.EQ.-1 ) NLEV2 = 17
C
      IF ( DELTA.LT.ZERO ) THEN
        PRINT *,' Subroutine ARRAY_POLAR_DRAW .'
        PRINT *,' DELTA = ', DELTA
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Find FMIN and FMAX
C
      FMAX = F( 8, 8 )
      FMIN = FMAX
      DO ITHE = 1, NTHE
        DO IRAD = 1, NRAD
          IF ( F( IRAD, ITHE ).GT.FMAX .AND.
     1         F( IRAD, ITHE ).LT.FLIMIT ) FMAX = F( IRAD, ITHE )
          IF ( F( IRAD, ITHE ).LT.FMIN ) FMIN = F( IRAD, ITHE )
        ENDDO
      ENDDO
C
C Modify FMIN and FMAX
C
      FMIN = FMIN - DELTA
      FMAX = FMAX + DELTA
C
C Define contour levels
C
      CALL RCLDR( NLEV2, FMIN, FMAX, CONTL )
      IF ( NLEV.EQ.-1 ) THEN
        FEXTRM = ABS( FMIN )
        IF ( ABS( FMAX ).GT.FEXTRM ) FEXTRM = ABS( FMAX )
        FMIN = (-1.0)*FEXTRM
        FMAX = FEXTRM
        CALL RCLDR( NLEV2, FMIN, FMAX, CONTL )
      ENDIF
C
C Fill in contour intervals
C
      CALL CONT_DRAW_POLAR( NRAD, NTHE, F, NLEV2, CONTL, IW )
C
      RETURN
      END
C*********************************************************************
