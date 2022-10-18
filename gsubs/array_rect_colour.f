C*********************************************************************
C subroutine ARRAY RECT COLOUR **************************************
C            ----- ----- ------ **************************************
C Steve Gibbons Mon Jan 22 13:51:42 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Colours in between (undrawn) contours for an array in polar coords C
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
C                  NLEV = -1 --> 17 contour levels with Andy's RBG   C
C                  colour scheme.                                    C
C     NPHI      : Number of phi grid nodes.                          C
C     NTHE      : Number of theta grid nodes.                        C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Data array with dimensions (NTHE,NPHI)             C
C     CONTL     : Dim (NLEV). Work array (to store contour values)   C
C     HUEPOS    : Colour for positive values                         C
C     HUENEG    : Colour for negative values                         C
C                                                                    C
C huepos and hueneg can take values in the interval [0,360] degrees. C
C 0 --> blue, 120 --> red and 240 --> green.                         C
C                                                                    C
C     CS        : 0.0 for monochrome and up to 1.0 for colour.       C
C     SCAL      : In interval ( 0.0, 1.0] 1.0 gives darkest colours  C
C                                                                    C
C   HUEPOS, HUENEG, CS and SCAL are only used if NLEV.NE.-1          C
C                                                                    C
C     VALMIN    : Minimum value to be contoured.                     C
C     VALMAX    : Maximum value to be contoured.                     C
C                                                                    C
C     DELTA     : Difference between highest contour level and       C
C                 maximum array value, or minimum array value and    C
C                 lowest contour level.                              C
C                                                                    C
C*********************************************************************
      SUBROUTINE ARRAY_RECT_COLOUR( F, NPHI, NTHE, NLEV, CONTL,
     1                              HUEPOS, HUENEG, CS, SCAL, DELTA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHI, NTHE, NLEV
      REAL    F( NPHI, NTHE ), CONTL( NLEV ),
     1        HUEPOS, HUENEG, CS, SCAL, DELTA
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IPHI, ITHE, NLEV2, ILEV
      REAL    RRED( 17 ), RRED1,
     1        RBLUE( 17 ), RBLUE1,
     2        RGREEN( 17 ), RGREEN1
      REAL    FMIN, FMAX, ZERO, FLIMIT, FEXTRM
      PARAMETER ( ZERO = 0.0, FLIMIT = 1.0e8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NLEV.EQ.-1 ) THEN
        NLEV2 = 17
      ELSE
        NLEV2 = NLEV
      ENDIF
C
      IF ( DELTA.LT.ZERO ) THEN
        PRINT *,' Subroutine ARRAY_RECT_COLOUR.'
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
        DO IPHI = 1, NPHI
          IF ( F( IPHI, ITHE ).GT.FMAX .AND.
     1         F( IPHI, ITHE ).LT.FLIMIT ) FMAX = F( IPHI, ITHE )
          IF ( F( IPHI, ITHE ).LT.FMIN ) FMIN = F( IPHI, ITHE )
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
C
C Select the colour scheme
C
      IF ( NLEV.EQ.-1 ) THEN
C
        FEXTRM = ABS( FMIN )
        IF ( ABS( FMAX ).GT.FEXTRM ) FEXTRM = ABS( FMAX )
        FMIN = (-1.0)*FEXTRM
        FMAX = FEXTRM
        CALL RCLDR( NLEV2, FMIN, FMAX, CONTL )
C
        RRED(   2 ) = 0.619608
        RRED(   3 ) = 0.666667
        RRED(   4 ) = 0.749020
        RRED(   5 ) = 0.819608
        RRED(   6 ) = 0.894118
        RRED(   7 ) = 1.000000
        RRED(   8 ) = 1.000000
        RRED(   9 ) = 1.000000
        RRED(  10 ) = 0.749020
        RRED(  11 ) = 0.647059
        RRED(  12 ) = 0.564706
        RRED(  13 ) = 0.466667
        RRED(  14 ) = 0.384314
        RRED(  15 ) = 0.329412
        RRED(  16 ) = 0.200000
        RRED(  17 ) = 0.094118
C
        RGREEN(   2 ) =  0.039216
        RGREEN(   3 ) =  0.156863
        RGREEN(   4 ) =  0.247059
        RGREEN(   5 ) =  0.349020
        RGREEN(   6 ) =  0.447059
        RGREEN(   7 ) =  0.549020
        RGREEN(   8 ) =  0.647059
        RGREEN(   9 ) =  0.756863
        RGREEN(  10 ) =  0.749020
        RGREEN(  11 ) =  0.647059
        RGREEN(  12 ) =  0.564706
        RGREEN(  13 ) =  0.466667
        RGREEN(  14 ) =  0.384314
        RGREEN(  15 ) =  0.329412
        RGREEN(  16 ) =  0.200000
        RGREEN(  17 ) =  0.188235
C
        RBLUE(   2 ) =  0.000000
        RBLUE(   3 ) =  0.039216
        RBLUE(   4 ) =  0.039216
        RBLUE(   5 ) =  0.039216
        RBLUE(   6 ) =  0.039216
        RBLUE(   7 ) =  0.156863
        RBLUE(   8 ) =  0.313726
        RBLUE(   9 ) =  0.466667
        RBLUE(  10 ) =  1.000000
        RBLUE(  11 ) =  1.000000
        RBLUE(  12 ) =  1.000000
        RBLUE(  13 ) =  1.000000
        RBLUE(  14 ) =  0.909804
        RBLUE(  15 ) =  0.839216
        RBLUE(  16 ) =  0.776471
        RBLUE(  17 ) =  0.713726
C
        DO ILEV = 2, 17
          RRED1 = RRED( ILEV )
          RBLUE1 = RBLUE( ILEV )
          RGREEN1 = RGREEN( ILEV )
          CALL PGSCR( ILEV, RRED1, RGREEN1, RBLUE1 )
        ENDDO
      ELSE
        CALL HLS_CONT_SET( NLEV2, CONTL, HUEPOS, HUENEG, CS, SCAL )
      ENDIF
C
C Fill in contour intervals
C
      CALL CONT_FILL_RECT( NPHI, NTHE, F, NLEV2, CONTL )
C
      RETURN
      END
C*********************************************************************
