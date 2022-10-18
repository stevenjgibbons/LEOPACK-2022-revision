
C*********************************************************************
C subroutine SIMPLE SPHERE PLOT **************************************
C            ------ ------ ---- **************************************
C Steve Gibbons Thu Nov  8 12:08:09 WET 2001                         C
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
C                                                                    C
C     NPHIRAD   : Number of phi points in the constant radius array  C
C     NTHERAD   : Number of the points in the constant radius array  C
C                                                                    C
C     IDEV      : Device number for output.                          C
C                                                                    C
C                  IDEV = 1 --> gif file (Landscape)'                C
C                         2 --> gif file. (Portrait)'                C
C                         3 --> ps file. (Landscape)'                C
C                         4 --> ps file. (Portrait)'                 C
C                         5 --> colour ps file. (Landscape)'         C
C                         6 --> colour ps file. (Portrait)'          C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     FRAD      : Array for constant radius plot.                    C
C                 Dim ( NPHIRAD, NTHERAD )                           C
C                                                                    C
C     RADV2     : Value of outer boundary radius.                    C
C                                                                    C
C     PAPWIDTH  : Paper width (inches)                               C
C                                                                    C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     EARCM     : Euler angle matrix Dim( 3, 3 ). (See EARCMC)       C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     STEM      : Filename stem (terminate with ' ')                 C
C                                                                    C
C*********************************************************************
      SUBROUTINE SIMPLE_SPHERE_PLOT( NLEV, NPHIRAD, NTHERAD, FRAD,
     1    RADV2, NLEVM, CONTL, PAPWIDTH, EARCM, HUEPOS, HUENEG, CS,
     2    SCAL, IDEV, STEM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NLEV, NPHIRAD, NTHERAD, NLEVM, IDEV
      REAL    FRAD( NPHIRAD, NTHERAD ), RADV2
      REAL    CONTL( NLEVM ), PAPWIDTH, HUEPOS, HUENEG, CS, SCAL
      DOUBLE PRECISION EARCM( 3, 3 )
      CHARACTER *(*) STEM
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NLEV2, ICONTOUR, IDEP, IPFLAG, ISTAT, ILEV
      REAL    C1V1, C1V2, C2V1, C2V2, COORD, VALMIN, VALMAX,
     1        DELTA, RATIO, PI,
     2        ZERO
      REAL    RRED( 17 ), RRED1,
     1        RBLUE( 17 ), RBLUE1,
     2        RGREEN( 17 ), RGREEN1
      PARAMETER ( ICONTOUR = 2, DELTA = 0.01, PI = 3.1415926,
     1            ZERO = 0.0 )
      CHARACTER *(1) CHDEV
      CHARACTER *(100) LINE
      REAL XLEFT, XRIGHT, YBOT, YTOP, X1, X2, Y1, Y2
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NLEV.EQ.-1 ) THEN
        IF ( NLEVM.GE.17 ) THEN
          NLEV2 = 17
        ELSE
          PRINT *,' Subroutine SIMPLE_SPHERE_PLOT.'
          PRINT *,' NLEVM = ', NLEVM
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDIF
C
      IF ( NLEV.GE.3 .AND. NLEV.LE.NLEVM ) THEN
        NLEV2 = NLEV
      ELSE
        IF ( NLEV.NE.-1 ) THEN
          PRINT *,' Subroutine SIMPLE_SPHERE_PLOT.'
          PRINT *,' NLEV = ',NLEV
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDIF
C
C Calculate maximum and minimum values
C
      VALMIN = 1.0e8
      VALMAX = -1.0e8
C
      CALL VALMIN_VALMAX_UPDATE( NPHIRAD, NTHERAD, FRAD,
     1                           VALMIN, VALMAX )
C
C Open the plotting device
C
      CALL PGPLOT_FILE_INIT( ISTAT, IDEV, STEM, CHDEV )
      RATIO = 1.0
      CALL PGPAP( PAPWIDTH, RATIO )
      CALL PGPAGE
C
C Now set the colours:
C RGB colour scheme only:
C
      IF ( NLEV.EQ.-1 ) THEN
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
C
      ENDIF
C
      XLEFT   = 0.10
      XRIGHT  = 0.90
      YBOT    = 0.10
      YTOP    = 0.90
      CALL PGSVP( XLEFT, XRIGHT, YBOT, YTOP )
      X1    = REAL( RADV2 )*(-1.0)
      X2    = REAL( RADV2 )
      Y1    = REAL( RADV2 )*(-1.0)
      Y2    = REAL( RADV2 )
      CALL PGSWIN( X1, X2, Y1, Y2 )
C
C Now begin the plotting sections
C First do constant radius plot
C
      IDEP   = 1
      IPFLAG = 3
      C1V1   = 0.0
      C1V2   = 2.0*PI
      C2V1   = 0.0
      C2V2   = PI
      COORD  = RADV2
C
      IF ( NLEV.EQ.-1 ) THEN
C
      CALL ARRAY_SPHER_PROJ__COLOUR2( FRAD, C1V1, C1V2, C2V1, C2V2,
     1            NPHIRAD, NTHERAD, NLEV2, CONTL, DELTA, ICONTOUR,
     2            VALMIN, VALMAX, IPFLAG, EARCM, COORD, IDEP )
C
      ELSE
C
      CALL ARRAY_SPHER_PROJ__COLOUR( FRAD, C1V1, C1V2, C2V1, C2V2,
     1            NPHIRAD, NTHERAD, NLEV2, CONTL, HUEPOS, HUENEG, CS,
     2            SCAL, DELTA, ICONTOUR, VALMIN, VALMAX,
     3            IPFLAG, EARCM, COORD, IDEP )
C
      ENDIF
C
C Now draw on spherical lines
C
      IF ( CHDEV.EQ.'P' ) THEN
        CALL PGSCI( 1 )
      ELSE
        CALL PGSCI( 0 )
      ENDIF
C
C Lines of longitude
C
c uncomment next two lines for phi = 0.0 line
      WRITE ( LINE, 171 ) RADV2, ZERO, ZERO, PI, 2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
c uncomment next two lines for phi = 90 degree line
      WRITE ( LINE, 171 ) RADV2, 0.5*PI, ZERO, PI, 2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 171 ) RADV2, PI, ZERO, PI, 2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 171 ) RADV2, 1.5*PI, ZERO, PI, 2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
C
C Lines of latitude
C
      WRITE ( LINE, 172 ) RADV2, 0.5*PI, 0.0, 2.0*PI,
     1                    3, 100, 1
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 172 ) RADV2, 0.7853981, 0.0, 2.0*PI, 
     1                   2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 172 ) RADV2, 2.3561945, 0.0, 2.0*PI, 
     1                   2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
C
C Spherical line
C
      WRITE ( LINE, 174 ) RADV2, 2, 200, 1
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
C
      CALL PGCLOS
C
 171  FORMAT ('CRP ',f12.7,f12.7,f12.7,f12.7,i4,i5,i4)
 172  FORMAT ('CRT ',f12.7,f12.7,f12.7,f12.7,i4,i5,i4)
 173  FORMAT ('CPT ',f12.7,f12.7,f12.7,f12.7,i4,i5,i4)
 174  FORMAT ('CRD ',f12.7,I4,I4,I4)
C
      RETURN
      END
C*********************************************************************

