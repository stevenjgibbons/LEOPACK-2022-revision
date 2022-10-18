
C*********************************************************************
C subroutine CUT OUT SPHERE PLOT 2 ***********************************
C            --- --- ------ ---- - ***********************************
C Steve Gibbons Thu Oct 11 09:07:19 WEST 2001                        C
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
C     NRADP1    : Number of radial nodes in FP1 array.               C
C     NTHEP1    : Number of theta nodes in FP1 array.                C
C                                                                    C
C     NRADP2    : Number of radial nodes in FP2 array.               C
C     NTHEP2    : Number of theta nodes in FP2 array.                C
C                                                                    C
C     NRADTHE   : Number of radial nodes in FTHE array.              C
C     NPHITHE   : Number of phi nodes in FTHE array.                 C
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
C     IOPT      : IOPT = 0 --> Simply modify VALMIN and VALMAX vals  C
C                              No plotting is done.                  C
C                 IOPT = 1 --> Plot using previously calculated      C
C                              values for VALMIN and VALMAX.         C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     FRAD      : Array for constant radius plot.                    C
C                 Dim ( NPHIRAD, NTHERAD )                           C
C                                                                    C
C     RADV1     : Value of inner boundary radius.                    C
C     RADV2     : Value of outer boundary radius.                    C
C                                                                    C
C     FP1       : Array for first constant phi plot.                 C
C                 Dim ( NRADP1, NTHEP1 )                             C
C                                                                    C
C     PHIV1     : First constant phi value.                          C
C                                                                    C
C     DTHETA    : Constant theta value.                              C
C                                                                    C
C     FP2       : Array for second constant phi plot.                C
C                 Dim ( NRADP1, NTHEP1 )                             C
C                                                                    C
C     PHIV2     : Second constant phi value.                         C
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
C     VALMIN    : Minimum value to be contoured.                     C
C     VALMAX    : Maximum value to be contoured.                     C
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
      SUBROUTINE CUT_OUT_SPHERE_PLOT2( NLEV, NPHIRAD, NTHERAD, FRAD,
     1    RADV1, RADV2, NRADP1, NTHEP1, FP1, PHIV1, NRADP2, NTHEP2,
     2    FP2, PHIV2, NRADTHE, NPHITHE, FTHE, DTHETA, NLEVM, CONTL,
     3    PAPWIDTH, EARCM, HUEPOS, HUENEG, CS, SCAL, IDEV, STEM,
     4    IOPT, VALMIN, VALMAX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NLEV, NPHIRAD, NTHERAD, NRADP1, NTHEP1, NRADP2, NTHEP2,
     1        NRADTHE, NPHITHE, NLEVM, IDEV, IOPT
      REAL    FRAD( NPHIRAD, NTHERAD ), RADV1, RADV2,
     1        FP1( NRADP1, NTHEP1 ), PHIV1, VALMIN,
     2        FP2( NRADP2, NTHEP2 ), PHIV2, VALMAX,
     3        FTHE( NRADTHE, NPHITHE ), DTHETA
      REAL    CONTL( NLEVM ), PAPWIDTH, HUEPOS, HUENEG, CS, SCAL
      DOUBLE PRECISION EARCM( 3, 3 )
      CHARACTER *(*) STEM
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NLEV2, ICONTOUR, IDEP, IPFLAG, ISTAT, ILEV
      REAL    C1V1, C1V2, C2V1, C2V2, COORD, ZERO, RTEMP,
     1        DELTA, C1B1, C1B2, C2B1, C2B2, RATIO, DTHE, DPHI, PI
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
          PRINT *,' Subroutine CUT_OUT_SPHERE_PLOT2.'
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
          PRINT *,' Subroutine CUT_OUT_SPHERE_PLOT2.'
          PRINT *,' NLEV = ',NLEV
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDIF
C
C Calculate maximum and minimum values
C
      IF ( IOPT.EQ.0 ) THEN
C
        CALL VALMIN_VALMAX_UPDATE( NPHIRAD, NTHERAD, FRAD,
     1                           VALMIN, VALMAX )
        CALL VALMIN_VALMAX_UPDATE( NRADP1, NTHEP1, FP1,
     1                           VALMIN, VALMAX )
        CALL VALMIN_VALMAX_UPDATE( NRADP2, NTHEP2, FP2,
     1                           VALMIN, VALMAX )
        CALL VALMIN_VALMAX_UPDATE( NRADTHE, NPHITHE, FTHE,
     1                           VALMIN, VALMAX )
        RETURN
      ENDIF
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
      DTHE = PI/FLOAT( NTHERAD - 1 )
      DPHI = 2.0*PI/FLOAT( NPHIRAD - 1 )
C
C We annihilate the part of the array corresponding to the
C cut-out:
C
      C1B1   = PHIV1 + DPHI*2.0
      C1B2   = PHIV2 - DPHI*2.0
      C2B1   = DTHE*2.0
      C2B2   = DTHETA - DTHE*2.0
C
      IDEP   = 1
      IPFLAG = 3
      C1V1   = 0.0
      C1V2   = 2.0*PI
      C2V1   = 0.0
      C2V2   = PI
      COORD  = RADV2
C
      CALL ARRAY_ANNIHILATE( NPHIRAD, NTHERAD, FRAD,
     1                       C1V1, C1V2, C2V1, C2V2,
     2                       C1B1, C1B2, C2B1, C2B2 )
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
C Now to constant phi 1 section
C
      IDEP   = 0
      IPFLAG = 1
      C1V1   = RADV1
      C1V2   = RADV2
      C2V1   = 0.0
      C2V2   = DTHETA
      COORD  = PHIV1
C
      IF ( NLEV.EQ.-1 ) THEN
C
      CALL ARRAY_SPHER_PROJ__COLOUR2( FP1, C1V1, C1V2, C2V1, C2V2,
     1            NRADP1, NTHEP1, NLEV2, CONTL, DELTA, ICONTOUR,
     2            VALMIN, VALMAX, IPFLAG, EARCM, COORD, IDEP )
C
      ELSE
C
      CALL ARRAY_SPHER_PROJ__COLOUR( FP1, C1V1, C1V2, C2V1, C2V2,
     1            NRADP1, NTHEP1, NLEV2, CONTL, HUEPOS, HUENEG, CS,
     2            SCAL, DELTA, ICONTOUR, VALMIN, VALMAX,
     3            IPFLAG, EARCM, COORD, IDEP )
C
      ENDIF
C
C Now to constant phi 2 section
C
      IDEP   = 0
      IPFLAG = 1
      C1V1   = RADV1
      C1V2   = RADV2
      C2V1   = 0.0
      C2V2   = DTHETA
      COORD  = PHIV2
C
      IF ( NLEV.EQ.-1 ) THEN
C
      CALL ARRAY_SPHER_PROJ__COLOUR2( FP2, C1V1, C1V2, C2V1, C2V2,
     1            NRADP2, NTHEP2, NLEV2, CONTL, DELTA, ICONTOUR,
     2            VALMIN, VALMAX, IPFLAG, EARCM, COORD, IDEP )
C
      ELSE
C
      CALL ARRAY_SPHER_PROJ__COLOUR( FP2, C1V1, C1V2, C2V1, C2V2,
     1            NRADP2, NTHEP2, NLEV2, CONTL, HUEPOS, HUENEG, CS,
     2            SCAL, DELTA, ICONTOUR, VALMIN, VALMAX,
     3            IPFLAG, EARCM, COORD, IDEP )
C
      ENDIF
C
C Now to constant theta section
C
      IDEP   = 0
      IPFLAG = 2
      C1V1   = RADV1
      C1V2   = RADV2
      C2V1   = PHIV1
      C2V2   = PHIV2
      COORD  = DTHETA
C
      IF ( (PI-DTHETA).GT.0.001 ) THEN
        IF ( NLEV.EQ.-1 ) THEN
C
        CALL ARRAY_SPHER_PROJ__COLOUR2( FTHE, C1V1, C1V2, C2V1, C2V2,
     1            NRADTHE, NPHITHE, NLEV2, CONTL, DELTA, ICONTOUR,
     2            VALMIN, VALMAX, IPFLAG, EARCM, COORD, IDEP )
C
        ELSE
C
        CALL ARRAY_SPHER_PROJ__COLOUR( FTHE, C1V1, C1V2, C2V1, C2V2,
     1            NRADTHE, NPHITHE, NLEV2, CONTL, HUEPOS, HUENEG, CS,
     2            SCAL, DELTA, ICONTOUR, VALMIN, VALMAX,
     3            IPFLAG, EARCM, COORD, IDEP )
C
        ENDIF
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
      WRITE ( LINE, 171 ) RADV2, ZERO, ZERO, PI, 2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 171 ) RADV2, 0.5*PI, ZERO, PI, 2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 171 ) RADV2, PI, ZERO, PI, 2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 171 ) RADV2, 1.5*PI, ZERO, PI, 2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
C
      RTEMP = MIN( PI, DTHETA )
      WRITE ( LINE, 171 ) RADV1, PHIV1, ZERO, RTEMP, 4, 100, 0
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 171 ) RADV1, PHIV2, ZERO, RTEMP, 4, 100, 0
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
C
      WRITE ( LINE, 171 ) RADV2, PHIV1, ZERO, RTEMP, 4, 100, 0
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 171 ) RADV2, PHIV1, DTHETA, PI, 2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 171 ) RADV2, PHIV2, ZERO, RTEMP, 4, 100, 0
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 171 ) RADV2, PHIV2, DTHETA, PI, 2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
C
C Lines of latitude
C
      WRITE ( LINE, 172 ) RADV2, 0.5*PI, PHIV2, PHIV1+2.0*PI,
     1                    3, 100, 1
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 172 ) RADV2, 0.7853981, PHIV2, PHIV1+2.0*PI, 
     1                   2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 172 ) RADV2, 2.3561945, PHIV2, PHIV1+2.0*PI, 
     1                   2, 100, 2
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      WRITE ( LINE, 172 ) RADV2, DTHETA, PHIV1, PHIV2, 4, 100, 1
      CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      IF ( (PI-DTHETA).GT.0.001 ) THEN
        WRITE ( LINE, 172 ) RADV1, DTHETA, PHIV1, PHIV2, 4, 100, 1
        CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      ENDIF
C
C Radial lines
C
      IF ( DTHETA.GT.0.4   ) THEN
        WRITE ( LINE, 173 ) ZERO, ZERO, RADV1, RADV2, 4, 40, 0
        CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
        IF ( (PI-DTHETA).GT.0.001 ) THEN
          WRITE ( LINE, 173 ) PHIV1, DTHETA, RADV1, RADV2, 4, 40, 0
          CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
          WRITE ( LINE, 173 ) PHIV2, DTHETA, RADV1, RADV2, 4, 40, 0
          CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
        ELSE
          WRITE ( LINE, 173 ) PHIV1, PI, RADV1, RADV2, 4, 40, 0
          CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
        ENDIF
      ENDIF
C
C Spherical line
C
      IF ( (PI-DTHETA).GT.0.001 ) THEN
        WRITE ( LINE, 174 ) RADV2, 2, 200, 1
        CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
      ENDIF
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

