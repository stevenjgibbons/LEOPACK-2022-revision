C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C Fri Jan 11 08:31:33 WET 2002                                       C
C                                                                    C
C Full Sphere Plot with Continents                                   C
C                                                                    C
C*********************************************************************
      PROGRAM continents_full_sphere_plot
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      REAL    THETA1, THETA2, PHI1, PHI2
      COMMON  / PARAMC / THETA1, THETA2, PHI1, PHI2
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRADMX, NTHMAX, NPHMAX, NLEVM, LHMAX, NHMAX,
     1        ISVMAX, NNDM
      PARAMETER ( NRADMX = 250, NTHMAX = 250, NPHMAX = 250,
     1            NLEVM = 20, LHMAX = 160, NHMAX = 3000,
     2            ISVMAX = NRADMX*NHMAX, NNDM = 6 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      REAL    FUNCRAD( NPHMAX, NTHMAX )
      REAL    CONTL( NLEVM )
C
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     2        IWORK( NNDM )
C
      DOUBLE PRECISION XARR( NRADMX ), VEC( ISVMAX ),
     1                 WORK1( NNDM ), WORK2( NNDM ),
     2                 COEFM( NNDM, NNDM )
C
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C     .
C     . General
C     .
      INTEGER          NLEV, LUIN, LU, IOP, IAS, NNDS,
     1                 IDEV, INARR( 3 ), NH, ICOMP, RMFLAG,
     2                 LH, NR, NR1, I, ILEN, NTHP, NPHI
      CHARACTER *(1)   ESCPCH
      CHARACTER *(180) STEM, LINE, FNAME, STEM2
      REAL             PI, ZERO, PAPWIDTH, RADV2
      DOUBLE PRECISION RAD
      DOUBLE PRECISION ALPHAD, BETAD, GAMMAD, EARCM( 3, 3 )
C     .
C Variables for setting the colour.
C
      REAL    HUEPOS, HUENEG, CS, SCAL
C     .
      PARAMETER ( PI=3.1415926, ZERO = 0.0, IAS = 0 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      ESCPCH = '*'
      LUIN   = 5
      PRINT *,' Enter ALPHAD, BETAD, GAMMAD.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) ALPHAD, BETAD, GAMMAD
      PRINT *,' ALPHAD = ', ALPHAD
      PRINT *,' BETAD  = ', BETAD
      PRINT *,' GAMMAD = ', GAMMAD
C
      CALL EARCMC( ALPHAD, BETAD, GAMMAD, EARCM )
      CALL EARMIR( EARCM )
C
      DO I = 1, 180
        STEM(I:I) = ' '
        STEM2(I:I) = ' '
      ENDDO
C
C First read output file-name stem
C
      LU     = 15
      PRINT *,' Enter filename stem.'
      CALL LINERD( LUIN, STEM, ESCPCH )
C
C Read integers file:
C
      PRINT *,' Enter integers file name.'
      CALL LINERD( LUIN, LINE, ESCPCH )
C
      DO I = 1, 180
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 50
        ENDIF
      ENDDO
 50   CONTINUE
      FNAME = LINE(1:ILEN)
      CALL BIHFRD( NH, NHMAX, MHT, MHL, MHM, LU, FNAME )
      INARR( 3 ) = NH
      LH = 0
      DO I = 1, NH
        IF ( MHL( I ).GT.LH ) LH = MHL( I )
      ENDDO
      IF ( LH.GT.LHMAX ) THEN
        PRINT *,' LH = ', LH
        PRINT *,' LHMAX = ', LHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now read in vectors file
C
      PRINT *,' Enter vector file name.'
      CALL LINERD( LUIN, LINE, ESCPCH )
C
      DO I = 1, 180
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 51
        ENDIF
      ENDDO
 51   CONTINUE
      FNAME = LINE(1:ILEN)
C
      CALL SVFRD( INARR, LU, NRADMX, VEC, FNAME )
      NR1 = INARR( 2 )
C
C Now read in radial spacings filename
C
      PRINT *,' Enter radial spacings file name.'
      CALL LINERD( LUIN, LINE, ESCPCH )
C
      DO I = 1, 180
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 52
        ENDIF
      ENDDO
 52   CONTINUE
      FNAME = LINE(1:ILEN)
C      
      CALL XARRRD( NR, NRADMX, XARR, LU, FNAME )
C
      IF ( NR1.NE.NR ) THEN
        PRINT *,' Solution vector and radial node '
        PRINT *,' file claim differing numbers of '
        PRINT *,' grid nodes. Program aborted.'
        STOP
      ENDIF
C
C Enter huepos, hueneg, csat, scal
C
      PRINT *,' Enter huepos, hueneg, csat, scal'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) HUEPOS, HUENEG, CS, SCAL
C
C Enter number of contour levels
C
      PRINT *,' Enter NLEV, IDEV, NNDS, PAPWIDTH, ICOMP'
      PRINT *,' IDEV: 1 --> gif file. (Landscape)'
      PRINT *,'       2 --> gif file. (Portrait)'
      PRINT *,'       3 --> ps file. (Landscape)'
      PRINT *,'       4 --> ps file. (Portrait)'
      PRINT *,'       5 --> colour ps file. (Landscape)'
      PRINT *,'       6 --> colour ps file. (Portrait)'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NLEV, IDEV, NNDS, PAPWIDTH, ICOMP
C
C Modification SJG Thu Jan 15 12:11:26 MET 2004
      IF ( ICOMP.LT.0 ) THEN
        ICOMP  = -ICOMP
        RMFLAG = 1
      ELSE
        RMFLAG = 0
      ENDIF
C END Modification SJG Thu Jan 15 12:11:26 MET 2004
C
      IF ( NLEV.GT.NLEVM ) THEN
        PRINT *,' NLEV  = ', NLEV 
        PRINT *,' NLEVM = ', NLEVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Enter resolutions and positions
C
      PRINT *,' Enter NTHP, NPHI, RADV2'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NTHP, NPHI, RADV2
C
      IF ( NTHP.GT.NTHMAX .OR. NPHI.GT.NPHMAX ) THEN
        PRINT *,' NTHMAX = ', NTHMAX
        PRINT *,' NPHMAX = ', NPHMAX
        PRINT *,' NTHP   = ', NTHP
        PRINT *,' NPHI   = ', NPHI
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Fill data arrays
C
      IOP = 0
      CALL RMATOP( FUNCRAD, ZERO, NPHMAX, NTHMAX, IOP )
C
C Enter the data arrays
C First constant radius data
C
      PHI1   = 0.0
      PHI2   = 2.0*PI
      THETA1 = 0.0
      THETA2 = PI
C
      RAD   = DBLE( RADV2 )
      CALL CONSTANT_R_RECT_EVAL( ICOMP, INARR, NNDS, IWORK,
     1          MHT, MHL, MHM, LH, IAS, RAD, VEC, XARR, WORK1,
     2          WORK2, COEFM, FUNCRAD, NPHI, NTHP, THETA1, THETA2,
     3          PHI1, PHI2 )
C
      DO I = 1, 180
        IF ( STEM(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 53
        ENDIF
      ENDDO
 53   CONTINUE
C
      STEM2(1:ILEN)        = STEM(1:ILEN)
      STEM2(ILEN+1:ILEN+1) = ' '
      STEM2                = STEM2(1:ILEN+6)
C
      CALL ARR_RM_MEAN( NPHI, NTHP, FUNCRAD, RMFLAG )
C
      CALL SIMPLE_SPHERE_PLOT_CONTINENTS( NLEV, NPHI, NTHP, FUNCRAD,
     1    RADV2, NLEVM, CONTL, PAPWIDTH, EARCM, HUEPOS, HUENEG, CS,
     2    SCAL, IDEV, STEM2, 21, 'coast.dat' )
C
      STOP
      END
C*********************************************************************

C*********************************************************************
C subroutine SIMPLE SPHERE PLOT CONTINENTS ***************************
C            ------ ------ ---- ---------- ***************************
C Steve Gibbons Fri Jan 11 08:36:56 WET 2002                         C
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
C     LUCONT    : Device number for continents file.                 C
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
C     FNCONT    : Filename containing coast data                     C
C                                                                    C
C*********************************************************************
      SUBROUTINE SIMPLE_SPHERE_PLOT_CONTINENTS( NLEV, NPHIRAD,
     1    NTHERAD, FRAD, RADV2, NLEVM, CONTL, PAPWIDTH, EARCM,
     2    HUEPOS, HUENEG, CS, SCAL, IDEV, STEM, LUCONT, FNCONT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NLEV, NPHIRAD, NTHERAD, NLEVM, IDEV, LUCONT
      REAL    FRAD( NPHIRAD, NTHERAD ), RADV2
      REAL    CONTL( NLEVM ), PAPWIDTH, HUEPOS, HUENEG, CS, SCAL
      DOUBLE PRECISION EARCM( 3, 3 )
      CHARACTER *(*) STEM, FNCONT
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NLEV2, ICONTOUR, IDEP, IPFLAG, ISTAT, ILEV, IREAD
      REAL    C1V1, C1V2, C2V1, C2V2, COORD, VALMIN, VALMAX,
     1        DELTA, RATIO, PI,
     2        ZERO
      REAL    RRED( 17 ), RRED1,
     1        RBLUE( 17 ), RBLUE1,
     2        RGREEN( 17 ), RGREEN1
      PARAMETER ( ICONTOUR = 2, DELTA = 0.01, PI = 3.1415926,
     1            ZERO = 0.0, IREAD = 1 )
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
      valmin = valmin - 0.000001
      valmax = valmax + 0.000001
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
C Now draw on continents
C
      CALL PGSLW( 3 )
      CALL CONTINENT_DRAW( LUCONT, FNCONT, EARCM, RADV2 )
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

C*********************************************************************
      SUBROUTINE  CONTINENT_DRAW( LUCONT, FNCONT, EARCM, RRAD )
      IMPLICIT NONE
C_____________________________________________________________________
      INTEGER          LUCONT
      CHARACTER *(*)   FNCONT
      REAL             RRAD
      DOUBLE PRECISION EARCM( 3, 3 )
C____________________________________________________________________C
      INTEGER          IREAD, ICONT, NDAT, IDAT
      REAL             RLAT, RLONG, RTHE, RPHI, PI, HALFPI, PIBYOE,
     1                 XWORLD, YWORLD, DEPTH, SRAD, RADEAR
      PARAMETER   ( IREAD = 1, PI = 3.141592 )
      CHARACTER *(200) LINE
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      RADEAR = 6371.0
      HALFPI = 0.5*PI
      PIBYOE = PI/180.0
C
      CALL FOPEN( LUCONT, FNCONT, IREAD )
C
 50   CONTINUE
      READ( LUCONT, 80, END=60 ) LINE
      READ ( LINE, * ) ICONT, NDAT
      print *,' continent ',ICONT,' ndata = ', NDAT
      DO IDAT = 1, NDAT
        READ ( LUCONT, * ) RLAT, RLONG
C
C Convert longitude and latitude to theta and phi
C
        RTHE = HALFPI - RLAT*PIBYOE
        RPHI = RLONG*PIBYOE
C
C Convert into world coordinates
C
        CALL SPHER_SAT_2_WORLD( RRAD, RTHE, RPHI, EARCM, XWORLD,
     1                          YWORLD, DEPTH )
C
        SRAD = SQRT( XWORLD*XWORLD + YWORLD*YWORLD )
C
        IF ( IDAT.EQ.1 ) THEN
          CALL PGMOVE( XWORLD, YWORLD )
        ELSE
C         .
C         . OK so we are not on the first point
C         .
          IF ( DEPTH.LT.0.0 ) THEN
            CALL PGMOVE( XWORLD, YWORLD )
          ELSE
            CALL PGSCI( 1 )
            CALL PGDRAW( XWORLD, YWORLD )
          ENDIF
C         .
        ENDIF
C
      ENDDO
      GOTO 50
 60   CONTINUE
      PRINT *,' All information read '
      CALL FCLOSE( LUCONT, FNCONT, 'Error' )
 80   FORMAT(A)
C
C Draw circle radius
C
c     LINE = 'CRD  6371.0   1  100 1'
c     WRITE ( LINE(4:12), 84 ) RADEAR
c     CALL SPHERICAL_LINE_DRAW( EARCM, LINE )
c84   FORMAT(f9.2)
C
      RETURN
      END
C*********************************************************************

C*********************************************************************
C Steve Gibbons 15 January 2004
C If RMFLAG is not equal to 1, return with no action.
C If RMFLAG is equal to 1, then calculate the mean of the real
C array and subtract it from each element.
C We also return with no action if N1 or N2 is less than 1.
C--------------------------------------------------
      SUBROUTINE ARR_RM_MEAN( N1, N2, ARR, RMFLAG )
      IMPLICIT NONE
      INTEGER        RMFLAG, N1, N2
      REAL           ARR( N1, N2 )
C
      INTEGER        I1, I2
      REAL           RTOTEL, RMEANV
C
      IF ( RMFLAG.NE.1 ) RETURN
      IF ( N1.LT.1 )     RETURN
      IF ( N2.LT.1 )     RETURN
C
      RTOTEL = 1.0/FLOAT( N1*N2 )
      RMEANV = 0.0
      DO I2 = 1, N2
        DO I1 = 1, N1
          RMEANV = RMEANV + ARR( I1, I2 )*RTOTEL
        ENDDO
      ENDDO
C
C Have now calculated the mean - now subtract it ...
C
      DO I2 = 1, N2
        DO I1 = 1, N1
          ARR( I1, I2 ) = ARR( I1, I2 ) - RMEANV
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
