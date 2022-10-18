C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C Fri Nov 16 11:32:01 WET 2001                                       C
C                                                                    C
C ARROWS_CONST_R3 - try and contour a function in Cartesian co-ords  C
C on a constant radius, (theta-phi) grid.                            C
C                                                                    C
C*********************************************************************
      PROGRAM ARROWS_CONST_R3
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, NTHMAX, NLEVM, LHMAX, NHMAX, ISVMAX, NNDM,
     1        NPMAX, NPHMAX
      PARAMETER ( NRMAX = 250, NTHMAX = 250, NPHMAX = 250, NLEVM = 20,
     1            LHMAX = 124, NHMAX = 6000, ISVMAX = NRMAX*NHMAX,
     2            NNDM = 6, NPMAX = (LHMAX+1)*(LHMAX+2)/2 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      REAL    F( NPHMAX, NTHMAX ),
     1        VC1( NPHMAX, NTHMAX ),
     2        VC2( NPHMAX, NTHMAX )
      REAL    CONTL( NLEVM )
C
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     2        IWORK( NNDM )
C
      DOUBLE PRECISION XARR( NRMAX ), VEC( ISVMAX ),
     1                 WORK1( NNDM ), WORK2( NNDM ),
     2                 COEFM( NNDM, NNDM )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C     .
C     . Variables which define our window:
C     .
      INTEGER          NTHE, NPHI
      REAL             THETA1, THETA2, PHI1, PHI2, LAT1, LAT2,
     1                 LONG1, LONG2
      COMMON  / PARAMC / THETA1, THETA2, PHI1, PHI2
C     .
C     . General
C     .
      INTEGER          NLEV, IWCONT, LUIN, LU, IOP, IAS, NNDS,
     1                 ISTAT, IDEV, INARR( 3 ), NH, ICOMP, ICONT,
     2                 LH, NR, NR1, I, ILEN, NPBA, ICONTOUR, RMFLAG,
     3                 INS, INW, IPS, IPW, INDNEG, INDPOS, IWARROW
      CHARACTER *(1)   CHDEV, ESCPCH
      CHARACTER *(180) STEM, LINE, FNAME
      REAL             PI, DELTA, ZERO, MINUS1, RLONG, RHEAD,
     1                 VALMIN, VALMAX, RWIDTH, RRATIO,
     2                 RNHUE, RNLIGHT, RNSAT, RPHUE, RPLIGHT, RPSAT
      DOUBLE PRECISION FRAD, RAD
      LOGICAL          OCOLOURS, OARROWS
C     .
C     . Variables for PGSVP
C     .
      REAL    XLEFT, XRIGHT, YBOT, YTOP
C
C Variables for setting the colour.
C
      REAL    HUEPOS, HUENEG, CS, SCAL
C     .
C     . Variables for PGSWIN
C     .
      REAL    X1, X2, Y1, Y2
C     .
      PARAMETER ( PI=3.1415926, ZERO = 0.0, MINUS1 = -1.0 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO I = 1, 180
        STEM(I:I) = ' '
      ENDDO
C
C First read output file-name stem
C
      ESCPCH = '*'
      LU     = 15
      LUIN   = 5
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
      CALL SVFRD( INARR, LU, NRMAX, VEC, FNAME )
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
      CALL XARRRD( NR, NRMAX, XARR, LU, FNAME )
C
      IF ( NR1.NE.NR ) THEN
        PRINT *,' Solution vector and radial node '
        PRINT *,' file claim differing numbers of '
        PRINT *,' grid nodes. Program aborted.'
        STOP
      ENDIF
C
      PRINT *,' Enter LONG1, LONG2, LAT1, LAT2, RAD.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) LONG1, LONG2, LAT1, LAT2, FRAD
C
      IF ( LONG2.LE.LONG1 .OR. LAT2.LE.LAT1 .OR.
     1     LAT1.LE.-90.0  .OR. LAT2.GE.90.0 .OR.
     2     FRAD.LT.0.0d0  .OR. FRAD.GT.1.0d0     ) THEN
        PRINT *,' LONG1 = ', LONG1
        PRINT *,' LONG2 = ', LONG2
        PRINT *,' LAT1  = ', LAT1
        PRINT *,' LAT2  = ', LAT2
        PRINT *,' FRAD   = ', FRAD
        PRINT *,' We need longitude and latitude in degrees.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RAD = XARR( 1 ) + ( XARR(NR) - XARR(1) )*FRAD
C
C Now convert long.s and lat.s into theta and phi
C
      PHI1 = LONG1*PI/180.0
      PHI2 = LONG2*PI/180.0
      THETA1 = 0.5*PI - LAT1*PI/180.0
      THETA2 = 0.5*PI - LAT2*PI/180.0
C
      PRINT *,' Enter NTHE, NPHI, NNDS.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NTHE, NPHI, NNDS
C
      IF ( NTHE.GT.NTHMAX .OR. NPHI.GT.NPHMAX .OR.
     1       NNDS.GT.NNDM     ) THEN
        PRINT *,' NNDS   = ', NNDS
        PRINT *,' NNDM   = ', NNDM
        PRINT *,' NTHE   = ', NTHE
        PRINT *,' NTHMAX = ', NTHMAX
        PRINT *,' NPHI   = ', NPHI
        PRINT *,' NPHMAX = ', NPHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      PRINT *,' Enter XLEFT, XRIGHT, YBOT, YTOP, RWIDTH, RRATIO'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) XLEFT, XRIGHT, YBOT, YTOP, RWIDTH, RRATIO
C
      IF (  XLEFT.LT.0.0 .OR. XLEFT.GT.1.0 .OR. XLEFT.GE.XRIGHT 
     1      .OR. XRIGHT.LT.0.0 .OR. XRIGHT.GT.1.0 .OR. YBOT.LT.0.0
     2      .OR. YBOT.GT.1.0 .OR. YBOT.GE.YTOP .OR. YTOP.LT.0.0
     3      .OR. YTOP.GT.1.0 )          THEN
        PRINT *,' XLEFT  = ', XLEFT 
        PRINT *,' XRIGHT = ', XRIGHT
        PRINT *,' YBOT   = ', YBOT
        PRINT *,' YTOP   = ', YTOP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Enter number of contour levels
C
      PRINT *,' Enter NLEV, IDEV, ICONT, ICOMP, IAS '
      PRINT *,' IDEV: 1 --> gif file. (Landscape)'
      PRINT *,'       2 --> gif file. (Portrait)'
      PRINT *,'       3 --> ps file. (Landscape)'
      PRINT *,'       4 --> ps file. (Portrait)'
      PRINT *,'       5 --> colour ps file. (Landscape)'
      PRINT *,'       6 --> colour ps file. (Portrait)'
      PRINT *,' ICONT: 1 --> colours only.'
      PRINT *,'        2 --> contours only.'
      PRINT *,'        3 --> colours and contours.'
      PRINT *,' ICOMP: 1 -> v_r - radial component of velocity '
      PRINT *,' 2 -> v_{theta} theta component of velocity '
      PRINT *,' 3 -> v_{phi} phicomponent of velocity '
      PRINT *,' 4 -> B_r - radial component of magnetic field. '
      PRINT *,' 5 -> B_{theta} theta component of magnetic field. '
      PRINT *,' 6 -> B_{phi} phicomponent of magnetic field. '
      PRINT *,' 7 -> Theta '
      PRINT *,' IAS: 0 -> all components.'
      PRINT *,' IAS: 1 -> only axisymmetric components.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NLEV, IDEV, ICONT, ICOMP, IAS
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
C Make sure these are valid
C
      IF ( ICONT.EQ.1 ) THEN
        OARROWS  = .FALSE.
        OCOLOURS = .TRUE.
      ENDIF
C
      IF ( ICONT.EQ.2 ) THEN
        OARROWS  = .TRUE.
        OCOLOURS = .FALSE.
      ENDIF
C
      IF ( ICONT.EQ.3 ) THEN
        OARROWS  = .TRUE.
        OCOLOURS = .TRUE.
      ENDIF
C
      IF ( ICONT.EQ.4 ) THEN
        OARROWS  = .TRUE.
        OCOLOURS = .FALSE.
      ENDIF
C
      IF ( ICONT.LT.1 .OR. ICONT.GT.4 ) THEN
        PRINT *,' ICONT = ', ICONT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NLEV.GT.NLEVM ) THEN
        PRINT *,' NLEV  = ', NLEV 
        PRINT *,' NLEVM = ', NLEVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Enter huepos, hueneg, csat, scal
C
      PRINT *,' Enter huepos, hueneg, csat, scal, iw'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) HUEPOS, HUENEG, CS, SCAL, IWCONT
C
C Enter parameters for positive contours
C (only effective if ICONT = 4)
C
      PRINT *,' Enter rphue, rpsat, rplight, ipw, ips'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) RPHUE, RPSAT, RPLIGHT, IPW, IPS
C
C Enter parameters for negative contours
C (only effective if ICONT = 4)
C
      PRINT *,' Enter rnhue, rnsat, rnlight, inw, ins'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) RNHUE, RNSAT, RNLIGHT, INW, INS
C
      PRINT *,' INW = ',INW,' INS = ',INS
      PRINT *,' IPW = ',IPW,' IPS = ',IPS
C
C Get information for arrow drawing.
C
      PRINT *,' Enter NPBA, RLONG, RHEAD, IWARROW.'
      PRINT *,' npba is the number of grid points between arrows.'
      PRINT *,' rlong is the length of the longest arrow.'
      PRINT *,' rhead is the size of the biggest arrowhead.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NPBA, RLONG, RHEAD, IWARROW
C
C Enter icontour, valmin, valmax
C
      PRINT *,' Enter icontour, valmin, valmax.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) ICONTOUR, VALMIN, VALMAX
C
C This inputs a device from the user.
C
      CALL PGPLOT_FILE_INIT( ISTAT, IDEV, STEM, CHDEV )
      PRINT *,' ISTAT = ', ISTAT
      PRINT *,' CHDEV = ', CHDEV
C
C Fill data array for colours
C
      IOP = 0
      CALL RMATOP( F, ZERO, NPHI, NTHE, IOP )
      CALL RMATOP( VC1, ZERO, NPHI, NTHE, IOP )
      CALL RMATOP( VC2, ZERO, NPHI, NTHE, IOP )
C
      CALL CONSTANT_R_RECT_EVAL( ICOMP, INARR, NNDS,
     1        IWORK, MHT, MHL, MHM, LH, IAS, RAD, VEC, XARR,
     2        WORK1, WORK2, COEFM, F, NPHI, NTHE, THETA1,
     3        THETA2, PHI1, PHI2 )
      CALL ARR_RM_MEAN( NPHI, NTHE, F, RMFLAG )
C
C Fill data array for arrows:
C iop = 1 --> puts phi component into VC1
C
      IOP = 3
      CALL CONSTANT_R_RECT_EVAL( IOP, INARR, NNDS,
     1           IWORK, MHT, MHL, MHM, LH, IAS, RAD, VEC, XARR,
     2           WORK1, WORK2, COEFM, VC1, NPHI, NTHE, THETA1,
     3           THETA2, PHI1, PHI2 )
C
C iop = 2 puts theta component into VC2
C
      IOP = 2
      CALL CONSTANT_R_RECT_EVAL( IOP, INARR, NNDS,
     1           IWORK, MHT, MHL, MHM, LH, IAS, RAD, VEC, XARR,
     2           WORK1, WORK2, COEFM, VC2, NPHI, NTHE, THETA1,
     3           THETA2, PHI1, PHI2 )
C
C However, v_theta is of the wrong sign for plotting
C so multiply all values by -1.0
C
      IOP = 1
      CALL RMATOP( VC2, MINUS1, NPHI, NTHE, IOP )
C
C Call PGPAGE to clear or initialise a new page.
C
      CALL PGPAP( RWIDTH, RRATIO )
      CALL PGPAGE
C
      CALL PGSVP( XLEFT, XRIGHT, YBOT, YTOP )
C
      X1 =    PHI1
      Y1 =    THETA1
      X2 =    PHI2
      Y2 =    THETA2
      PRINT *,' (X1,Y1) = (', X1,',',Y1,')'
      PRINT *,' (X2,Y2) = (', X2,',',Y2,')'
C
      CALL PGSWIN( X1, X2, Y1, Y2 )
C 
C Begin batch to buffer
C 
      CALL PGBBUF
C 
      DELTA = 0.01
      IF ( OCOLOURS )
     1      CALL ARRAY_RECT_COLOUR_OPTION( F, NPHI, NTHE, NLEV,
     2                         CONTL, HUEPOS, HUENEG, CS, SCAL,
     3                         DELTA, ICONTOUR, VALMIN, VALMAX )
C
      INDNEG = 43
      INDPOS = 44
      CALL PGSHLS( INDNEG, RNHUE, RNLIGHT, RNSAT )
      CALL PGSHLS( INDPOS, RPHUE, RPLIGHT, RPSAT )
C
      IF ( ICONT.EQ.4 ) 
     1      CALL ARRAY_RECT_DRAW_OPTION_MOD( F, NPHI, NTHE, NLEV,
     2                    CONTL, DELTA, ICONTOUR, VALMIN, VALMAX,
     3                    INW, IPW, INS, IPS, INDNEG, INDPOS )

      IF ( CHDEV.EQ.'P' ) THEN
        CALL PGSCI( 1 )
      ELSE
        IF ( OCOLOURS ) THEN
          CALL PGSCI( 0 )
        ELSE
          CALL PGSCI( 1 )
        ENDIF
      ENDIF
C
      CALL PGSLW( 1 )
C
      IF ( OARROWS )
     1      CALL ARRAY_RECT_ARROWS( VC1, VC2, NPBA, RLONG, RHEAD,
     2                              IWARROW, NPHI, NTHE )
C
      CALL PGSLW( IWCONT )
      CALL CONTINENT_DRAW( 71, 'coast.dat' )
C
      IF ( ICONT.EQ.4 ) THEN
        CALL PGSLW( 2 )
        CALL PGMOVE( X1      , Y1-DELTA )
        CALL PGDRAW( X1      , Y2+DELTA )
        CALL PGDRAW( X2      , Y2+DELTA )
        CALL PGDRAW( X2      , Y1-DELTA )
        CALL PGDRAW( X1      , Y1-DELTA )
      ENDIF
C
C End batch to buffer
C 
      CALL PGEBUF
C 
C Close the graphics devices.
C
      CALL PGCLOS
      STOP
      END
C*********************************************************************

C*********************************************************************
      SUBROUTINE CONTINENT_DRAW( LUCONT, FNCONT )
      IMPLICIT NONE
C_____________________________________________________________________
      INTEGER          LUCONT
      CHARACTER *(*)   FNCONT
C____________________________________________________________________C
      REAL             THETA1, THETA2, PHI1, PHI2
      COMMON  / PARAMC / THETA1, THETA2, PHI1, PHI2
C____________________________________________________________________C
      INTEGER          IREAD, ICONT, NDAT, IDAT
      REAL             RLAT, RLONG, RTHE, RPHI, PI, HALFPI, PIBYOE,
     1                 XWORLD, YWORLD, TWOPI, CLONG, SLONG, BLONG,
     2                 OLDXWORLD, OLDRPHI, XW2
      PARAMETER   ( IREAD = 1, PI = 3.141592 )
      CHARACTER *(200) LINE
      LOGICAL   DRAW
C____________________________________________________________________C
      HALFPI = 0.5*PI
      TWOPI  = 2.0*PI
      PIBYOE = PI/180.0
C
      CALL FOPEN( LUCONT, FNCONT, IREAD )
C
 50   CONTINUE
      READ( LUCONT, 80, END=60 ) LINE
      READ ( LINE, * ) ICONT, NDAT
c     print *,' continent ',ICONT,' ndata = ', NDAT
      DO IDAT = 1, NDAT
        READ ( LUCONT, * ) RLAT, RLONG
C
        IF ( IDAT.GT.1 ) THEN
          OLDXWORLD = XW2
          OLDRPHI   = RPHI
        ENDIF
C
C Convert longitude and latitude to theta and phi
C
        RTHE  = HALFPI - RLAT*PIBYOE
        RPHI  = RLONG*PIBYOE
        SLONG = SIN( RPHI )
        CLONG = COS( RPHI )
C
        BLONG = ATAN2( SLONG, CLONG )
        IF ( BLONG.LT.0.0 ) THEN
          XW2 = BLONG + TWOPI
        ELSE
          XW2 = BLONG
        ENDIF
        IF ( BLONG.LT.PHI1 ) THEN
          XWORLD = BLONG + TWOPI
        ELSE
          XWORLD = BLONG
        ENDIF
C
C xw2 is the longitude in radians between 0 and 2*pi
C
        YWORLD = RTHE
C
c       print *,' rphi = ', rphi,' xw = ',xworld,
c    1              ' p1 = ', phi1,' p2 = ', phi2
C
        IF ( IDAT.EQ.1 ) THEN
          DRAW = .FALSE.
        ELSE
          DRAW = .TRUE.
          IF ( XWORLD.LT.PHI1 .OR. XWORLD.GT.PHI2 .OR.
     1         YWORLD.GT.THETA1 .OR. YWORLD.LT.THETA2 ) DRAW = .FALSE.
          IF ( (XW2-OLDXWORLD)*(RPHI-OLDRPHI).LT.0.0 ) DRAW = .FALSE.
        ENDIF
C
        IF ( DRAW ) THEN
          CALL PGDRAW( XWORLD, YWORLD )
        ELSE
          CALL PGMOVE( XWORLD, YWORLD )
        ENDIF
C
      ENDDO
      GOTO 50
 60   CONTINUE
      PRINT *,' All information read '
      CALL FCLOSE( LUCONT, FNCONT, 'Error' )
 80   FORMAT(A)
C
      RETURN
      END
C
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
