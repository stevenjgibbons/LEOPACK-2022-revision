C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C Mon Jun 25 10:01:49 WEST 2001                                      C
C                                                                    C
C EQ_MERID1 - contour either a meridian or equatorial section        C
C             of a function.                                         C
C                                                                    C
C*********************************************************************
      PROGRAM ps_2plot_z_eq_merid2
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, NTMAX, NLEVM, LHMAX, NHMAX, ISVMAX, NNDM,
     1        NPMAX
      PARAMETER ( NRMAX = 250, NTMAX = 250, NLEVM = 20,
     1            LHMAX = 160, NHMAX = 3000, ISVMAX = NRMAX*NHMAX,
     2            NNDM = 6, NPMAX = (LHMAX+1)*(LHMAX+2)/2 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      REAL    F( NRMAX, NTMAX )
      REAL    CONTL( NLEVM )
C
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     2        IWORK( NNDM )
C
      DOUBLE PRECISION XARR( NRMAX ), VEC( ISVMAX ),
     1                 WORK1( NNDM ), WORK2( NNDM ),
     2                 COEFM( NNDM, NNDM ), WORK3( NTMAX )
C
      DOUBLE PRECISION PA( NPMAX, NTMAX ), DPA( NPMAX, NTMAX )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C     .
C     . General
C     .
      INTEGER          NLEV, IW, LUIN, LU, IOP, IAS, NNDS, IRUN,
     1                 ISTAT, IDEV, INARR( 3 ), NH, ICOMP, ICONT,
     2                 LH, NR, NR1, I, ILEN, IVIEW, ICOMP1, ICOMP2
      INTEGER          ICONTOUR, ICONTOUR1, ICONTOUR2, RMFLAG,
     1                 RMFLAG1, RMFLAG2
      CHARACTER *(1)   CHDEV, ESCPCH
      CHARACTER *(180) STEM, LINE, FNAME
      REAL             PI, DELTA, ZERO, COORD,
     1                 VALMIN, VALMIN1, VALMIN2,
     2                 VALMAX, VALMAX1, VALMAX2
      DOUBLE PRECISION PHI, VMISS, THE, ZED, S
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
      PARAMETER ( PI=3.1415926, ZERO = 0.0, VMISS = 9.9d9 )
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
      PRINT *,' Enter NRAD, NTHE, NNDS.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NRAD, NTHE, NNDS
C
      IF ( NTHE.GT.NTMAX .OR. NRAD.GT.NRMAX .OR.
     1       NNDS.GT.NNDM     ) THEN
        PRINT *,' NNDS  = ', NNDS
        PRINT *,' NNDM  = ', NNDM
        PRINT *,' NTHE  = ', NTHE
        PRINT *,' NTMAX = ', NTMAX
        PRINT *,' NRAD  = ', NRAD
        PRINT *,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
c     PRINT *,' Enter XLEFT, XRIGHT, YBOT, YTOP'
c     CALL LINERD( LUIN, LINE, ESCPCH )
c     READ ( LINE, * ) XLEFT, XRIGHT, YBOT, YTOP
C
c     IF (  XLEFT.LT.0.0 .OR. XLEFT.GT.1.0 .OR. XLEFT.GE.XRIGHT 
c    1      .OR. XRIGHT.LT.0.0 .OR. XRIGHT.GT.1.0 .OR. YBOT.LT.0.0
c    2      .OR. YBOT.GT.1.0 .OR. YBOT.GE.YTOP .OR. YTOP.LT.0.0
c    3      .OR. YTOP.GT.1.0 )          THEN
c       PRINT *,' XLEFT  = ', XLEFT 
c       PRINT *,' XRIGHT = ', XRIGHT
c       PRINT *,' YBOT   = ', YBOT
c       PRINT *,' YTOP   = ', YTOP
c       PRINT *,' Program aborted.'
c       STOP
c     ENDIF
C
      PRINT *,' Enter TFIRST, TLAST, IVIEW, COORD'
      PRINT *,' iview = 1  --> meridian section.'
      PRINT *,'                COORD = const. phi value.'
      PRINT *,' iview = 2  --> equatorial section.'
      PRINT *,'                COORD = const. the value.'
      PRINT *,' iview = 3  --> constant z section.'
      PRINT *,'                COORD = const. z value.'
      PRINT *,' COORD will be multiplied by PI.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) TFIRST, TLAST, IVIEW, COORD
C
      IF ( IVIEW.NE.1 .AND. IVIEW.NE.2 .AND.
     1      IVIEW.NE.3  ) THEN
        PRINT *,' IVIEW = ', IVIEW
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      TFIRST  = TFIRST * PI
      TLAST   = TLAST  * PI
      IF ( IVIEW.EQ.1 .OR. IVIEW.EQ.2 )
     1                 COORD   = COORD  * PI
C
      IF ( IVIEW.EQ.1 ) PHI = DBLE( COORD )
      IF ( IVIEW.EQ.2 ) THE = DBLE( COORD )
      IF ( IVIEW.EQ.3 ) ZED = DBLE( COORD )
C
C Set inner and outer boundaries in simple cases:
C
      IF ( IVIEW.EQ.1 .OR. IVIEW.EQ.2 ) THEN
        RFIRST = REAL( XARR( 1 ) )
        RLAST  = REAL( XARR( NR ) )
        XARR( 1 ) = DBLE( REAL( XARR( 1 ) ) )
        XARR( NR ) = DBLE( REAL( XARR( NR ) ) )
      ENDIF
C
C Set inner and outer circles for constant z.
C
      IF ( IVIEW.EQ.3 ) THEN
C
C       First check that range of ZED supplied
C
        IF ( DABS( ZED ).GE.XARR( NR ) ) THEN
          PRINT *,' Outer radius  = ', XARR( NR )
          PRINT *,' ZED requested = ', ZED
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
C First do possibility of zed. gt. ri
C
        IF ( DABS( ZED ).GE.XARR( 1 )  ) THEN
          RFIRST = 0.01
          S      = DSQRT( XARR( NR )*XARR( NR ) - ZED*ZED )
          RLAST  = REAL( S )
        ELSE
          S      = DSQRT( XARR(  1 )*XARR(  1 ) - ZED*ZED )
          RFIRST = REAL( S )
          S      = DSQRT( XARR( NR )*XARR( NR ) - ZED*ZED )
          RLAST  = REAL( S )
        ENDIF
      ENDIF
C
C Enter number of contour levels
C
      PRINT *,' Enter NLEV, IDEV, ICONT, ICOMP1, ICOMP2, IAS '
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
      READ ( LINE, * ) NLEV, IDEV, ICONT, ICOMP1, ICOMP2, IAS
C
C Modification SJG Thu Jan 15 12:11:26 MET 2004
      IF ( ICOMP1.LT.0 ) THEN
        ICOMP1  = -ICOMP1
        RMFLAG1 = 1
      ELSE
        RMFLAG1 = 0
      ENDIF
      IF ( ICOMP2.LT.0 ) THEN
        ICOMP2  = -ICOMP2
        RMFLAG2 = 1
      ELSE
        RMFLAG2 = 0
      ENDIF
C END Modification SJG Thu Jan 15 12:11:26 MET 2004
C
C Make sure these are valid
C
      IF ( ICONT.LT.1 .OR. ICONT.GT.3 ) THEN
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
      READ ( LINE, * ) HUEPOS, HUENEG, CS, SCAL, IW
C
C Enter icontour1, valmin1, valmax1
C
      PRINT *,' Enter icontour1, valmin1, valmax1.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) ICONTOUR1, VALMIN1, VALMAX1
C
C Enter icontour2, valmin2, valmax2
C
      PRINT *,' Enter icontour2, valmin2, valmax2.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) ICONTOUR2, VALMIN2, VALMAX2
C
C This inputs a device from the user.
C
      CALL PGPLOT_FILE_INIT( ISTAT, IDEV, STEM, CHDEV )
      PRINT *,' ISTAT = ', ISTAT
      PRINT *,' CHDEV = ', CHDEV
C
C Call PGPAGE to clear or initialise a new page.
C
      CALL PGPAGE
C
      DO IRUN = 1, 2
C
      IF ( IRUN.EQ.1 ) THEN
        ICOMP    = ICOMP1
        RMFLAG   = RMFLAG1
        XLEFT    = 0.02
        XRIGHT   = 0.48
        YBOT     = 0.20
        YTOP     = 0.80
        ICONTOUR = ICONTOUR1
        VALMIN   = VALMIN1
        VALMAX   = VALMAX1
      ELSE
        ICOMP    = ICOMP2
        RMFLAG   = RMFLAG2
        XLEFT    = 0.52
        XRIGHT   = 0.98
        YBOT     = 0.20
        YTOP     = 0.80
        ICONTOUR = ICONTOUR2
        VALMIN   = VALMIN2
        VALMAX   = VALMAX2
      ENDIF
C
C Fill data array
C
      IOP = 0
      CALL RMATOP( F, ZERO, NRAD, NTHE, IOP )
C
      IF ( IVIEW.EQ.1 ) 
     1           CALL MERID_SEC_POLAR_EVAL( ICOMP, INARR, NNDS,
     2           IWORK, MHT, MHL, MHM, LH, IAS, PHI, VEC, XARR,
     3           WORK1, WORK2, COEFM, WORK3, PA, DPA, VMISS, F )
C
      IF ( IVIEW.EQ.2 ) 
     1           CALL EQ_SEC_POLAR_EVAL( ICOMP, INARR, NNDS,
     2           IWORK, MHT, MHL, MHM, LH, IAS, THE, VEC, XARR,
     3           WORK1, WORK2, COEFM, VMISS, F )
C
      IF ( IVIEW.EQ.3 )
     1           CALL CONSTANT_Z_POLAR_EVAL( ICOMP, INARR, NNDS,
     2           IWORK, MHT, MHL, MHM, LH, IAS, ZED, VEC, XARR,
     3           WORK1, WORK2, COEFM, VMISS, F )
C
      CALL ARR_RM_MEAN( NRAD, NTHE, F, RMFLAG )
C
      CALL PGSVP( XLEFT, XRIGHT, YBOT, YTOP )
C
      CALL EXTREME_COORDS_FIND( X1, Y1, X2, Y2 )
      PRINT *,' (X1,Y1) = (', X1,',',Y1,')'
      PRINT *,' (X2,Y2) = (', X2,',',Y2,')'
C
      CALL PGSWIN( X1, X2, Y1, Y2 )
C 
C Begin batch to buffer
C 
      CALL PGBBUF
C 
      DELTA = 0.00001
      IF ( ICONT.EQ.1 .OR. ICONT.EQ.3 )
     1      CALL ARRAY_POLAR_COLOUR_OPTION( F, NRAD, NTHE, NLEV,
     2                       CONTL, HUEPOS, HUENEG, CS, SCAL, DELTA,
     3                       ICONTOUR, VALMIN, VALMAX )
c    1      CALL ARRAY_POLAR_COLOUR( F, NRAD, NTHE, NLEV, CONTL,
c    2                         HUEPOS, HUENEG, CS, SCAL, DELTA )
C
      IF ( CHDEV.EQ.'P' ) THEN
        CALL PGSCI( 1 )
      ELSE
        CALL PGSCI( 1 )
      ENDIF
C
      IF ( ICONT.EQ.2 .OR. ICONT.EQ.3 )
     1  CALL ARRAY_POLAR_DRAW_OPTION( F, NRAD, NTHE, NLEV, CONTL,
     2                    DELTA, IW, ICONTOUR, VALMIN, VALMAX )
c    1      CALL ARRAY_POLAR_DRAW( F, NRAD, NTHE, NLEV,
c    2                             CONTL, DELTA, IW )
C
      CALL POLAR_BOX_DRAW
C
C End batch to buffer
C 
      CALL PGEBUF
C 
      ENDDO
C 
C Close the graphics devices.
C
      CALL PGCLOS
      STOP
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
