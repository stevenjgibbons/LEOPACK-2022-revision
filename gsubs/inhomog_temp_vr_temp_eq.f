C*********************************************************************
C subroutine INHOMOG TEMP VR and TEMP EQuatorial section plot ********
C            ------- ---- --     ---- --                      ********
C Steve Gibbons Mon May  7 12:41:48 WEST 2001                        C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NLEV      : Number of contour levels.                          C
C     NPHIP     : Number of phi points for plot.                     C
C     NRADP     : Number of rad points for plot.                     C
C     NNDS      : Number of nodes for interpolation.                 C
C                                                                    C
C     INARR     : Indexing for homogeneous solution vector. Dim (3)  C
C     MHT       : Harmonic type: Dim ( ) (Homogeneous vector)        C
C     MHL       : Harmonic degree: Dim ( ) (Homogeneous vector)      C
C     MHM       : Harmonic order: Dim ( ) (Homogeneous vector)       C
C                                                                    C
C     INARRC    : Indexing for inhomog. vector. Dim (3)              C
C     MHTC      : Harmonic type: Dim ( ) (Inhomogeneous vector)      C
C     MHLC      : Harmonic degree: Dim ( ) (Inhomogeneous vector)    C
C     MHMC      : Harmonic order: Dim ( ) (Inhomogeneous vector)     C
C                                                                    C
C     NPBA      : Number of grid points between arrows               C
C     IDEV      : 1 = .gif, 5 = .ps                                  C
C     ITHICK    : Thickness of line.                                 C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Homogeneous part of vector. Dim ( )                C
C     VECC      : Inhomogeneous part of vector. Dim ( )              C
C     XARR      : Radial spacings array. Dim ( )                     C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
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
C     RLONG     : Length of the longest arrow.                       C
C     RHEAD     : Size of the longest arrow-head.                    C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     STEM      : Filename stem (terminate with ' ')                 C
C                                                                    C
C*********************************************************************
      SUBROUTINE INHOMOG_TEMP_VR_TEMP_EQ( NLEV, NPHIP, NRADP, NNDS,
     1         INARR, MHT, MHL, MHM, INARRC, MHTC, MHLC, MHMC, NPBA,
     2         IDEV, ITHICK, VEC, VECC, XARR, PAPWIDTH, HUEPOS,
     3         HUENEG, CS, SCAL, RLONG, RHEAD, STEM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NLEV, NPHIP, NRADP, NNDS, INARR( 3 ), MHT( * ),
     1        MHL( * ), MHM( * ), INARRC( 3 ), MHTC( * ), MHLC( * ),
     2        MHMC( * ), NPBA, IDEV, ITHICK
      DOUBLE PRECISION VEC( * ), VECC( * ), XARR( * )
      REAL    PAPWIDTH, HUEPOS, HUENEG, CS, SCAL, RLONG, RHEAD
      CHARACTER *(*) STEM
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER  NRMAX, NTMAX, NLEVM, NNDM, NR, ISTAT, IOP,
     1         ICOMP, IAS, LH
      PARAMETER ( NRMAX = 250, NTMAX = 250, NLEVM = 20,
     1            NNDM = 6, IAS = 0 )
      REAL     VRADA( NRMAX, NTMAX ),
     1         VPHIA( NRMAX, NTMAX ),
     2         TEMPA( NRMAX, NTMAX ),
     3         CONTL( NLEVM )
      REAL     PI, ZERO, RATIO, VMISS
      PARAMETER ( PI = 3.1415926, ZERO = 0.0, VMISS = 1.0e9 )
      CHARACTER *(1) CHDEV
C
      INTEGER          IWORK( NNDM )
      DOUBLE PRECISION WORK1( NNDM ), WORK2( NNDM ),
     1                 COEFM( NNDM, NNDM ), THE, DPI
      PARAMETER ( DPI = 3.1415926d0 )
C
      REAL XLEFT, XRIGHT, YBOT, YTOP, X1, X2, Y1, Y2, DELTA,
     1     X, Y, FJUST, ANGLE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
C_____________________________________________________________________
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check dimensions of input variables
C
      IF ( NLEV.GT.NLEVM .OR. NRADP.GT.NRMAX .OR.
     1     NPHIP.GT.NTMAX .OR. NNDS.GT.NNDM ) THEN
        PRINT *,' Subroutine INHOMOG_TEMP_VR_TEMP_EQ'
        PRINT *,' NLEV   = ',NLEV ,' NLEVM  = ', NLEVM
        PRINT *,' NRADP  = ',NRADP,' NRMAX  = ', NRMAX
        PRINT *,' NPHIP  = ',NPHIP,' NTMAX  = ', NTMAX
        PRINT *,' NNDS   = ',NNDS ,' NNDM   = ', NNDM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Our input values seem o.k. so let's proceed
C
      LH     = 3000
C
C Note that EQ_SEC_POLAR_EVAL doesn't really need to know LH
C so we just give LH a junk value ...
C
      NR     = INARR( 2 )
      NRAD   = NRADP
      NTHE   = NPHIP
      RFIRST = REAL( XARR(  1 ) )
      RLAST  = REAL( XARR( NR ) )
      TFIRST = 0.0
      TLAST  = 2.0*PI
C
      RATIO = 0.5
      CALL PGPLOT_FILE_INIT( ISTAT, IDEV, STEM, CHDEV )
      CALL PGPAP( PAPWIDTH, RATIO )
      CALL PGPAGE
      CALL PGSLW( ITHICK )
C
C Fill data array for colouring:
C
      IOP = 0
      CALL RMATOP( VRADA, ZERO, NRAD, NTHE, IOP )
      CALL RMATOP( VPHIA, ZERO, NRAD, NTHE, IOP )
      CALL RMATOP( TEMPA, ZERO, NRAD, NTHE, IOP )
C
      THE   = 0.5d0*DPI
      ICOMP = 1
      CALL EQ_SEC_POLAR_EVAL( ICOMP, INARR, NNDS, IWORK, MHT, MHL,
     1                        MHM, LH, IAS, THE, VEC, XARR,
     3                        WORK1, WORK2, COEFM, VMISS, VRADA )
      ICOMP = 3
      CALL EQ_SEC_POLAR_EVAL( ICOMP, INARR, NNDS, IWORK, MHT, MHL,
     1                        MHM, LH, IAS, THE, VEC, XARR,
     3                        WORK1, WORK2, COEFM, VMISS, VPHIA )
      ICOMP = 7
      CALL EQ_SEC_POLAR_EVAL( ICOMP, INARR, NNDS, IWORK, MHT, MHL,
     1                        MHM, LH, IAS, THE, VEC, XARR,
     3                        WORK1, WORK2, COEFM, VMISS, TEMPA )
      CALL EQ_SEC_POLAR_EVAL( ICOMP, INARRC, NNDS, IWORK, MHTC, MHLC,
     1                        MHMC, LH, IAS, THE, VECC, XARR,
     3                        WORK1, WORK2, COEFM, VMISS, TEMPA )
C
C Left hand plot is contours of V_r
C
      XLEFT   = 0.05
      XRIGHT  = 0.45
      YBOT    = 0.05
      YTOP    = 0.85
      X1 = RLAST*(-1.0)
      X2 = RLAST
      Y1 = RLAST*(-1.0)
      Y2 = RLAST
      CALL PGSVP( XLEFT, XRIGHT, YBOT, YTOP )
      CALL PGSWIN( X1, X2, Y1, Y2 )
      CALL PGBBUF
C
      DELTA = 0.01
      CALL ARRAY_POLAR_COLOUR( VRADA, NRAD, NTHE, NLEV, CONTL,
     1                         HUEPOS, HUENEG, CS, SCAL, DELTA )
      CALL PGSCI( 1 )
      CALL POLAR_BOX_DRAW
C
      CALL PGEBUF
      XLEFT   = 0.05
      XRIGHT  = 0.45
      YBOT    = 0.90
      YTOP    = 0.95
      X1 = 0.0
      X2 = 1.0
      Y1 = 0.0
      Y2 = 1.0
      CALL PGSVP( XLEFT, XRIGHT, YBOT, YTOP )
      CALL PGSWIN( X1, X2, Y1, Y2 )
      CALL PGBBUF
C
      CALL PGSCH( 1.5 )
      X      = 0.5
      Y      = 0.5
      ANGLE  = 0.0
      FJUST  = 0.5
      CALL PGPTXT( X, Y, ANGLE, FJUST, 'Radial velocity' )
C
      CALL PGEBUF
C
      IF ( CHDEV.EQ.'P' ) THEN
        CALL PGSCI( 1 )
      ELSE
        CALL PGSCI( 0 )
      ENDIF
C
C Right hand plot is temperature contours with lines
C of flow:
C 
      XLEFT   = 0.55
      XRIGHT  = 0.95
      YBOT    = 0.05
      YTOP    = 0.85
      X1 = RLAST*(-1.0)
      X2 = RLAST
      Y1 = RLAST*(-1.0)
      Y2 = RLAST
      CALL PGSVP( XLEFT, XRIGHT, YBOT, YTOP )
      CALL PGSWIN( X1, X2, Y1, Y2 )
      CALL PGBBUF
C
      DELTA = 0.01
      CALL ARRAY_POLAR_COLOUR( TEMPA, NRAD, NTHE, NLEV, CONTL,
     1                         HUEPOS, HUENEG, CS, SCAL, DELTA )
      CALL ARRAY_POLAR_ARROWS( VRADA, VPHIA, NPBA, RLONG, RHEAD,
     1                         ITHICK )
C
      CALL PGSCI( 1 )
      CALL POLAR_BOX_DRAW
C
      CALL PGEBUF
      XLEFT   = 0.55
      XRIGHT  = 0.95
      YBOT    = 0.90
      YTOP    = 0.95
      X1 = 0.0
      X2 = 1.0
      Y1 = 0.0
      Y2 = 1.0
      CALL PGSVP( XLEFT, XRIGHT, YBOT, YTOP )
      CALL PGSWIN( X1, X2, Y1, Y2 )
      CALL PGBBUF
C
      CALL PGSCH( 1.5 )
      X      = 0.5
      Y      = 0.5
      ANGLE  = 0.0
      FJUST  = 0.5
      CALL PGPTXT( X, Y, ANGLE, FJUST, 'Temperature and flow' )
C
      CALL PGEBUF
C
C Close the graphics devices.
C
      CALL PGCLOS
C
      RETURN
      END
C*********************************************************************
