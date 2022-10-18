C*********************************************************************
C subroutine ARRAY POLAR ARROWS **************************************
C            ----- ----- ------ **************************************
C Steve Gibbons Wed Jan 24 13:17:34 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Draws arrows upon a plot in polar coordinates. The arrays VC1 and  C
C VC2 (both with dim {NRAD, NTHE}) contain the r and theta vel.s     C
C respectively.                                                      C
C The coordinates are defined by a COMMON block:                     C
C                                                                    C
C      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST   C
C                                                                    C
C It consists of two INTEGERs                                        C
C NRAD (number of equally spaced radial grid nodes) and              C
C NTHE (number of equally spaced phi grid nodes), and four REAL      C
C variables RFIRST, RLAST, TFIRST, TLAST                             C
C                                                                    C
C RFIRST is the inner boundary radius                                C
C RLAST is the outer boundary radius                                 C
C                                                                    C
C TFIRST is the first phi point, in radians.                         C
C TLAST is the last phi point, in radians.                           C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NPBA      : Number of points between arrows.                   C
C     IW        : Width of line.                                     C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     VC1       : Dim ( NRAD, NTHE ). Radial velocity component.     C
C     VC2       : Dim ( NRAD, NTHE ). Theta velocity component.      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ARRAY_POLAR_ARROWS( VC1, VC2, NPBA, RLONG, RHEAD,
     1                               IW )
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NPBA, IW
      REAL             VC1( NRAD, NTHE ), VC2( NRAD, NTHE ),
     1                 RLONG, RHEAD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IRAD, ITHE
      REAL             VRAD, VTHE, VRMAX, VTMAX,
     1                 XSTART, XEND, YSTART, YEND,
     2                 R1, T1, FACR, FACT, RLOW, RLIM,
     3                 REND, TEND, RMAG, MGMAX, VMAG
      PARAMETER ( RLOW = 0.000001, RLIM = 1.0e8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NPBA.LT.1 ) THEN
        PRINT *,' Subroutine ARRAY_POLAR_ARROWS.'
        PRINT *,' NPBA = ', NPBA
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C First, we calculate the maximum values of
C VRAD and VTHE (not taking sign into account)
C
      VRMAX = 0.0
      VTMAX = 0.0
      MGMAX = 0.0
C
      DO ITHE = 1, NTHE
        DO IRAD = 1, NRAD
          VRAD = VC1( IRAD, ITHE )
          VTHE = VC2( IRAD, ITHE )
          VMAG = SQRT( VRAD*VRAD + VTHE*VTHE )
          IF ( ABS( VRAD ).GT.VRMAX .AND.
     1         ABS( VRAD ).LT.RLIM   ) VRMAX = ABS( VRAD )
          IF ( ABS( VTHE ).GT.VTMAX.AND.
     1         ABS( VTHE ).LT.RLIM    ) VTMAX = ABS( VTHE )
          IF ( ABS( VMAG ).GT.MGMAX.AND.
     1         ABS( VMAG ).LT.RLIM    ) MGMAX = ABS( VMAG )
        ENDDO
      ENDDO
C
      IF (  VRMAX.LE.RLOW  ) THEN
        FACR = 0.0
      ELSE
        FACR = RLONG/VRMAX
      ENDIF
C
      IF (  VTMAX.LE.RLOW  ) THEN
        FACT = 0.0
      ELSE
        FACT = RLONG/VTMAX
      ENDIF
C
      DO ITHE = 1, NTHE - 1
        IF ( ITHE/NPBA*NPBA.NE.ITHE ) GOTO 50
        T1 = TFIRST + (TLAST-TFIRST)*REAL(ITHE-1)/REAL(NTHE-1)
        DO IRAD = 1, NRAD - 1
          R1 = RFIRST + (RLAST-RFIRST)*REAL(IRAD-1)/REAL(NRAD-1)
          IF ( IRAD/NPBA*NPBA.NE.IRAD ) GOTO 51
          VRAD = VC1( IRAD, ITHE )
          VTHE = VC2( IRAD, ITHE )
          VMAG = SQRT( VRAD*VRAD + VTHE*VTHE )
C         .
C         . OK - so we are going to put an arrow
C         . at this point.
C         .
          XSTART = R1*COS( T1 )
          YSTART = R1*SIN( T1 )
C         .
          REND   = R1 + VRAD*FACR
          TEND   = T1 + VTHE*FACT
C         .
          XEND = REND*COS( TEND )
          YEND = REND*SIN( TEND )
C         .
          RMAG = RHEAD*VMAG/MGMAX
          CALL PGSCH( RMAG )
          CALL PGSLW( IW )
C         .
          CALL PGARRO( XSTART, YSTART, XEND, YEND )
C         .
 51       CONTINUE
        ENDDO
 50     CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************
