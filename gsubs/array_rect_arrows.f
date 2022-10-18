C*********************************************************************
C subroutine ARRAY RECT ARROWS **************************************
C            ----- ----- ------ **************************************
C Steve Gibbons Wed Jan 24 13:17:34 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Requires COMMON block / PARAMC / THETA1, THETA2, PHI1, PHI2        C
C                                                                    C
C Draws arrows upon a plot in polar coordinates. The arrays VC1 and  C
C VC2 (both with dim {NPHI, NTHE}) contain the r and theta vel.s     C
C respectively.                                                      C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NPHI      : Number of phi grid nodes.                          C
C     NTHE      : Number of theta grid nodes.                        C
C     NPBA      : Number of points between arrows.                   C
C     IW        : Width of line.                                     C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     VC1       : Dim ( NPHI, NTHE ). Radial velocity component.     C
C     VC2       : Dim ( NPHI, NTHE ). Theta velocity component.      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ARRAY_RECT_ARROWS( VC1, VC2, NPBA, RLONG, RHEAD,
     1                               IW, NPHI, NTHE )
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      REAL    THETA1, THETA2, PHI1, PHI2
      COMMON  / PARAMC / THETA1, THETA2, PHI1, PHI2
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NPBA, IW, NPHI, NTHE
      REAL             VC1( NPHI, NTHE ), VC2( NPHI, NTHE ),
     1                 RLONG, RHEAD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IPHI, ITHE
      REAL             VPHI, VTHE, VRMAX, VTMAX,
     1                 XSTART, XEND, YSTART, YEND,
     2                 P1, T1, FACR, FACT, RLOW,
     3                 PEND, TEND, RMAG, MGMAX, VMAG
      PARAMETER ( RLOW = 0.000001 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NPBA.LT.1 ) THEN
        PRINT *,' Subroutine ARRAY_RECT_ARROWS.'
        PRINT *,' NPBA = ', NPBA
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C First, we calculate the maximum values of
C VPHI and VTHE (not taking sign into account)
C
      VRMAX = ABS( VC1( 8, 8 )  )
      VTMAX = ABS( VC2( 1, 1 )  )
      MGMAX = SQRT( VRMAX*VRMAX + VTMAX*VTMAX )
C
      DO ITHE = 1, NTHE
        DO IPHI = 1, NPHI
          VPHI = VC1( IPHI, ITHE )
          VTHE = VC2( IPHI, ITHE )
          VMAG = SQRT( VPHI*VPHI + VTHE*VTHE )
          IF ( ABS( VPHI ).GT.VRMAX ) VRMAX = ABS( VPHI )
          IF ( ABS( VTHE ).GT.VTMAX ) VTMAX = ABS( VTHE )
          IF ( ABS( VMAG ).GT.MGMAX ) MGMAX = ABS( VMAG )
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
        T1 = THETA1 + (THETA2-THETA1)*REAL(ITHE-1)/REAL(NTHE-1)
        DO IPHI = 1, NPHI - 1
          P1 = PHI1 + (PHI2-PHI1)*REAL(IPHI-1)/REAL(NPHI-1)
          IF ( IPHI/NPBA*NPBA.NE.IPHI ) GOTO 51
          VPHI = VC1( IPHI, ITHE )
          VTHE = VC2( IPHI, ITHE )
          VMAG = SQRT( VPHI*VPHI + VTHE*VTHE )
C         .
C         . OK - so we are going to put an arrow
C         . at this point.
C         .
          XSTART = P1
          YSTART = T1
C         .
          PEND   = P1 + VPHI*FACR*(PHI2-PHI1)/ABS(PHI2-PHI1)
          TEND   = T1 + VTHE*FACT*(THETA2-THETA1)/ABS(THETA2-THETA1)
C         .
          XEND = PEND
          YEND = TEND
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
