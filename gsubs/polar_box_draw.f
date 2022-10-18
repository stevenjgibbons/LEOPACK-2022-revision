C*********************************************************************
C subroutine POLAR_BOX_DRAW ******************************************
C            ----- --- ---- ******************************************
C Steve Gibbons Mon Jan 22 11:59:27 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Draws inner and outer boundaries and (if necessary) bounding lines C
C of constant theta.                                                 C
C                                                                    C
C The following COMMON block must be declared.                       C
C                                                                    C
C      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST   C
C                                                                    C
C It consists of two INTEGERs                                        C
C NRAD (number of equally spaced radial grid nodes) and              C
C NTHE (number of equally spaced theta grid nodes), and four REAL    C
C variables RFIRST, RLAST, TFIRST, TLAST                             C
C                                                                    C
C RFIRST is the inner boundary radius                                C
C RLAST is the outer boundary radius                                 C
C                                                                    C
C TFIRST is the first theta point, in radians.                       C
C TLAST is the last theta point, in radians.                         C
C                                                                    C
C*********************************************************************
      SUBROUTINE POLAR_BOX_DRAW
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFAC
      REAL    PI, TWOPI, RLOW, ZERO, DIFF, FAC, FAC2
      LOGICAL LINEDRAW
      PARAMETER ( PI=3.1415926, RLOW = 0.001, ZERO = 0.0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      TWOPI = 2*PI
      DIFF  = TLAST-TFIRST
      FAC   = DIFF/TWOPI
      IFAC  = INT( FAC )
      FAC2  = REAL( IFAC )
      IF ( ABS( FAC - FAC2 ).LT.RLOW    .OR.
     1     ABS( FAC - FAC2 - 1.0).LT.RLOW  ) THEN
        LINEDRAW = .FALSE.
      ELSE
        LINEDRAW = .TRUE.
      ENDIF
C
      IF ( LINEDRAW ) THEN
        CALL CONST_THE_LINE_DRAW( TFIRST )
        CALL CONST_THE_LINE_DRAW( TLAST  )
      ENDIF
C
C Draw circle lines
C
      IF ( RFIRST.GT.ZERO ) CALL CONST_RAD_ARC_DRAW( RFIRST )
      CALL CONST_RAD_ARC_DRAW( RLAST )
C
      RETURN
      END
C*********************************************************************
