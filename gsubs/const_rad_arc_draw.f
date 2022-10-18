C*********************************************************************
C subroutine CONST_RAD_ARC_DRAW **************************************
C            ------------------ **************************************
C Steve Gibbons Fri Jan 19 12:54:31 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Simply draws an arc at radius RAD from the first theta point       C
C to the last.                                                       C
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
C                                                                    C
C*********************************************************************
      SUBROUTINE CONST_RAD_ARC_DRAW( RAD )
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      REAL    RAD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHE
      REAL    XWORLD, YWORLD, THE
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameter:-
C
      IF ( RAD.LT.RFIRST .OR. RAD.GT.RLAST ) THEN
        PRINT *,' Subroutine CONST_RAD_ARC_DRAW.'
        PRINT *,' RAD    = ', RAD
        PRINT *,' RFIRST = ', RFIRST
        PRINT *,' RLAST  = ', RLAST
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      XWORLD = RAD*COS( TFIRST )
      YWORLD = RAD*SIN( TFIRST )
C
      CALL PGMOVE( XWORLD, YWORLD )
C
      DO ITHE = 1, NTHE
        THE = TFIRST +
     1          (TLAST-TFIRST)*REAL(ITHE - 1)/REAL( NTHE - 1 )
        XWORLD = RAD*COS( THE )
        YWORLD = RAD*SIN( THE )
        CALL PGDRAW( XWORLD, YWORLD )
      ENDDO
C
      RETURN
      END
C*********************************************************************
