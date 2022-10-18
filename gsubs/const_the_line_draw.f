C*********************************************************************
C subroutine CONST_THE_LINE_DRAW *************************************
C            ------------------- *************************************
C Steve Gibbons Fri Jan 19 12:54:31 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Simply draws an line at theta = THE from the first radial point    C
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
      SUBROUTINE CONST_THE_LINE_DRAW( THE )
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      REAL    THE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRAD
      REAL    XWORLD, YWORLD, RAD
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameter:-
C
c     IF ( THE.LT.TFIRST .OR. THE.GT.TLAST ) THEN
c       PRINT *,' Subroutine CONST_THE_LINE_DRAW.'
c       PRINT *,' THE    = ', THE
c       PRINT *,' TFIRST = ', TFIRST
c       PRINT *,' TLAST  = ', TLAST
c       PRINT *,' Program aborted.'
c       STOP
c     ENDIF
C
      XWORLD = RFIRST*COS( THE )
      YWORLD = RFIRST*SIN( THE )
C
      CALL PGMOVE( XWORLD, YWORLD )
C
      DO IRAD = 1, NRAD
        RAD = RFIRST +
     1          (RLAST-RFIRST)*REAL(IRAD - 1)/REAL( NRAD - 1 )
        XWORLD = RAD*COS( THE )
        YWORLD = RAD*SIN( THE )
        CALL PGDRAW( XWORLD, YWORLD )
      ENDDO
C
      RETURN
      END
C*********************************************************************
