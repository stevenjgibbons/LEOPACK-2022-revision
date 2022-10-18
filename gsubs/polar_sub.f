C*********************************************************************
C subroutine POLAR_SUB ***********************************************
C            --------- ***********************************************
C Steve Gibbons Fri Jan 19 12:35:10 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C This is simply to assist the routine PGCONX.                       C
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
C Note that Z is the value of the contour passed in by PGCONX.       C
C It is not referred to by POLAR_SUB.                                C
C                                                                    C
C*********************************************************************
      SUBROUTINE POLAR_SUB( VISIBLE, RAD, THE, Z )
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER VISIBLE
      REAL RAD, THE, Z
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      REAL XWORLD, YWORLD, R1, T1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      R1 = RFIRST +
     1     (RLAST-RFIRST)*(RAD - 1.0)/REAL( NRAD - 1 )
      T1 = TFIRST +
     1     (TLAST-TFIRST)*(THE - 1.0)/REAL( NTHE - 1 )
C
      XWORLD = R1*COS( T1 )
      YWORLD = R1*SIN( T1 )
C
      IF ( VISIBLE.EQ.0 ) THEN
        CALL PGMOVE( XWORLD, YWORLD )
      ELSE
        CALL PGDRAW( XWORLD, YWORLD )
      ENDIF
C
      RETURN
      END
C*********************************************************************
