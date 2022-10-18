C*********************************************************************
C subroutine EXTREME CO-ORDinateS FIND *******************************
C            ------- -- ---     - ---- *******************************
C Steve Gibbons Mon Jan 22 18:06:01 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Calculates the necessary mappings for the inclusion of every point C
C for given rad and theta possibilities. Returns co-ords (X1, Y1 )   C
C for the bottom left hand corner and (X2, Y2) for the top right     C
C hand corner of a box.                                              C
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
      SUBROUTINE EXTREME_COORDS_FIND( X1, Y1, X2, Y2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      REAL    X1, Y1, X2, Y2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRAD, ITHE
      REAL    T1, R1, X, Y
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      X1 =  1000.0
      Y1 =  1000.0
      X2 = -1000.0
      Y2 = -1000.0
C
      DO ITHE = 1, NTHE
        T1 = TFIRST + (TLAST-TFIRST)*REAL(ITHE-1)/REAL( NTHE - 1 )
        DO IRAD = 1, NRAD
          R1 = RFIRST + (RLAST-RFIRST)*REAL(IRAD-1)/REAL( NRAD - 1 )
          X  = R1*COS( T1 )
          Y  = R1*SIN( T1 )
          IF ( X.LT.X1 ) X1 = X
          IF ( X.GT.X2 ) X2 = X
          IF ( Y.LT.Y1 ) Y1 = Y
          IF ( Y.GT.Y2 ) Y2 = Y
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
