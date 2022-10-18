C*********************************************************************
C subroutine Euler Angles of Rotation Matrix Inversion Routine *******
C            -     -         -        -      -         -       *******
C Steve Gibbons Sat Mar 17 15:38:29 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C All this does is invert a double precision 3 x 3 matrix!!          C
C It is a front end to some rather irritating LAPACK routines for    C
C such a simple task!                                                C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C                                                                    C
C     EARCM     : Euler angle matrix Dim( 3, 3 ).                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE EARMIR( EARCM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION EARCM( 3, 3 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          M, N, LDA, IPIV( 3 ), INFO, LWORK
      DOUBLE PRECISION WORK( 3 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      M      = 3
      N      = 3
      LDA    = 3
      LWORK  = 3
C
      CALL DGETRF( M, N, EARCM, LDA, IPIV, INFO )
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine EARMIR.'
        PRINT *,' DGETRF returned INFO = ', INFO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      CALL DGETRI( N, EARCM, LDA, IPIV, WORK, LWORK, INFO )
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine EARMIR.'
        PRINT *,' DGETRI returned INFO = ', INFO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
