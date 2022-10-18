C*********************************************************************
C subroutine SHELL sort version 1. Based on numerical Recipes shell **
C            -----              -                             ^^^^^ **
C Adapted by Steve Gibbons Wed Jun 14 11:43:31 BST 2000              C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N         : Length of array to be sorted.                      C
C     IFLAG     :   1 --> sort by value                              C
C                   2 --> sort by magnitude                          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     A         : Array of dimension N to be sorted.                 C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHELL1( N, A, IFLAG )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N, IFLAG
      DOUBLE PRECISION A( N )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER i,j,inc
      DOUBLE PRECISION V
      LOGICAL OMAG
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IFLAG.NE.1 .AND. IFLAG.NE.2 ) THEN
         PRINT *,' Subroutine SHELL1 '
         PRINT *,' IFLAG must be set to 1 or 2.'
         PRINT *,' Program aborted,'
         STOP
      ENDIF
      IF ( IFLAG.EQ.1 ) OMAG = .FALSE.
      IF ( IFLAG.EQ.2 ) OMAG = .TRUE.
C-------------------------------------- numerical recipes part
      INC = 1
 1    INC = 3*INC + 1
      IF ( INC.LE.N ) GOTO 1
 2    CONTINUE
        INC = INC/3
        DO 11 I = INC + 1, N
          V = A( I )
          J = I
 3        if( OMAG .AND. ( ABS( A(J-INC) ) .GT. ABS(V) )  .OR.
     1        .NOT. OMAG .AND. ( A(J-INC) .GT. V )     ) THEN
            A( J ) = A( J-INC )
            J = J - INC
            if( J.LE.INC ) GOTO 4
          GOTO 3
          ENDIF
 4        A(J) = V
 11     CONTINUE
      IF ( INC.GT.1 ) GOTO 2
      RETURN
      END
C*****************************************************************

