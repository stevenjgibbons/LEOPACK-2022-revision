C*********************************************************************
C subroutine Simple Banded Matrix Vector Multiplication Routine ******
C            -      -      -      -      -              -       ******
C Steve Gibbons Tue Oct 24 15:24:28 MET DST 2000                     C
C____________________________________________________________________C
C                                                                    C
C This routine is basically a specialisation of the BLAS level 2     C
C routine DGBMV. It multiplies a banded matrix A, stored in LAPACK   C
C format (see MATIND) by a vector X to give a vector Y.              C
C                                                                    C
C A has dimensions ( LDA, N ), KL lower diagonals and KU upper       C
C diagonals. A is n x n in all cases. If more general, use DGBMV!    C
C                                                                    C
C   y := A.x                                                         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N         : Main dimension of matrix A.                        C
C     LDA       : Leading dimension of A.                            C
C     KL        : Number of lower diagonals in matrix A.             C
C     KU        : Number of upper diagonals in matrix A.             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix dimension ( LDA, * )                        C
C                 Element a_{ij} is stored in                        C
C                 A( ku + 1 + i - j, j )                             C
C                                                                    C
C     X         : Vector of length n.                                C
C     Y         : Vector of length n.                                C
C                                                                    C
C____________________________________________________________________C
C 
C*********************************************************************
      SUBROUTINE SBMVMR( N, LDA, KL, KU, A, X, Y )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N, LDA, KL, KU
      DOUBLE PRECISION A( LDA, * ), X( * ), Y( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER KUP1, K, I, J
      DOUBLE PRECISION ZERO, TEMP
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First zero the vector Y
C
      DO I = 1, N
        Y( I ) = ZERO
      ENDDO
C
      KUP1 = KU + 1
      DO J = 1, N
        IF ( X( J ).NE.ZERO ) THEN
          TEMP = X( J )
          K    = KUP1 - J
          DO I = MAX( 1, J - KU ), MIN( N, J + KL )
            Y( I ) = Y( I ) + TEMP*A( K + I, J )
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
      END
C*********************************************************************
