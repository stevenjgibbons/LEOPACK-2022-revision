C*********************************************************************
C double precision function Element of Banded Matrix Vector Mult. ****
C                           -          -      -      -      -     ****
C Steve Gibbons Tue Oct 24 17:06:42 MET DST 2000                     C
C____________________________________________________________________C
C                                                                    C
C A is a banded matrix stored in LAPACK format (i.e. element a_{ij}  C
C is in A( ku + 1 + i - j, j ) ...). X is a vector of length n.      C
C A is n x n with ku upper diagonals and kl lower diagonals.         C
C ebmvm returns the value of y( i ) where y:= A . x                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     I         : Element of y vector required.                      C
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
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION EBMVM( I, N, LDA, KL, KU, A, X )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER I, N, LDA, KL, KU
      DOUBLE PRECISION A( LDA, * ), X( * ), EBMVM
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER J, K
      DOUBLE PRECISION ZERO, TEMP
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      EBMVM = ZERO
C
      DO J = MAX( 1, I - KL ), MIN( N, I + KU )
        TEMP = X( J )
        IF ( TEMP.NE.ZERO ) THEN
          K = KU + 1 - J + I
          EBMVM = EBMVM + TEMP*A( K, J )
        ENDIF
      ENDDO
C
      RETURN
      END
C*********************************************************************
