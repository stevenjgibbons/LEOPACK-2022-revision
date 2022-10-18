C*********************************************************************
C subroutine Banded Matrix Addition 2 Banded Matrix ******************
C            -      -      -        - -      -      ******************
C Steve Gibbons Fri Dec 10 12:01:28 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C A pretty specialised routine which trivially adds the double       C
C precision multiple FACB * the matrix B to FACA * the matrix A.     C
C Both A and B have 2*KL + 1 width but A has an additional KLE       C
C subdiagonals for the benefit of LAPACK routines.                   C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     KL        : Number of lower diagonals in banded matrix.        C
C     N2        : Second dimension of matrices A and B.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Double precision matrix dim ( 3*KL + 1, N2 )       C
C     B         : Double precision matrix dim ( 2*KL + 1, N2 )       C
C     FACA      : Multiplication factor                              C
C     FACB      : Multiplication factor                              C
C                                                                    C
C  A( i + kl, j ) := FACA * A( i + kl, j )  +  FACB * B( i, j )      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE BMA2BM( KL, N2, A, B, FACA, FACB )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KL, N2
      DOUBLE PRECISION A( 3*KL + 1, N2 ), B( 2*KL + 1, N2 ),
     1                 FACA, FACB
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, J, N1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      N1 = 2*KL + 1
      DO J = 1, N2
        DO I = 1, N1
          A( I + KL, J ) = FACA*A( I + KL, J ) + FACB*B( I, J )
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
