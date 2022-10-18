C*********************************************************************
C subroutine Double Precision Matrix Boundary Value Fix **************
C            -      -         -      -        -     -   **************
C Steve Gibbons 22.3.99                                              C
C____________________________________________________________________C
C Given a harmonic IH and a radial grid node IRN, dpmbvf will        C
C replace the corresponding row of the matrix with an equation which C
C enforces a boundary condition by relating the values of the        C
C function at adjacent grid points.                                  C
C                                                                    C
C All entries in the corresponding row are filled with 0.0d0         C
C                                                                    C
C Then, for j = 1, NNT, the matrix element corresponding to the      C
C row of harmonic IH, grid node IRN and the column of harmonic IH,   C
C grid node JRN( j ) is filled with the value VALS( j )              C
C  Note that in cases of banded matrices ( imf=1,2 ), care must      C
C  be taken with the values of JRN. If they are out of bounds,       C
C  the routine MATIND will trap the error.                           C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of harmonics (all types)                    C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C                                                                    C
C     IMF       : Matrix format flag.                                C
C                                                                    C
C          imf = 1; Matrix is in LAPACK banded format                C
C                   ie element a_{i,j} is stored in                  C
C                   A( kle + ku + 1 + i - j , j )                    C
C                                                                    C
C          imf = 2; Matrix is banded but with element a_{i,j}        C
C                   stored in A( kl + 1 + j - i , i ).               C
C                                                                    C
C          imf = 3; Matrix is square - ie a_{i,j} is stored          C
C                   in A( i, j ).                                    C
C     NNT       : Number of new terms                                C
C     JRN       : Array dimension ( * ) - see above                  C
C     IRN       : Grid node at which alteration is to be done        C
C     IH        : Number of harmonic to be altered                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     AMAT      : Matrix. Dimensions ( N1, N2 )                      C
C              Will generally be banded due to the nature of the     C
C              numerical scheme. KL, KU and KLE parameterise this.   C
C     VALS      : Array *(*) - see above for information.            C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT, 
     1                    IRN, IH, JRN, AMAT, VALS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, NR, NH, KL, KU, KLE, IMF, NNT, 
     1        JRN( * ), IRN, IH
      DOUBLE PRECISION AMAT( N1, N2 ), VALS( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR, IC, IROW, ICOL, J
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of N2 against NH and NR etc.
C
      IF ( N2.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPMBVF, bad array size'
        PRINT *,' N2 (second array dimension) = ',N2
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C First zero row of matrix. The procedure is different
C depending on the IMF flag
C
      IR = ( IRN - 1 )*NH + IH
C
      IF ( IMF.EQ.1 ) THEN
        DO IC = IR - KL, IR + KU
          IF ( IC.GT.N2 .OR. IC.LT.1 ) GOTO 500
          IROW = KLE + KU + 1 + IR - IC
          AMAT( IROW, IC ) = ZERO
 500    CONTINUE
        ENDDO
      ENDIF
C
      IF ( IMF.EQ.2 ) THEN
        DO IC = IR - KL, IR + KU
          IF ( IC.GT.N2 .OR. IC.LT.1 ) GOTO 501
          ICOL = IR
          IROW = KL + 1 + IC - IR
          AMAT( IROW, ICOL ) = ZERO
 501    CONTINUE
        ENDDO
      ENDIF
C
      IF ( IMF.EQ.3 ) THEN
        DO IC = 1, N2
          AMAT( IR, IC ) = ZERO
        ENDDO
      ENDIF
C
C Now loop around the values which need entering
C
      DO J = 1, NNT
        IC = ( JRN( J ) - 1 )*NH + IH
        CALL MATIND (IR,IC,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
        AMAT( IROW, ICOL ) = VALS( J )
      ENDDO
C
      RETURN
      END
C*********************************************************************
