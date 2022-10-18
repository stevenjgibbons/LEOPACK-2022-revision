C*********************************************************************
C subroutine MATrix INDices ******************************************
C            ---    ---     ******************************************
C Steve Gibbons                                                      C
C____________________________________________________________________C
C                                                                    C
C  Tues. Nov 3rd. 1998                                               C
C  -------------------                                               C
C In order that subroutines which write matrices can be adapted to   C
C write in different storage formats, this routine takes IR and IC   C
C (the ACTUAL coordinates of the 'physical' matrix) and calculates   C
C IROW and ICOL (the row and column number in which the number is    C
C to be placed). Hence, the same code will be able to write to a     C
C matrix of any format.                                              C
C  The code also checks that the bounds of the matrix have not       C
C been exceeded.                                                     C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Row from theoretical matrix.                       C
C     IC        : Column from theoretical matrix.                    C
C     IMF       : Matrix format flag.                                C
C                                                                    C
C          imf = 1; Matrix is in LAPACK banded format                C
C                   ie element a_{i,j} is stored in                  C
C                   A( kle + ku + 1 + i - j , j )                    C
C                                                                    C
C          imf = 2; Matrix is banded but with element a_{i,j}        C
C                   stored in A( kl + 1 + j - i , i ). The matrices  C
C                   in Dave Gubbins' dynamo codes are in this        C
C                   format.                                          C
C                                                                    C
C          imf = 3; Matrix is square - ie a_{i,j} is stored          C
C                   in A( i, j ).                                    C
C                                                                    C
C     KL        : Number of lower diagonals in band matrix.          C
C     KU        : Number of upper diagonals in band matrix.          C
C     KLE       : Number of additional lower diagonals               C
C                  in band matrix. This should be set to zero for    C
C                  most applications but should be equal to KL if    C
C                  the matrix is to be LU decomposed by a LAPACK or  C
C                  NAG subroutine.                                   C
C                                                                    C
C  ( Note - KL, KU and KLE are not referenced if IMF is set to 3.)   C
C                                                                    C
C     N1        : Leading dimension of matrix A.                     C
C     N2        : Second dimension of matrix A.                      C
C                                                                    C
C       ( Note - i.e. A is dimensioned A( N1, N2 ) ... )             C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Integer                                                           C
C  -------                                                           C
C     IROW      : 'i' coordinate in the final matrix.                C
C     ICOL      : 'j' coordinate in the final matrix.                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MATIND (IR,IC,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR,IC,IMF,KL,KU,KLE,N1,N2,IROW,ICOL
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NJ
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C  First check for errors on the definition of IMF
C
      IF ( IMF.NE.1 .AND. IMF.NE.2 .AND. IMF.NE.3 ) THEN
         PRINT *,' Subroutine MATIND '
         PRINT *,' IMF = ', IMF
         PRINT *,' Illegal value. Program aborted.'
         STOP
      ENDIF
C
C  Check that N1 is equal to (KLE + KL + KU + 1) for
C  the banded matrices.
C
      IF ( IMF.EQ.1 .OR. IMF.EQ.2 ) THEN
         NJ = KLE + KL + KU + 1
         IF ( N1.NE.NJ ) THEN
           PRINT *,' Subroutine MATIND '
           PRINT *,' KL  = ', KL 
           PRINT *,' KU  = ', KU 
           PRINT *,' KLE = ', KLE
           PRINT *,' N1  = ', N1 
           PRINT *,' Illegal value. Program aborted.'
           STOP
         ENDIF
      ENDIF
C
C  Do case IMF = 1. So this is a LAPACK format band matrix.
C
      IF ( IMF.EQ.1 ) THEN
         ICOL = IC
         IROW = KLE + KU + 1 + IR - IC
      ENDIF
C
C  Do case IMF = 2. So this is a DG format band matrix.
C
      IF ( IMF.EQ.2 ) THEN
         ICOL = IR
         IROW = KL + 1 + IC - IR
      ENDIF
C
C  Do case IMF = 3. So this is a square matrix.
C
      IF ( IMF.EQ.3 ) THEN
         ICOL = IC
         IROW = IR
      ENDIF
C
C  Final error check on IROW and ICOL
C
      IF ( IROW.LT.1 .OR. IROW.GT.N1 ) THEN
         PRINT *,' Subroutine MATIND '
         PRINT *,' IMF  = ',IMF
         PRINT *,' IROW = ', IROW
         PRINT *,' 1st. dimension = ', N1
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( ICOL.LT.1 .OR. ICOL.GT.N2 ) THEN
         PRINT *,' Subroutine MATIND '
         PRINT *,' IMF  = ',IMF
         PRINT *,' ICOL = ', ICOL
         PRINT *,' 2nd. dimension = ', N2
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C
      RETURN
      END
C*********************************************************************
