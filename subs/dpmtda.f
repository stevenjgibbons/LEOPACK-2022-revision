C*********************************************************************
C subroutine Double Precision Matrix Time Derivative Add *************
C            -      -         -      -    -          -   *************
C Steve Gibbons 22.3.99                                              C
C____________________________________________________________________C
C Adds the term to the double precision matrix resulting from        C
C the time derivatives for both velocity and temperature.            C
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
C                                                                    C
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C                                                                    C
C  mht and mhl define the spherical harmonics present                C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     AMAT      : Matrix. Dimensions ( N1, N2 )                      C
C              Will generally be banded due to the nature of the     C
C              numerical scheme. KL, KU and KLE parameterise this.   C
C                                                                    C
C     FACT      : Multiplication factor for temperature term         C
C     FACV      : Multiplication factor for velocity term            C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     ORD       : *(2). Determines the order of accuracy. Options :  C
C             SS - Strictly second order                             C
C             SF - Strictly fourth order                             C
C             O5 - Optimum accuracy for bandwidth 5; this gives      C
C                  Fourth order accuracy for 1st and 2nd derivatives C
C                  and second order accuracy for 3rd and 4th der.s   C
C             O7 - Optimum accuracy for bandwidth 7; this gives      C
C                  Sixth order accuracy for 1st and 2nd derivatives  C
C                  and fourth order accuracy for 3rd and 4th der.s   C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPMTDA ( N1, N2, NR, NH, KL, KU, KLE, IMF, MHT,
     1                    MHL, RI, RO, AMAT, FACT, FACV, ORD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, NR, NH, KL, KU, KLE, IMF,
     1        MHT( NH ), MHL( NH )
      DOUBLE PRECISION RI, RO, AMAT( N1, N2 ), FACT, FACV
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      EXTERNAL LAPTS
      INTEGER IIH, INTPAR( 2 )
      DOUBLE PRECISION DPPARS( 1 ), VEC( 1 ), FAC, TOL
      PARAMETER ( TOL = 1.0d-6 )
C nb - neither dppars or vec are referred to by LAPTS
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of N2 against NH and NR etc.
C
      IF ( N2.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPMTDA, bad array size'
        PRINT *,' N2 (second array dimension) = ',N2
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
      DO IIH = 1, NH
        IF ( ABS(FACV).LT.TOL ) GOTO 500
        IF ( MHT( IIH ).EQ.1 ) THEN
c
         INTPAR( 1 ) = MHL( IIH )
         INTPAR( 2 ) = 1
         FAC = (-1.0d0)*FACV
         CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                IIH, INTPAR, LAPTS, ORD, RI, RO, AMAT,
     2                DPPARS, VEC, FAC )
c
        ENDIF
c
        IF ( MHT( IIH ).EQ.2 ) THEN
c
         INTPAR( 2 ) = 0
         FAC = FACV
         CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                IIH, INTPAR, LAPTS, ORD, RI, RO, AMAT,
     2                DPPARS, VEC, FAC )
c
        ENDIF
c
 500    CONTINUE
        IF ( ABS(FACT).LT.TOL ) GOTO 501
        IF ( MHT( IIH ).EQ.3 ) THEN
c
         INTPAR( 2 ) = 0
         FAC = FACT
         CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                IIH, INTPAR, LAPTS, ORD, RI, RO, AMAT,
     2                DPPARS, VEC, FAC )
c
        ENDIF
 501    CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************

