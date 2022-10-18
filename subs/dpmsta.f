C*********************************************************************
C subroutine Double Precision Matrix Stratification Term Add *********
C            -      -         -      -              -    -   *********
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C Adds the term to the double precision matrix resulting from        C
C the V . grad ( t_0 ( r ) ) term in the heat equation.              C
C                                                                    C
C The general form is t_0 ( r ) = CB1 * r + CB2 / r^2 and so in      C
C effect, we are adding  ( CB1 + CB2 / r^3 ) to the poloidal         C
C velocity term.                                                     C
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
C                                                                    C
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C     MHM       : Dimension ( NH )                                   C
C     MHC       : Dimension ( NH )                                   C
C                                                                    C
C  mht, mhl, mhm and mhc define the spherical harmonics present      C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C   MHT( I ) = 4 for a poloidal magnetic field vector                C
C   MHT( I ) = 5 for a poloidal magnetic field vector                C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     AMAT      : Matrix. Dimensions ( N1, N2 )                      C
C              Will generally be banded due to the nature of the     C
C              numerical scheme. KL, KU and KLE parameterise this.   C
C                                                                    C
C     FAC       : Multiplication factor for term                     C
C                                                                    C
C     CB1       : See above                                          C
C     CB2       : See above                                          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPMSTA ( N1, N2, NR, NH, KL, KU, KLE, IMF, MHT,
     1                MHL, MHM, MHC, RI, RO, AMAT, FAC, CB1, CB2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, NR, NH, KL, KU, KLE, IMF,
     1        MHT( NH ), MHM( NH ), MHC( NH ), MHL( NH )
      DOUBLE PRECISION RI, RO, AMAT( N1, N2 ), FAC, CB1, CB2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      CHARACTER *(2) ORD
      EXTERNAL STRTS
      INTEGER IIH, IOH, LI, MI, ICSI, INTPAR( 1 )
      DOUBLE PRECISION DPPARS( 2 ), VEC( 1 ), TOL
      PARAMETER ( TOL = 1.0d-8 )
C nb - vec is not referred to by STRTS
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of N2 against NH and NR etc.
C
      IF ( N2.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPMSTA, bad array size'
        PRINT *,' N2 (second array dimension) = ',N2
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FAC ).LT.TOL ) RETURN
C
      DPPARS( 1 ) = CB1
      DPPARS( 2 ) = CB2
C
C Since the stratification term does not involve any
C derivatives of the velocity, ORD is irrelevant and
C may be arbitrarily set.
C
      ORD = 'SS'
C
C Loop around the poloidal velocity harmonics and seek
C the equivalent temperature harmonics
C
      DO IIH = 1, NH
        IF ( MHT( IIH ).EQ.1 ) THEN
c
c       . we have a poloidal velocity harmonic
c       . so loop around to find equivalent
c       . temperature harmonic.
c
         LI   = MHL( IIH )
         MI   = MHM( IIH )
         ICSI = MHC( IIH )
         DO IOH = 1, NH
           IF (    MHT( IOH ).EQ.3       .AND.
     1             MHL( IOH ).EQ.LI      .AND.
     2             MHM( IOH ).EQ.MI      .AND.
     3             MHC( IOH ).EQ.ICSI    ) THEN
              INTPAR( 1 ) = LI
              CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                   IOH, INTPAR, STRTS, ORD, RI, RO, AMAT,
     2                   DPPARS, VEC, FAC )
              GOTO 500
           ENDIF
         ENDDO
c
        ENDIF
 500    CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************

