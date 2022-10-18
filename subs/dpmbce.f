C*********************************************************************
C subroutine Double Precision Matrix Boundary Condition Enforce ******
C            -      -         -      -        -         -       ******
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
C     MHM       : Dimension ( NH )                                   C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     AMAT      : Matrix. Dimensions ( N1, N2 )                      C
C              Will generally be banded due to the nature of the     C
C              numerical scheme. KL, KU and KLE parameterise this.   C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     BCIVEL    : Inner boundary velocity condition                  C
C     BCOVEL    : Outer boundary velocity condition                  C
C      bcivel and bcovel can be set to 'SF' stress free and          C
C      'NS' no slip.                                                 C
C                                                                    C
C     BCITHE    : Inner boundary temperature condition               C
C     BCOTHE    : Outer boundary temperature condition               C
C      bcithe and bcothe can be set to 'HF' heat flux and            C
C      'TM' - temperature.                                           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPMBCE ( N1, N2, NR, NH, KL, KU, KLE, IMF, MHT,
     1                    MHL, MHM, RI, RO, AMAT,
     2                    BCIVEL, BCOVEL, BCITHE, BCOTHE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, NR, NH, KL, KU, KLE, IMF,
     1        MHT( NH ), MHM( NH ), MHL( NH )
      DOUBLE PRECISION RI, RO, AMAT( N1, N2 )
      CHARACTER *(2) BCIVEL, BCOVEL, BCITHE, BCOTHE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRN, JRN( 5 ), NNT, IH, IPTT
      DOUBLE PRECISION H, VALS( 5 ), FAC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the size of the matrix
C
      IF ( N2.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPMBCE, bad array size'
        PRINT *,' N2 (second array dimension) = ',N2
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C check character string variables ...
C
      IF ( BCITHE.NE.'TM' .AND. BCITHE.NE.'HF' ) THEN
        PRINT *,' BCITHE = ', BCITHE
        STOP
      ENDIF
      IF ( BCOTHE.NE.'TM' .AND. BCOTHE.NE.'HF' ) THEN
        PRINT *,' BCOTHE = ', BCOTHE
        STOP
      ENDIF
      IF ( BCIVEL.NE.'NS' .AND. BCIVEL.NE.'SF' ) THEN
        PRINT *,' BCIVEL = ', BCIVEL
        STOP
      ENDIF
      IF ( BCOVEL.NE.'NS' .AND. BCOVEL.NE.'SF' ) THEN
        PRINT *,' BCOVEL = ', BCOVEL
        STOP
      ENDIF
C
      H = ( RO - RI )/DBLE( NR - 1 )
C
C Begin loop around all the harmonics
C
      DO IH = 1, NH
        IPTT = MHT( IH )
c       .
c       . Consider case of poloidal harmonics
c       .
        IF ( IPTT.EQ.1 ) THEN
c          .
c          . First enforce the impenetrable condition
c          . INNER BOUNDARY
           IRN = 1
           NNT = 1
           JRN( 1 ) = 1
           VALS( 1 ) = 1.0d0
           CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
c          .
c          . OUTER BOUNDARY
           IRN = NR
           NNT = 1
           JRN( 1 ) = NR
           VALS( 1 ) = 1.0d0
           CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
c          .
c          . Now enforce no-slip at inner boundary
c          .
           IF ( BCIVEL.EQ.'NS' ) THEN
             IRN = 2
             NNT = 5
             JRN( 1 ) = 1
             JRN( 2 ) = 2
             JRN( 3 ) = 3
             JRN( 4 ) = 4
             JRN( 5 ) = 5
             VALS( 1 ) = -25.0d0
             VALS( 2 ) = 48.0d0
             VALS( 3 ) = -36.0d0
             VALS( 4 ) = 16.0d0
             VALS( 5 ) = -3.0d0
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
c          . Now enforce stress-free at inner boundary
c          .
           IF ( BCIVEL.EQ.'SF' ) THEN
             IRN = 2
             NNT = 5
             JRN( 1 ) = 1
             JRN( 2 ) = 2
             JRN( 3 ) = 3
             JRN( 4 ) = 4
             JRN( 5 ) = 5
             VALS( 1 ) = 35.0d0
             VALS( 2 ) = -104.0d0
             VALS( 3 ) = 114.0d0
             VALS( 4 ) = -56.0d0
             VALS( 5 ) = 11.0d0
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
c          . Now enforce no-slip at outer boundary
c          .
           IF ( BCOVEL.EQ.'NS' ) THEN
             IRN = NR - 1
             NNT = 5
             JRN( 1 ) = NR - 4
             JRN( 2 ) = NR - 3
             JRN( 3 ) = NR - 2
             JRN( 4 ) = NR - 1
             JRN( 5 ) = NR
             VALS( 1 ) = 3.0d0
             VALS( 2 ) = -16.0d0
             VALS( 3 ) = 36.0d0
             VALS( 4 ) = -48.0d0
             VALS( 5 ) = 25.0d0
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
c          . Now enforce stress-free at outer boundary
c          .
           IF ( BCOVEL.EQ.'SF' ) THEN
             IRN = NR - 1
             NNT = 5
             JRN( 1 ) = NR - 4
             JRN( 2 ) = NR - 3
             JRN( 3 ) = NR - 2
             JRN( 4 ) = NR - 1
             JRN( 5 ) = NR
             VALS( 1 ) = 11.0d0
             VALS( 2 ) = -56.0d0
             VALS( 3 ) = 114.0d0
             VALS( 4 ) = -104.0d0
             VALS( 5 ) = 35.0d0
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
        ENDIF
c       .
c       . Consider case of toroidal harmonics
c       .
        IF ( IPTT.EQ.2 ) THEN
c          .
c          . Now enforce no-slip at inner boundary
c          .
           IF ( BCIVEL.EQ.'NS' ) THEN
             IRN = 1
             NNT = 1
             JRN( 1 ) = 1
             VALS( 1 ) = 1.0d0
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
c          . Now enforce stress-free at inner boundary
c          .
           IF ( BCIVEL.EQ.'SF' ) THEN
             IRN = 1
             NNT = 5
             JRN( 1 ) = 1
             JRN( 2 ) = 2
             JRN( 3 ) = 3
             JRN( 4 ) = 4
             JRN( 5 ) = 5
             FAC = RI/(25.0d0*RI + 12.0d0*H)
             VALS( 1 ) = -1.0d0
             VALS( 2 ) = 48.0d0*FAC
             VALS( 3 ) = -36.0d0*FAC
             VALS( 4 ) = 16.0d0*FAC
             VALS( 5 ) = -3.0d0*FAC
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
c          . Now enforce no-slip at outer boundary
c          .
           IF ( BCOVEL.EQ.'NS' ) THEN
             IRN = NR
             NNT = 1
             JRN( 1 ) = NR
             VALS( 1 ) = 1.0d0
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
c          . Now enforce stress-free at outer boundary
c          .
           IF ( BCOVEL.EQ.'SF' ) THEN
             IRN = NR
             NNT = 5
             JRN( 1 ) = NR - 4
             JRN( 2 ) = NR - 3
             JRN( 3 ) = NR - 2
             JRN( 4 ) = NR - 1
             JRN( 5 ) = NR
             FAC = RO/(12.0d0*H - 25.0d0*RO)
             VALS( 1 ) = 3.0d0*FAC
             VALS( 2 ) = -16.0d0*FAC
             VALS( 3 ) = 36.0d0*FAC
             VALS( 4 ) = -48.0d0*FAC
             VALS( 5 ) = -1.0d0
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
        ENDIF
c       .
c       . Consider case of temperature harmonics
c       .
        IF ( IPTT.EQ.3 ) THEN
c          .
c          . Now enforce fixed temperature at inner boundary
c          .
           IF ( BCITHE.EQ.'TM' ) THEN
             IRN = 1
             NNT = 1
             JRN( 1 ) = 1
             VALS( 1 ) = 1.0d0
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
c          . Now enforce fixed heat flux at inner boundary
c          .
           IF ( BCITHE.EQ.'HF' ) THEN
             IRN = 1
             NNT = 5
             JRN( 1 ) = 1
             JRN( 2 ) = 2
             JRN( 3 ) = 3
             JRN( 4 ) = 4
             JRN( 5 ) = 5
             FAC = 1.0d0/(12.0d0*H)
             VALS( 1 ) = -25.0d0*FAC
             VALS( 2 ) = 48.0d0*FAC
             VALS( 3 ) = -36.0d0*FAC
             VALS( 4 ) = 16.0d0*FAC
             VALS( 5 ) = -3.0d0*FAC
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
c          . Now enforce fixed temperature at outer boundary
c          .
           IF ( BCOTHE.EQ.'TM' ) THEN
             IRN = NR
             NNT = 1
             JRN( 1 ) = NR
             VALS( 1 ) = 1.0d0
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .
c          . Now enforce fixed heat flux at outer boundary
c          .
           IF ( BCOTHE.EQ.'HF' ) THEN
             IRN = NR
             NNT = 5
             JRN( 1 ) = NR
             JRN( 2 ) = NR - 1
             JRN( 3 ) = NR - 2
             JRN( 4 ) = NR - 3
             JRN( 5 ) = NR - 4
             FAC = 1.0d0/(12.0d0*H)
             VALS( 1 ) = 25.0d0*FAC
             VALS( 2 ) = -48.0d0*FAC
             VALS( 3 ) = 36.0d0*FAC
             VALS( 4 ) = -16.0d0*FAC
             VALS( 5 ) = 3.0d0*FAC
             CALL DPMBVF ( N1, N2, NR, NH, KL, KU, KLE, IMF, NNT,
     1                    IRN, IH, JRN, AMAT, VALS )
           ENDIF
c          .

        ENDIF
c       .
      ENDDO
C
      RETURN
      END
C*********************************************************************

