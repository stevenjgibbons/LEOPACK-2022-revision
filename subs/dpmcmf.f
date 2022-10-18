C*********************************************************************
C subroutine Double Precision Matrix Convection Matrix Form **********
C            -      -         -      -          -      -    **********
C Steve Gibbons 27.6.99                                              C
C____________________________________________________________________C
C                                                                    C
C This routine merely calls all the DPM*** routines which fill in    C
C the appropriate terms in the convection matrix. You must be        C
C careful in supplying the coefficients. This routine ADDS all the   C
C terms to the vector, regardless of the expression in the Navier    C
C stokes and heat equations.                                         C
C                               (vectors denoted by capital letters) C
C If the inputs are  V and Theta   then we add to the RVEC,          C
C                                                                    C
C  V . R ( cfb1  +  cfb2 / r^3 )  to the temperature term,           C
C                                                                    C
C  facvt0 [ V . Grad ( Theta_0 ) ] to the temperature term,          C
C  facv0t [ V_0 . Grad ( Theta ) ] to the temperature term,          C
C                                                                    C
C  cfd Nabla^2 [ Theta ] to the temperature term,                    C
C                                                                    C
C  facvv0 curl [ V . grad ( V_0 ) ] to the momentum equation.        C
C  facv0v curl [ V_0 . grad ( V ) ] to the momentum equation.        C
C                                                                    C
C  cfg curl [ K \times V ] to the momentum equation.                 C
C                                                                    C
C  cfh curl [ Theta R ] to the momentum equation.                    C
C                                                                    C
C  cfi curl [ Nabla^2 [ V ] ] to the momentum equation.              C
C                                                                    C
C  cfa dTheta / dt to the temperature term with drift rate CVAL.     C
C  cfe dV / dt to the  momentum equation with drift rate CVAL.       C
C                                                                    C
C  tfac*cfa dTheta / dt to the temperature term - no drift speed     C
C  tfac*cfe dV / dt to the  momentum equation - no drift speed       C
C                                                                    C
C  Routine checks that either CVAL or TFAC is zero - it makes no     C
C  sense if both are non-zero.                                       C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     LVEC      : Length of the solution vector.                     C
C     NR        : Number of radial grid nodes.                       C
C     NNH       : Number of harmonics in matrix vector.              C
C     NOH       : Number of harmonics in old solution vector.        C
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
C     LH        : Level of harmonics.                                C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C                                                                    C
C     MOHT       : Dimension ( NOH )                                 C
C     MOHL       : Dimension ( NOH )                                 C
C     MOHM       : Dimension ( NOH )                                 C
C     MOHC       : Dimension ( NOH )                                 C
C                                                                    C
C  moht, mohl, mohm and mohc define the spherical harmonics present  C
C  in the old solution vector. For spherical harmonic number I;      C
C                                                                    C
C   MOHT( I ) = 1 for a poloidal velocity vector                     C
C   MOHT( I ) = 2 for a toroidal velocity vector                     C
C   MOHT( I ) = 3 for a temperature / codensity term                 C
C                                                                    C
C   MOHL( I ) = spherical harmonic degree, l                         C
C                                                                    C
C   MOHM( I ) = spherical harmonic order, m                          C
C                                                                    C
C   MOHC( I ) = 1 for a cosine dependence in phi and                 C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C     MNHT       : Dimension ( NNH )                                 C
C     MNHL       : Dimension ( NNH )                                 C
C     MNHM       : Dimension ( NNH )                                 C
C     MNHC       : Dimension ( NNH )                                 C
C                                                                    C
C  mnht, mnhl, mnhm and mnhc define the spherical harmonics present  C
C  in the new solution vector. For spherical harmonic number I;      C
C                                                                    C
C   MNHT( I ) = 1 for a poloidal velocity vector                     C
C   MNHT( I ) = 2 for a toroidal velocity vector                     C
C   MNHT( I ) = 3 for a temperature / codensity term                 C
C                                                                    C
C   MNHL( I ) = spherical harmonic degree, l                         C
C                                                                    C
C   MNHM( I ) = spherical harmonic order, m                          C
C                                                                    C
C   MNHC( I ) = 1 for a cosine dependence in phi and                 C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C     IDHSF     : Different harmonic set flag.                       C
C                  = 1 if the harmonics which constitute OLDVEC      C
C                      are identical to those which constitute the   C
C                      matrix vector.                                C
C                                                                    C
C                  = 2 if the harmonics which constitute OLDVEC      C
C                      are different to those which constitute the   C
C                      matrix vector.                                C
C                                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     AMAT      : Matrix. Dimensions ( N1, N2 )                      C
C              Will generally be banded due to the nature of the     C
C              numerical scheme. KL, KU and KLE parameterise this.   C
C                                                                    C
C     OLDVEC    : DP vector of dimension ( LVEC ) - old vector       C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Gaussian weights evaluated by the routine          C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH+ 1 )*( LH+ 2 )/2 , NTHPTS }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     DPPARS    : Array containing parameters. Arranged as followe   C
C                                                                    C
C                  DPPARS(  1 ) = CVAL                               C
C                  DPPARS(  2 ) = CFA                                C
C                  DPPARS(  3 ) = CFB1                               C
C                  DPPARS(  4 ) = CFB2                               C
C                  DPPARS(  5 ) = FACVT0                             C
C                  DPPARS(  6 ) = FACV0T                             C
C                  DPPARS(  7 ) = CFD                                C
C                  DPPARS(  8 ) = CFE                                C
C                  DPPARS(  9 ) = FACVV0                             C
C                  DPPARS( 10 ) = FACV0V                             C
C                  DPPARS( 11 ) = CFG                                C
C                  DPPARS( 12 ) = CFH                                C
C                  DPPARS( 13 ) = CFI                                C
C                  DPPARS( 14 ) = TFAC                               C
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
C Working variables :-                                               C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     QST       : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector                           C
C                  Has dimensions (  LH*(LH+2) , 3).                 C
C     VF1       : Vector Function. An array of dimensions            C
C                  ( NPHPTS, NTHPTS, 3) which contain the R, THETA   C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IPHI, ITHETA, 2 ) is the Theta      C
C                  compontent of the vector at (iphi, itheta).       C
C     VF2       : Second vector Function. See above                  C
C     VF3       : Third vector Function. See above                   C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C     SF        : Scalar Function. Dimensions ( NTHPTS, NPHPTS )     C
C     SHC       : Spherical Harmonic Coefficients.                   C
C                  Has the usual ordering.                           C
C                   SHC( l*l ) = P_l^0                               C
C                   SHC( l*l + 2*m - 1 ) = P_l^mc ( m non zero )     C
C                   SHC( l*l + 2*m ) = P_l^ms ( m non zero )         C
C                  Dimension is ( LH*( LH + 2) )                     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPMCMF ( N1, N2, LVEC, NR, NNH, NOH, KL, KU, KLE,
     1           IMF, LH, NTHPTS, NPHPTS, MOHT, MOHL, MOHM, MOHC,
     2           MNHT, MNHL, MNHM, MNHC, IDHSF, RI, RO, AMAT, OLDVEC,
     3           GAUX, GAUW, PA, DPA, DPPARS,
     4           QST, VF1, VF2, VF3, FTF1, FTF2, FTF3, SF, SHC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, LVEC, NR, NNH, NOH, KL, KU, KLE, IMF, LH,
     1        NTHPTS, NPHPTS, IDHSF,
     2        MOHT( NOH ), MOHL( NOH ), MOHM( NOH ), MOHC( NOH ),
     3        MNHT( NNH ), MNHL( NNH ), MNHM( NNH ), MNHC( NNH )
      DOUBLE PRECISION RI, RO, AMAT( N1, N2 ), OLDVEC( LVEC ),
     1                 GAUX( NTHPTS ), GAUW( NTHPTS ),
     2                 PA( (LH+1)*(LH+2)/2 , NTHPTS ),
     3                 DPA( (LH+1)*(LH+2)/2 , NTHPTS )
      DOUBLE PRECISION DPPARS( * )
      CHARACTER *(2) ORD
      DOUBLE PRECISION 
     1                 QST( LH*(LH + 2), 3 ),
     2                 VF1( NPHPTS, NTHPTS, 3),
     3                 VF2( NPHPTS, NTHPTS, 3),
     4                 VF3( NPHPTS, NTHPTS, 3)
      DOUBLE PRECISION FTF1( 2*NPHPTS ), FTF2( 2*NPHPTS ),
     1                 FTF3( 2*NPHPTS )
      DOUBLE PRECISION SF( NTHPTS, NPHPTS ),
     4                 SHC( LH*(LH+2) )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP
      DOUBLE PRECISION CVAL, CFA, CFB1, CFB2, FACVT0, FACV0T, CFD,
     1                 CFE, FACVV0, FACV0V, CFG, CFH, CFI, ZERO,
     2                 X1, X2, FAC, FACT, TFAC, FACV, TOL
      LOGICAL POWT
      PARAMETER ( ZERO = 0.0d0, TOL = 1.0d-6 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of N2 against NNH and NR etc.
C
      IF ( N2.NE.NNH*NR ) THEN
        PRINT *,' Subroutine DPMCMF, bad array size'
        PRINT *,' N2   = ',N2
        PRINT *,' NNH = ',NNH,' NR = ',NR
        STOP
      ENDIF
C
C Check value of LVEC against NOH and NR etc.
C
      IF ( LVEC.NE.NOH*NR ) THEN
        PRINT *,' Subroutine DPMCMF, bad array size'
        PRINT *,' LVEC = ',LVEC
        PRINT *,' NOH = ',NOH,' NR = ',NR
        STOP
      ENDIF
C
C Check value of IMF
C
      IF ( IMF.NE.1 ) THEN
        PRINT *,' Subroutine DPMCMF, IMF = ', IMF
        PRINT *,' IMF = 1 is the only current option.'
        STOP
      ENDIF
C
C Check the value of N1 for particular values of IMF
C
      IF ( IMF.EQ.1 ) THEN
        IF ( N1.NE.(KL+KU+KLE) ) THEN
           PRINT *,' Subroutine DPMCMF, N1 = ', N1
           PRINT *,' KL  = ', KL
           PRINT *,' KU  = ', KU
           PRINT *,' KLE = ', KLE
           STOP
        ENDIF
      ENDIF
C
C Check the IDHSF flag
C
      IF ( IDHSF.NE.1 .AND. IDHSF.NE.2 ) THEN
        PRINT *,' Subroutine DPMCMF, IDHSF = ', IDHSF
        PRINT *,' IDHSF = 1 and 2 are the only options.'
        STOP
      ENDIF
C
C Check the IDHSF flag against different numbers
C of harmonics - doesn't bother checking the sets
C 
      IF ( IDHSF.EQ.1 .AND. NNH.NE.NOH ) THEN
        PRINT *,' Subroutine DPMCMF, IDHSF = ', IDHSF
        PRINT *,' NNH     = ', NNH
        PRINT *,' NOH     = ', NOH
        STOP
      ENDIF
C
C
      CVAL   = DPPARS(  1 ) 
      CFA    = DPPARS(  2 ) 
      CFB1   = DPPARS(  3 ) 
      CFB2   = DPPARS(  4 ) 
      FACVT0 = DPPARS(  5 ) 
      FACV0T = DPPARS(  6 ) 
      CFD    = DPPARS(  7 ) 
      CFE    = DPPARS(  8 ) 
      FACVV0 = DPPARS(  9 ) 
      FACV0V = DPPARS( 10 ) 
      CFG    = DPPARS( 11 ) 
      CFH    = DPPARS( 12 ) 
      CFI    = DPPARS( 13 ) 
      TFAC   = DPPARS( 14 ) 
C
C Check the values of TFAC and CVAL
C
      IF ( ABS( CVAL ).GT.TOL  .AND.
     1     ABS( TFAC ).GT.TOL ) THEN
        PRINT *,'Subroutine DPMCMF.'
        PRINT *,' CVAL  = ', CVAL
        PRINT *,' TFAC  = ', TFAC
        STOP
      ENDIF
C
C Check the values of NTHPTS and NPHPTS
C
      CALL POWTWO( NPHPTS, POWT )
      IF ( 2*LH.GT.NPHPTS .OR. 2*LH.GT.NTHPTS .OR.
     1     ( .NOT. POWT )                    ) THEN
        PRINT *,'Subroutine DPMCMF. Fault with nphpts or nthpts.'
        PRINT *,'LH = ',LH,' NTHPTS= ',NTHPTS,' NPHPTS= ',NPHPTS
        STOP
      ENDIF
C
C Calculate Gauss points and weights
C
      X1 = -1.0d0
      X2 = 1.0d0
      CALL GAUWTS ( X1, X2, GAUX, GAUW, NTHPTS, NTHPTS )
      CALL SCHNLA ( PA, DPA, GAUX, LH, LH, NTHPTS, NTHPTS)
C
C Zero the array AMAT
C
      IOP = 0
      CALL MATOP( AMAT, ZERO, N1, N2, IOP )
C
C Now, do the stratification term ...
C
      FAC = 1.0d0
      CALL DPMSTA ( N1, N2, NR, NNH, KL, KU, KLE, IMF, MNHT,
     1              MNHL, MNHM, MNHC, RI, RO, AMAT, FAC, CFB1, CFB2 )
C
C Now, do the non-linear term V. Grad ( Theta )
C
      IF ( IDHSF.EQ.1 ) THEN
C
         CALL DPMVGT ( N1, N2, NR, NNH, KL, KU, KLE, IMF, MNHT, MNHL,
     1                 MNHM, MNHC, LH, NPHPTS, NTHPTS, RI, RO, AMAT,
     2                 FACV0T, FACVT0, GAUX, GAUW, PA, DPA,
     3                 OLDVEC, ORD, FTF1, VF1, VF2, SF, SHC )
C
      ENDIF
C
      IF ( IDHSF.EQ.2 ) THEN
C
         CALL DPMV0T ( N1, N2, LVEC, NR, NOH, NNH, KL, KU, KLE,
     1                 IMF, MOHT, MOHL, MOHM, MOHC, MNHT, MNHL,
     1                 MNHM, MNHC, LH, NPHPTS, NTHPTS, RI, RO,
     2                 AMAT, FACV0T, GAUX, GAUW, PA, DPA,
     3                 OLDVEC, ORD, FTF1, VF1, VF2, SF, SHC )
C
         CALL DPMVT0 ( N1, N2, LVEC, NR, NOH, NNH, KL, KU, KLE,
     1                 IMF, MOHT, MOHL, MOHM, MOHC, MNHT, MNHL,
     1                 MNHM, MNHC, LH, NPHPTS, NTHPTS, RI, RO,
     2                 AMAT, FACVT0, GAUX, GAUW, PA, DPA,
     3                 OLDVEC, ORD, FTF1, VF1, VF2, SF, SHC )
C
      ENDIF
C
C Now do thermal diffusivity
C
      FAC = CFD
      CALL DPML3A ( N1, N2, NR, NNH, KL, KU, KLE, IMF, MNHT,
     1              MNHL, RI, RO, AMAT, FAC, ORD )
C
C Now do non-linear inertial term
C
      IF ( IDHSF.EQ.1 ) THEN
C
         CALL DPMVV ( N1, N2, NR, NNH, KL, KU, KLE, IMF, MNHT, MNHL,
     1                MNHM, MNHC, LH, NPHPTS, NTHPTS, RI, RO, AMAT,
     2                FACV0V, FACVV0, GAUX, GAUW, PA, DPA,
     3                OLDVEC, ORD, FTF1, FTF2, FTF3, VF1, VF2,
     4                VF3, QST )
C
      ENDIF
C
      IF ( IDHSF.EQ.2 ) THEN
C
         CALL DPMV0V ( N1, N2, LVEC, NR, NOH, NNH, KL, KU, KLE,
     1                 IMF, MOHT, MOHL, MOHM, MOHC, MNHT, MNHL,
     2                 MNHM, MNHC, LH, NPHPTS, NTHPTS, RI, RO, AMAT,
     3                 FACV0V, GAUX, GAUW, PA, DPA,
     4                 OLDVEC, ORD, FTF1, FTF2, FTF3, VF1, VF2,
     5                 VF3, QST )
C
         CALL DPMVV0 ( N1, N2, LVEC, NR, NOH, NNH, KL, KU, KLE,
     1                 IMF, MOHT, MOHL, MOHM, MOHC, MNHT, MNHL,
     2                 MNHM, MNHC, LH, NPHPTS, NTHPTS, RI, RO, AMAT,
     3                 FACVV0, GAUX, GAUW, PA, DPA,
     4                 OLDVEC, ORD, FTF1, FTF2, FTF3, VF1, VF2,
     5                 VF3, QST )
C
      ENDIF
C
C Now add curl ( K \times V ) term
C
      FAC = CFG
      CALL DPMCTA ( N1, N2, NR, NNH, KL, KU, KLE, IMF, MNHT,
     1              MNHL, MNHM, MNHC, RI, RO, AMAT, FAC,
     2              LH, NTHPTS, NPHPTS, GAUX, GAUW, QST,
     3              FTF1, FTF2, FTF3, VF1, PA, DPA, ORD )
C
C Now add thermal buoyancy term
C
      FAC = CFH
      CALL DPMBTA ( N1, N2, NR, NNH, KL, KU, KLE, IMF, MNHT,
     1              MNHL, MNHM, MNHC, RI, RO, AMAT, FAC )
C
C Now add viscous diffusivities
C
      FAC = CFI
      CALL DPML1A ( N1, N2, NR, NNH, KL, KU, KLE, IMF, MNHT,
     1              MNHL, RI, RO, AMAT, FAC, ORD )
      CALL DPML2A ( N1, N2, NR, NNH, KL, KU, KLE, IMF, MNHT,
     1              MNHL, RI, RO, AMAT, FAC, ORD )
C
C Now add time derivatives
C
      CALL DPMDFT ( N1, N2, NR, NNH, KL, KU, KLE, IMF, MNHT,
     1              MNHL, MNHM, MNHC, RI, RO, AMAT, CFA, CFE,
     2              ORD, CVAL )
C
      FACT = CFA*TFAC
      FACV = CFE*TFAC
      CALL DPMTDA ( N1, N2, NR, NNH, KL, KU, KLE, IMF, MNHT,
     1              MNHL, RI, RO, AMAT, FACT, FACV, ORD )
C
      RETURN
      END
C*********************************************************************
