C*********************************************************************
C subroutine Double Precision Vector Right hand Vector Form **********
C            -      -         -      -          -      -    **********
C Steve Gibbons 27.6.99                                              C
C____________________________________________________________________C
C                                                                    C
C This routine merely calls all the DPV*** routines which fill in    C
C the appropriate terms in the right hand vector. You must be        C
C careful in supplying the coefficients. This routine ADDS all the   C
C terms to the vector, regardless of the expression in the Navier    C
C stokes and heat equations.                                         C
C                               (vectors denoted by capital letters) C
C If the inputs are  V and Theta   then we add to the RVEC,          C
C                                                                    C
C  V . R ( cfb1  +  cfb2 / r^3 )  to the temperature term,           C
C                                                                    C
C  cfc [ V . Grad ( Theta ) ] to the temperature term,               C
C                                                                    C
C  cfd Nabla^2 [ Theta ] to the temperature term,                    C
C                                                                    C
C  cff curl [ V . grad ( V ) ] to the momentum equation.             C
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
C     NDIM      : Length of the solution vector.                     C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of harmonics (all types)                    C
C                                                                    C
C     LH        : Level of harmonics.                                C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
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
C     RVEC      : DP vector of dimension ( NDIM )                    C
C     OLDVEC    : DP vector of dimension ( NDIM ) - old vector       C
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
C                  DPPARS(  5 ) = CFC                                C
C                  DPPARS(  6 ) = CFD                                C
C                  DPPARS(  7 ) = CFE                                C
C                  DPPARS(  8 ) = CFF                                C
C                  DPPARS(  9 ) = CFG                                C
C                  DPPARS( 10 ) = CFH                                C
C                  DPPARS( 11 ) = CFI                                C
C                  DPPARS( 12 ) = TFAC                               C
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
C     SV1       : Workspace of dimension ( NDIM )                    C
C     SV2       : Workspace of dimension ( NDIM )                    C
C     SV3       : Workspace of dimension ( NDIM )                    C
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
C     DSHC      : The derivatives of the above at the radial nodes.  C
C                  Dimension is ( LH*( LH + 2) )                     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPVRVF ( NDIM, NR, NH, LH, NTHPTS, NPHPTS, MHT, MHL,
     1                    MHM, MHC, RI, RO, RVEC, OLDVEC, GAUX, GAUW,
     2                    PA, DPA, DPPARS, ORD, SV1, SV2, SV3, QST,
     3                    VF1, VF2, VF3, FTF1, FTF2, FTF3, SF,
     4                    SHC, DSHC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM, NR, NH, LH, NTHPTS, NPHPTS, 
     1        MHT( NH ), MHL( NH ), MHM( NH ), MHC( NH )
      DOUBLE PRECISION RI, RO, RVEC( NDIM ), OLDVEC( NDIM ),
     1                 GAUX( NTHPTS ), GAUW( NTHPTS ),
     2                 PA( (LH+1)*(LH+2)/2 , NTHPTS ),
     3                 DPA( (LH+1)*(LH+2)/2 , NTHPTS )
      DOUBLE PRECISION DPPARS( * )
      CHARACTER *(2) ORD
      DOUBLE PRECISION SV1( NDIM ), SV2( NDIM ), SV3( NDIM ),
     1                 QST( LH*(LH + 2), 3 ),
     2                 VF1( NPHPTS, NTHPTS, 3),
     3                 VF2( NPHPTS, NTHPTS, 3),
     4                 VF3( NPHPTS, NTHPTS, 3)
      DOUBLE PRECISION FTF1( 2*NPHPTS ), FTF2( 2*NPHPTS ),
     1                 FTF3( 2*NPHPTS )
      DOUBLE PRECISION SF( NTHPTS, NPHPTS ),
     4                 SHC( LH*(LH+2) ), DSHC( LH*(LH+2) )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP
      DOUBLE PRECISION CVAL, CFA, CFB1, CFB2, CFC, CFD, CFE, CFF,
     1                 CFG, CFH, CFI, ZERO, X1, X2, FAC,
     2                 TFAC, FACT, FACV, DPONE, TOL
      LOGICAL POWT
      PARAMETER ( ZERO = 0.0d0, DPONE = 1.0d0, TOL = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of NDIM against NH and NR etc.
C
      IF ( NDIM.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPVRVF, bad array size'
        PRINT *,' NDIM = ',NDIM
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
      CVAL  = DPPARS(  1 ) 
      CFA   = DPPARS(  2 ) 
      CFB1  = DPPARS(  3 ) 
      CFB2  = DPPARS(  4 ) 
      CFC   = DPPARS(  5 ) 
      CFD   = DPPARS(  6 ) 
      CFE   = DPPARS(  7 ) 
      CFF   = DPPARS(  8 ) 
      CFG   = DPPARS(  9 ) 
      CFH   = DPPARS( 10 ) 
      CFI   = DPPARS( 11 ) 
      TFAC  = DPPARS( 12 ) 
C
C Check the values of TFAC and CVAL
C
      IF ( ABS( CVAL ).GT.TOL  .AND.
     1     ABS( TFAC ).GT.TOL ) THEN
        PRINT *,'Subroutine DPVRVF.'
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
        PRINT *,'Subroutine DPVRVF. Fault with nphpts or nthpts.'
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
C Zero the array RVEC
C
      IOP = 0
      CALL VECOP(  RVEC, ZERO, NDIM, IOP )
C
C Now, do the stratification term ...
C
      FAC = 1.0d0
      CALL DPVSTA ( NDIM, NR, NH, MHT, MHL, MHM, MHC, RI, RO,
     1              RVEC, OLDVEC, FAC, CFB1, CFB2 )
C
C Now, do the non-linear term V. Grad ( Theta )
C
      FAC = CFC
      CALL DPVVTA ( LH, NDIM, NTHPTS, NPHPTS, NR, NH, MHT, MHL,
     1              MHM, MHC, ORD, RVEC, OLDVEC, RI, RO, FAC, SF,
     2              QST, VF1, VF2, FTF1, FTF2, FTF3,
     3              GAUX, GAUW, PA, DPA, SHC, DSHC )
C
C Now do thermal diffusivity
C
      FAC = CFD
      CALL DPVL3A ( NDIM, NR, NH, MHT, MHL, RI, RO,
     1              RVEC, OLDVEC, FAC, ORD )
C
C Now do non-linear inertial term
C
      FAC = CFF
      CALL DPVVVA ( LH, NDIM, NTHPTS, NPHPTS, NR, NH, MHT, MHL,
     1              MHM, MHC, ORD, RVEC, OLDVEC, RI, RO, FAC,
     2              SV1, SV2, SV3, QST, VF1, VF2, VF3, FTF1,
     3              FTF2, FTF3, GAUX, GAUW, PA, DPA )
C
C Now add curl ( K \times V ) term
C
      FAC = CFG
      CALL DPVCTA ( LH, NDIM, NTHPTS, NPHPTS, NR, NH, MHT, MHL,
     1              MHM, MHC, ORD, RVEC, OLDVEC, RI, RO, FAC,
     2              QST, SV1, SV2, VF1, FTF1, FTF2, FTF3,
     3              GAUX, GAUW, PA, DPA )
C
C Now add thermal buoyancy term
C
      FAC = CFH
      CALL DPVBTA ( NDIM, NR, NH, MHT, MHL, MHM, MHC,
     1              RVEC, OLDVEC, FAC )
C
C Now add viscous diffusivities
C
      FAC = CFI
      CALL DPVL1A ( NDIM, NR, NH, MHT, MHL, RI, RO,
     1              RVEC, OLDVEC, FAC, ORD )
      CALL DPVL2A ( NDIM, NR, NH, MHT, MHL, RI, RO,
     1              RVEC, OLDVEC, FAC, ORD )
C
C Now add time derivatives
C
      CALL DPVDFT ( NDIM, NR, NH, MHT, MHL, MHM, MHC, RI, RO,
     1              RVEC, OLDVEC, CFA, CFE, ORD, CVAL )
C
C
C
      FACT = CFA*TFAC
      FACV = CFE*TFAC
      CALL DPVTDA ( NDIM, NR, NH, MHT, MHL, RI, RO, RVEC,
     1              OLDVEC, FACT, FACV, ORD, DPONE, DPONE )
C
      RETURN
      END
C*********************************************************************
