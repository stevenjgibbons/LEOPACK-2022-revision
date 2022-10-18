C*********************************************************************
C subroutine Kinematic Dynamo Critical Magnetic Reynolds Number ******
C            -         -      -        -        -        -      ******
C Steve Gibbons Wed Apr  4 11:52:24 MET DST 2001                     C
C____________________________________________________________________C
C                                                                    C
C Iterates over RM for a fixed velocity to find a critical magnetic  C
C Reynolds number. Steady or oscillatory solutions.                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     MAXNVI   : Maximum number of vector interactions.              C
C                This is an upper limit for NVI and the dimension    C
C                of the arrays IHNALP, IHNBET, IHNGAM, CVI and       C
C                TVHI.                                               C
C                                                                    C
C     IHNALP   : Number of alpha harmonics. Dim ( * )                C
C     IHNBET   : Number of beta harmonics. Dim ( * )                 C
C     IHNGAM   : Number of gamma harmonics. Dim ( * )                C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARR     : Int. parameter array corresponding to VI.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR      See INDFUN for details       C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     IN0       : Int. parameter array corresponding to VEC0.        C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 IN0( 1 ) = IFORMF0                                 C
C                 IN0( 2 ) = NRR      See INDFUN for details         C
C                 IN0( 3 ) = NH0                                     C
C                                                                    C
C     MHT      : Array length ( * ) - atleast length NH              C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C  MHT( ih ) must only equal 4 and 5 for this routine ...            C
C                                                                    C
C     MHL      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. degree, l.                             C
C     MHM      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHP      : Array length ( * ) - atleast length NH              C
C                  Pointer array to finite difference coefficients.  C
C                                                                    C
C     MT0      : Array length ( * ) - atleast length NH0             C
C                  See above for key. (corresponds to VEC0).         C
C     ML0      : Array length ( * ) - atleast length NH0             C
C                  Sph. harm. degree, l.                             C
C     MM0      : Array length ( * ) - atleast length NH0             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MP0      : Array length ( * ) - atleast length NH0             C
C                  Pointer array to finite difference coefficients.  C
C                  MP0( ih ) = is, which is the 4th index of         C
C                  array SVFDC0 - indicates f.d. scheme used.        C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences for mag. field.               C
C                                                                    C
C     NBN0      : Number of nodes on each side of point for          C
C                  central differences for velocity.                 C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NTHP      : The number of theta points.                        C
C     NPHP      : The number of phi points.                          C
C                                                                    C
C     KL        : Number of lower diagonals in matrix                C
C                                                                    C
C     NCFM      : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     NCFM0     : Leading dim of SVFDC0. See SVFDCF.                 C
C                  (Must be atleast 2*NBN0 + 1 )                     C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C      (NDRVM must be atleast 2 and NDRVS must be 2 for atleast      C
C       grid nodes IR = 2, NR - 1 ).                                 C
C                                                                    C
C     N1        : First dimension of matrix, A.                      C
C     N2        : Second dimension of matrix, A.                     C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     NDCS0     : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC0.                       C
C                                                                    C
C     NEV       : Number of requested eigenvalues.                   C
C                                                                    C
C     NCV       : Length of Arnoldi factorisation.                   C
C     NCVM      : Maximum length of Arnoldi factorisation.           C
C                                                                    C
C ncv must be less than or equal to N2 and greater than or equal     C
C to (nev+2).                                                        C
C                                                                    C
C     MXIT      : Maximum number of Arnoldi iterations allowed.      C
C                                                                    C
C     IPIV      : Dim ( N2 ). Array for pivotting.                   C
C                                                                    C
C     MHIBC     : Inner boundary condition. Dim ( NDCS )             C
C                                                                    C
C  MHIBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( ih ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                                                                    C
C     MHOBC     : Outer boundary condition. Dim ( NDCS )             C
C                                                                    C
C  MHOBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( ih ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                                                                    C
C     NOITM     : Maximum number of iterations allowed.              C
C     LULOG     : Logical unit number of output file.                C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C     SVFDC0    : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM0, NR, NDRVM0+1, NDCS0 ).        C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C     CVI       : Coefficients of vector interaction. Dim ( * )      C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     FTF1      : Work array - dim. (2*NPHP)                         C
C     FTF2      : Work array - dim. (2*NPHP)                         C
C     FTF3      : Work array - dim. (2*NPHP)                         C
C     VF1       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF2       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF3       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     QST       : Work array - dim. ( LH*(LH+2) , 3)                 C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 )                             C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C     DIAGEL    : Arbitrary constant to be put in the diagonal       C
C                 element of redundant nodes. Must be non-zero to    C
C                 prevent singularity of matrix and must not be      C
C                 close to a part of the eigenspectra which is of    C
C                 interest. (For growth rate problems, DIAGEL will   C
C                 ideally be very negative.)                         C
C                                                                    C
C     DRSV      : Real part of the applied shift.                    C
C                                                                    C
C     DR        : Returned containing real part of eigenvalues.      C
C     DI        : Returned containing imag. part of eigenvalues.     C
C     D3        : Returned containing direct residuals.              C
C                                                                    C
C     WORKEV    : Dimension ( 3*NCVM ). Work array.                  C
C     WORKD     : Dimension ( 3*N2 ). Work array.                    C
C     WORKL     : Dimension ( 3*NCVM*NCVM + 6*NCVM ). Work array.    C
C                                                                    C
C     ARTOL     : Stopping criterion for Arnoldi iteration.          C
C                 (Make this small).                                 C
C                                                                    C
C                                                                    C
C     V         : Dimension ( N2, NCVM ). Returned eigenvectors.     C
C                                                                    C
C     VECRE     : Dimension ( N2 ). Returned real part of eigenvec.  C
C     VECIM     : Dimension ( N2 ). Returned imag part of eigenvec.  C
C                                                                    C
C     VEC0      : Dimension ( * ). Vector containing velocity.       C
C                                                                    C
C     DRM1      : Lower bound for Rmc                                C
C     DRM2      : Upper bound for Rmc                                C
C                                                                    C
C  A critical Rm value must lie within these bounds!                 C
C                                                                    C
C     DRMC      : Returned critical magnetic Reynolds number.        C
C     DOSC      : Returned critical oscillation frequency.           C
C     DTOL      : Stopping criterion for real part of growth rate.   C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     SELECT    : Dimension ( NCVM ). Work array.                    C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. ( MAXNVI )   C
C     CHST      : *(1) Solution type. Returned                       C
C                                                                    C
C                 'S' means stationary and converged                 C
C                 'O' means oscillatory and converged                C
C                 's' means stationary but not converged to req. tol C
C                 'o' means oscillatory + not converged to req. tol  C
C                 'N' means there didn't appear to be a root.        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE KDCMRN( MAXNVI, IHNALP, IHNBET, IHNGAM, NR, INARR,
     1       IN0, MHT, MHL, MHM, MHP, MT0, ML0, MM0, MP0, NBN, NBN0,
     2       LH, NTHP, NPHP, KL, NCFM, NCFM0, NDRVM, SVFDC, SVFDC0,
     3       A, N1, N2, NDCS, NDCS0, XARR, DIAGEL, NEV, NCV, NCVM,
     4       MXIT, DRSV, SELECT, DR, DI, D3, WORKEV, WORKD, WORKL,
     5       IPIV, ARTOL, VECRE, VECIM, V, TVHI, CVI, GAUX, GAUW,
     6       PA, DPA, FTF1, FTF2, FTF3, VF1, VF2, VF3, QST, MHIBC,
     7       MHOBC, VEC0, CHST, DRM1, DRM2, DRMC, DOSC, DTOL,
     8       NOITM, LULOG )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER MAXNVI, IHNALP( MAXNVI ), IHNBET( MAXNVI ),
     1        IHNGAM( MAXNVI ), NR, INARR( 3 ), IN0( 3 ), MHT( * ),
     2        MHL( * ), MHM( * ), MHP( * ), MT0( * ), ML0( * ),
     3        MM0( * ), MP0( * ), NBN, NBN0, LH, NTHP, NPHP, KL,
     4        NCFM, NCFM0, NDRVM, N1, N2, NDCS, NEV, NCV, NCVM, MXIT,
     5        IPIV( * ), MHIBC( NDCS ), MHOBC( NDCS ), NDCS0,
     6        NOITM, LULOG
      DOUBLE PRECISION CVI( MAXNVI ), GAUX( NTHP ), GAUW( NTHP ),
     1        DRM1, DRM2, DRMC, DOSC, DTOL
      DOUBLE PRECISION PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     1                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     2                 FTF1( 2*NPHP ), FTF2( 2*NPHP ), FTF3( 2*NPHP )
      DOUBLE PRECISION VF1( NPHP, NTHP, 3), VF2( NPHP, NTHP, 3),
     1                 VF3( NPHP, NTHP, 3 ), QST( LH*( LH + 2), 3 )
      DOUBLE PRECISION A( N1, N2 ), XARR( NR ), DIAGEL,
     1                 DR( NCVM ), DI( NCVM ), D3( NCVM ), DRSV,
     2                 WORKEV( 3*NCVM ), WORKD( 3*N2 ),
     3                 WORKL( 3*NCVM*NCVM + 6*NCVM ), ARTOL
      DOUBLE PRECISION VECRE( N2 ), VECIM( N2 ), V( N2, NCVM ),
     1                 SVFDC( NCFM, NR, NDRVM+1, NDCS ),
     2                 SVFDC0( NCFM0, NR, NDRVM+1, NDCS0 ),
     3                 VEC0( * )
      CHARACTER *(4) TVHI( MAXNVI )
      CHARACTER *(1) CHST
      LOGICAL SELECT( NCVM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NH, NH0, I, NVI, NOIT, NCE, NEGRP, NEGIP, NEG,
     1        ILN, NDRVS, IRN, IMF
      DOUBLE PRECISION DZERO, DLOW, FAC, CMTD, RM, FF1, FF2,
     1        GRMX, X, DC, DM, FH, FL, VAL1, XH, XL
      PARAMETER ( DZERO = 0.0d0, DLOW = 1.0d-7, IMF = 1 )
      LOGICAL   FIRST
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NH   = INARR( 3 )
      NH0  = IN0( 3 )
C
      ILN = 2
      IRN = NR - 1
C
      IF ( DTOL.LT.DZERO ) THEN
        PRINT *,' Subroutine KDCMRN'
        PRINT *,' DTOL = ', DTOL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NVI = 0
      CALL VCPCC( NVI, MAXNVI, IHNALP, IHNBET, IHNGAM, NH0,
     1    MT0, ML0, MM0, NH, MHT, MHL, MHM, NH, MHT, MHL, MHM,
     2    LH, NTHP, NPHP, LH, TVHI, CVI, GAUX, GAUW, PA, DPA,
     3    FTF1, FTF2, FTF3, VF1, VF2, VF3, QST )
C
      WRITE ( LULOG, * ) 'Number of vector interactions = ', NVI
C
C Begin iterations to find critical Rayleigh number
C
      NOIT = 0
 50   CONTINUE
      NOIT = NOIT + 1
C
      IF ( NOIT.GT.NOITM ) THEN
C
        WRITE ( LULOG, 22 ) RM, GRMX, DOSC
 22   FORMAT('Unsolved: RM= ',1PD16.7,' Real= ',1PD16.7,
     1         ' Imag= ',1PD16.7)
C
        DRMC = RM
        IF ( DOSC.NE.DZERO ) THEN
          CHST = 'o'
          CALL EVECEX( N2, NCE, NEGRP, V, VECRE )
          CALL EVECEX( N2, NCE, NEGIP, V, VECIM )
        ELSE
          CHST = 's'
          CALL EVECEX( N2, NCE, NEGRP, V, VECRE )
        ENDIF
        GOTO 60
      ENDIF
C
C Choose value for Magnetic Reynolds number
C
      IF ( NOIT.EQ.1 ) RM = DRM1
      IF ( NOIT.EQ.2 ) RM = DRM2
      IF ( NOIT.GT.2 ) RM = X
C
C Zero matrix A
C
      I = 0
      CALL MATOP( A, DZERO, N1, N2, I )
C
      FAC   = RM
      NDRVS = 2
      CALL RV0MFA( NR, N1, N2, KL, KL, KL, IMF, ILN, IRN,
     1    ILN, IRN, INARR, MHT, MHL, MHM, MHP, MHT, MHL, MHM, IN0,
     2    MT0, ML0, MM0, MP0, NBN, NDCS, NDRVS, NDRVM, NCFM, NBN0,
     3    NDCS0, NDRVS, NDRVM, NCFM0, A, FAC, XARR, VEC0, SVFDC,
     4    SVFDC0, NVI, IHNALP, IHNBET, IHNGAM, TVHI, CVI )
C
C Add diffusion coefficients
C
      FAC = 1.0d0
      I   = 4
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, I, MHT,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILN, IRN, SVFDC, XARR, NDCS, A, N1,
     3           N2, IMF, KL, KL, KL )
C
      I   = 5
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, I, MHT,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILN, IRN, SVFDC, XARR, NDCS, A, N1,
     3           N2, IMF, KL, KL, KL )
C
      CMTD = 1.0d0
      CALL MFSEPS( NR, INARR, MHT, MHL, MHM, MHP, NBN, KL,
     1             NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR,
     2             DIAGEL, MHIBC, MHOBC, NEV, NCV, NCVM, MXIT, NCE,
     3             DRSV, SELECT, DR, DI, D3, WORKEV, WORKD, WORKL,
     4             IPIV, ARTOL, VECRE, V, VECIM, CMTD )
C
C We should now have solved the eigenvalue problem.
C
      IF ( NCE.LT.1 ) THEN
        PRINT *,' Subroutine KDCMRN.'
        PRINT *,' NCE = ', NCE
        PRINT *,' Problem with eigensolver.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      WRITE ( LULOG, 17 )
      WRITE ( LULOG, 16 ) NOIT, NCE
      WRITE ( LULOG, 19 ) RM
      WRITE ( LULOG, 17 )
      GRMX   = -1.0d6
      NEGRP  = -1
      NEGIP  = -1
      FIRST  = .TRUE.
      DO I = 1, NCE
        NEG  = -1
        WRITE ( LULOG, 18 ) I, DR( I ), DI( I ), D3( I )
        IF ( DR( I ).GE.GRMX ) THEN
          NEG   = I
          GRMX  = DR( I )
          DOSC  = DI( I )
        ENDIF
        IF ( NEG.GE.0 ) THEN
          IF ( DABS( DOSC ).GE.DZERO ) THEN
            IF ( FIRST ) THEN
             NEGRP = NEG
             NEGIP = NEG + 1
             FIRST = .FALSE.
            ELSE
             FIRST = .TRUE.
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C
 16   FORMAT('Iteration ',I4,': ',I3,' eigenvalues converged.')
 19   FORMAT('Rm value = ',f20.7)
 17   FORMAT('----------------------------------------------')
 18   FORMAT('Eval ',I2,' (',f16.7,',',f16.7,') Res = ',f16.7)
C
C OK - we are ready
C
      IF ( DABS( GRMX ).LT.DTOL ) THEN
C
        WRITE ( LULOG, 21 ) RM, GRMX, DOSC
 21   FORMAT('RM= ',1PD16.7,' Real= ',1PD16.7,' Imag= ',1PD16.7)
C
        DRMC = RM
        IF ( DOSC.NE.DZERO ) THEN
          CHST = 'O'
          CALL EVECEX( N2, NCE, NEGRP, V, VECRE )
          CALL EVECEX( N2, NCE, NEGIP, V, VECIM )
        ELSE
          CHST = 'S'
          CALL EVECEX( N2, NCE, NEGRP, V, VECRE )
        ENDIF
        GOTO 60
      ENDIF
C
C We have not converged: so let's reiterate
C First iteration
C
      IF ( NOIT.EQ.1 ) THEN
        FF1 = GRMX
      ENDIF
C
C Second iteration
C
      IF ( NOIT.EQ.2 ) THEN
        FF2   = GRMX
        VAL1  = FF1*FF2
C
C val1 will be positive if growth rate at both limits
C is of the same sign - return if this is the case.
C
        IF ( VAL1.GE.DZERO ) THEN
          CHST = 'N'
          GOTO 60
        ENDIF
C
C OK - we appear to have Critical Reynolds number between
C our guesses - so lets set out our low and high values
C of x and f(x) ....
C
        IF ( FF1.LT.DZERO ) THEN
          FL = FF1
          XL = DRM1
          FH = FF2
          XH = DRM2
        ELSE
          FL = FF2
          XL = DRM2
          FH = FF1
          XH = DRM1
        ENDIF
C
      ENDIF
C
C Past 2 iterations
C
      IF ( NOIT.GT.2 ) THEN
        IF ( GRMX.LT.DZERO ) THEN
          XL = X
          FL = GRMX
        ELSE
          XH = X
          FH = GRMX
        ENDIF
      ENDIF
C
C Now set our new value of X according to
C the values of XL, FL, XH and FH
C
      IF ( NOIT.GE.2 ) THEN
C
C  In principle, the following scenario will not occur ???
C
        IF ( DABS( XL - XH ).LT.DLOW ) THEN
          PRINT *,' Subroutine KDCMRN.'
          PRINT *,' XL = ',XL,' XH = ', XH
          PRINT *,' This will lead to an undefined DM.'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
C  In principle, the following scenario will not occur ???
C
        IF ( DABS( FL - FH ).LT.DLOW ) THEN
          PRINT *,' Subroutine KDCMRN.'
          PRINT *,' FL = ',FL,' FH = ', FH
          PRINT *,' This will lead to an undefined DM.'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
        DM = (FH-FL)/(XH-XL)
C
        DC = FL - DM*XL
C
        X = (-1.0d0)*DC/DM
C
      ENDIF
C
      GOTO 50
 60   CONTINUE
C     .
      RETURN
      END
C*********************************************************************

