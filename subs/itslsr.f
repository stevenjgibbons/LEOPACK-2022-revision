C*********************************************************************
C subroutine Inhomogeneous Thermal boundary Solution Linear          *
C            -             -                -        -               *
C                                              Stability Routine     *
C                                              -         -           *
C Steve Gibbons Wed Jun  7 12:48:21 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Using BLCNRS (or equivalently), we have solved for a steady flow   C
C contained in the vector VEC0, indexed by IN0, MHT0, MHL0, MHM0,    C
C MHP0. In addition, the temperature T is given by                   C
C                                                                    C
C  T = \Theta + T_0, where T_0 is a temperature which satisfies      C
C the boundary conditions, and \Theta has homogeneous boundaries.    C
C                                                                    C
C For each harmonic, IH, MHI0( IH ) contains an integer number,      C
C IITH. Then there are 3 numbers, ca, cb and cc with                 C
C                                                                    C
C  CA is stored in CAFIT( 1, IITH )                                  C
C  CB is stored in CAFIT( 2, IITH )                                  C
C  CC is stored in CAFIT( 3, IITH )                                  C
C                                                                    C
C and then the radial function associated with harmonic IH which     C
C gives the inhomogeneous temperature is                             C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C                                                                    C
C An instability with indices INARR, MHT, MHL, MHM and MHP is solved C
C for with the eigensystem. (Note these harmonic sets must be        C
C pre-calculated).                                                   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Our equation is                                                    C
C                                                                    C
C    CA sigma T =     CD \nabla^2 T                                  C
C                     - CC v . \nabla ( T_0 )                        C
C                     - CC v_0 . \nabla ( T )                        C
C                     + v . ( CB1 r + CB2 r^{-2} , 0 , 0 )           C
C                                                                    C
C    CE sigma T =     CI \nabla^2 \curl v                            C
C                     - CG \curl ( k \times v )                      C
C                     + CH \curl (    T   {\bm r } )                 C
C                     + CF \curl ( v \times \curl v_0 )              C
C                     + CF \curl ( v_0 \times \curl v )              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IN0       : Int. parameter array corresponding to VEC0.        C
C                                                                    C
C                 IN0( 1 ) = IFORMF                                  C
C                 IN0( 2 ) = NRR      See INDFUN for details         C
C                 IN0( 3 ) = NH0                                     C
C                                                                    C
C     MHT0     : Array length ( * ) - atleast length NH0             C
C                                                                    C
C         MHT0( IH ) = 1 --> harmonic is poloidal velocity           C
C         MHT0( IH ) = 2 --> harmonic is toroidal velocity           C
C         MHT0( IH ) = 3 --> harmonic is temperature.                C
C         MHT0( IH ) = 4 --> harmonic is poloidal magnetic field     C
C         MHT0( IH ) = 5 --> harmonic is toroidal magnetic field     C
C                                                                    C
C     MHL0     : Array length ( * ) - atleast length NH0             C
C                  Sph. harm. degree, l.                             C
C     MHM0     : Array length ( * ) - atleast length NH0             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHI0     : Dim (*). See above. Inhomog. boundary array.        C
C                                                                    C
C     MHP0     : Array length ( * ) - atleast length NH0             C
C                  Pointer array to finite difference coefficients.  C
C                  MPI( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C     N1        : First dimension of A matrix. Must equal 3*KL+1     C
C     N2        : Second dimension of A matrix. Length of vector.    C
C     KL        : Number of diagonal elements in A matrix.           C
C                                                                    C
C     MHT      : Array length ( * ) - atleast length NH              C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. degree, l.                             C
C     MHM      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHP      : Array length ( * ) - atleast length NH              C
C                  Pointer array to finite difference coefficients.  C
C                  MPI( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     MHTR     : Array length ( * ) - atleast length NH              C
C                 Type of function in equation rows.                 C
C                 Under normal circumstances, MHTR is formed by      C
C                                                                    C
C                 CALL CINDSW( NH, MHT, MHTR )                       C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NCFM      : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C      (NDRVM must be atleast 4 and NDRVS must be 4 for atleast      C
C       grid nodes IR = 2, NR - 1 ).                                 C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     NTHP      : The number of theta points.                        C
C     NPHP      : The number of phi points.                          C
C                                                                    C
C     MMAX      : Maximum sph. harmonic order, m.                    C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
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
C     IPIV      : Dim ( N2 ). Array for pivotting.                   C
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
C     NCE       : Number of converged eigenvalues.                   C
C     IEV       : Eigenvalue with largest real part.                 C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC0      : Solution vector, v_0. Dim (*)                      C
C                                                                    C
C     PARAM     : Array dim ( 10 ) containing these param.s          C
C                                                                    C
C  On input:                                                         C
C                                                                    C
C            PARAM(  1 ) = CA                                        C
C            PARAM(  2 ) = CB1                                       C
C            PARAM(  3 ) = CB2                                       C
C            PARAM(  4 ) = CC                                        C
C            PARAM(  5 ) = CD                                        C
C            PARAM(  6 ) = CE                                        C
C            PARAM(  7 ) = CF                                        C
C            PARAM(  8 ) = CG                                        C
C            PARAM(  9 ) = CH                                        C
C            PARAM( 10 ) = CI                                        C
C                                                                    C
C     A         : Work array. Dim ( N1, N2 )                         C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     CAFIT     : Coeff arr. for inhomog. temp. Dim( 3, NITHMX ).    C
C                                                                    C
C     SBRVEC    : Work array - dim. (N2) Critical eigenvector is     C
C                 returned in SBRVEC (real part)                     C
C                                                                    C
C     DR        : Dim (NCVM). Real parts of eigenvalues.             C
C     DI        : Dim (NCVM). Imag parts of eigenvalues.             C
C     D3        : Dim (NCVM). Direct residuals.                      C
C                                                                    C
C     DRSV      : Real part of shift.                                C
C     ARTOL     : Convergence criterion for ARPACK.                  C
C                                                                    C
C     WORKEV    : Dim (3*NCVM). Work array.                          C
C     WORKD     : Dim (3*N2).  Work array.                           C
C     RESID     : Dim (N2).  Work array.                             C
C     WVEC      : Dim (N2). Work array.                              C
C     WORKL     : Dim (3*NCVM*NCVM + 6*NCVM). Work array.            C
C     V         : Dim (N2,NCVM). Contains eigenvectors on return.    C
C     W2        : Work array - dim. (N2) Critical eigenvector is     C
C                 returned in SBRVEC (imag part)                     C
C                                                                    C
C     GRR       : Greatest real eigenvalue part.                     C
C     GRI       : Greatest imag eigenvalue part.                     C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     SELECT    : Dimension ( NCVM ). Work array.                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITSLSR( IN0, MHT0, MHL0, MHM0, MHI0, MHP0, NR, N1,
     1  N2, KL, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, NCFM, NDRVM,
     2  NDCS, NTHP, NPHP, MMAX, LH, MHIBC, MHOBC, IPIV, NEV, NCV,
     3  NCVM, MXIT, VEC0, PARAM, A, SVFDC, XARR, GAUX, GAUW, PA, DPA,
     4  CAFIT, SBRVEC, DR, DI, D3, DRSV, ARTOL, WORKEV, WORKD, RESID,
     5  WVEC, WORKL, V, W2, SELECT, NCE, IEV, GRR, GRI )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IN0( * ), MHT0( * ), MHL0( * ), MHM0( * ), MHP0( * ),
     1        MHI0( * ), NR, N1, N2, KL, MHT( * ), MHM( * ), IEV,
     2        MHP( * ), MHTR( * ), NBN, NCFM, NDRVM, NDCS, NTHP,
     3        NPHP, MMAX, LH, MHIBC( NDCS ), MHOBC( NDCS ), NCE,
     4        IPIV( * ), NEV, NCV, NCVM, MXIT, MHL( * ), INARR( * )
      DOUBLE PRECISION VEC0( * ), PARAM( * ), A( N1, N2), SBRVEC(*),
     1        DR( NCVM ), DI( NCVM ), D3( NCVM ), DRSV, ARTOL,
     2        WORKEV( 3*NCVM ), WORKD( 3*N2 ), RESID( N2 ), WVEC( N2),
     3        WORKL( 3*NCVM*NCVM + 6*NCVM ), V( N2, NCVM ), W2( N2 )
      DOUBLE PRECISION CAFIT( 3, * ), GAUX( NTHP ), GAUW( NTHP ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 SVFDC( NCFM, NR, NDRVM+1, NDCS ), XARR( * )
      DOUBLE PRECISION GRR, GRI
      LOGICAL SELECT( NCVM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
C Array defining integers:
C
      INTEGER NPHMAX, NTHMAX, LHMAX, MAXNVI
      PARAMETER ( NPHMAX = 128, NTHMAX = 64, LHMAX = 62,
     1            MAXNVI = 200000 )
C
      INTEGER KKA01( MAXNVI ), KKB01( MAXNVI ), KKG01( MAXNVI ),
     1        KKA10( MAXNVI ), KKB10( MAXNVI ), KKG10( MAXNVI ),
     2        KKA02( MAXNVI ), KKB02( MAXNVI ), KKG02( MAXNVI ),
     3        KKA20( MAXNVI ), KKB20( MAXNVI ), KKG20( MAXNVI )
C
      DOUBLE PRECISION PAR8( 8 ), F1( 2*NPHMAX ), F2( 2*NPHMAX ),
     1                 F3( 2*NPHMAX ), SHC( LHMAX*(LHMAX + 2) ),
     2                 CVI10( MAXNVI ), CVI20( MAXNVI ),
     3                 CVI01( MAXNVI ), CVI02( MAXNVI )
C
      CHARACTER *(4) TVHI01( MAXNVI ), TVHI02( MAXNVI ),
     1               TVHI10( MAXNVI ), TVHI20( MAXNVI )
C
      DOUBLE PRECISION QST( LHMAX*(LHMAX + 2), 3 ),
     1                 VF1( NPHMAX, NTHMAX, 3),
     2                 VF2( NPHMAX, NTHMAX, 3),
     3                 VF3( NPHMAX, NTHMAX, 3),
     4                 SF( NPHMAX, NTHMAX )
C
C ------------------------------
C Ordinary working variables ...
C ------------------------------
C
      INTEGER NH, ICLS, ILN, IRN, ILNT, IRNT, NVI10, NVI01,
     1        NVI20, NVI02, NH0, IMF, NDRVS, ITGN
      DOUBLE PRECISION ZERO, CA, CB1, CB2, CC, CD, CE, CF, CG, CH,
     1        CI, DIAGEL, FAC, RESTOL, DSRSV
      LOGICAL OSBR
      PARAMETER ( ZERO = 0.0d0, ICLS = 1, DIAGEL = -200.0d0,
     1            IMF = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Unload initial values of PARAM ...
C
      CA      = PARAM(  1 )
      CB1     = PARAM(  2 )
      CB2     = PARAM(  3 )
      CC      = PARAM(  4 )
      CD      = PARAM(  5 )
      CE      = PARAM(  6 )
      CF      = PARAM(  7 )
      CG      = PARAM(  8 )
      CH      = PARAM(  9 )
      CI      = PARAM( 10 )
C
      NDRVS = 4
      ITGN  = NR/2
C
      ILN   = 2
      IRN   = NR - 1
      ILNT  = 3
      IRNT  = NR - 2
C
      NH    = INARR( 3 )
      NH0   = IN0( 3 )
C
      IF ( NPHP.GT.NPHMAX .OR. LH.GT.LHMAX .OR. NTHP.GT.NTHMAX ) THEN
        PRINT *,' Subroutine ITSLSR.'
        PRINT *,' LH   = ', LH,  ' LHMAX  = ', LHMAX
        PRINT *,' NPHP = ', NPHP,' NPHMAX = ', NPHMAX
        PRINT *,' NTHP = ', NTHP,' NTHMAX = ', NTHMAX
        PRINT *,' Recompile routine with higher dimensions.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      PAR8(  1 ) = ZERO
      PAR8(  2 ) = CB1
      PAR8(  3 ) = CB2
      PAR8(  4 ) = CD
      PAR8(  5 ) = ZERO
      PAR8(  6 ) = CG
      PAR8(  7 ) = CH
      PAR8(  8 ) = CI
C
C We now calculate the vector interactions which
C are used to form the matrix non-linear terms.
C
      NVI10 = 0
      NVI01 = 0
      NVI20 = 0
      NVI02 = 0
C
C First: Scalar product ( v_0 . \nabla \Theta )
C
      CALL VSPCC( NVI01, MAXNVI, KKA01, KKB01, KKG01, NH0, MHT0, MHL0,
     1            MHM0, NH, MHT, MHL, MHM, NH, MHTR, MHL, MHM, LH,
     2            NTHP, NPHP, MMAX, TVHI01, CVI01, GAUX, GAUW, PA,
     3            DPA, F1, VF1, VF2, SF, SHC )
C
C Then: Scalar product ( v . \nabla \Theta_0 )
C
      CALL VSPCC( NVI10, MAXNVI, KKA10, KKB10, KKG10, NH, MHT, MHL,
     1            MHM, NH0, MHT0, MHL0, MHM0, NH, MHTR, MHL, MHM, LH,
     2            NTHP, NPHP, MMAX, TVHI10, CVI10, GAUX, GAUW, PA,
     3            DPA, F1, VF1, VF2, SF, SHC )
C
C Now: Inertial term ( - \curl ( v_0 \times \curl v )   )
C
      CALL VCCPCC( NVI02, MAXNVI, KKA02, KKB02, KKG02, NH0, MHT0, MHL0,
     1             MHM0, NH, MHT, MHL, MHM, NH, MHTR, MHL, MHM, LH,
     2             NTHP, NPHP, MMAX, TVHI02, CVI02, GAUX, GAUW, PA, DPA,
     3             F1, F2, F3, VF1, VF2, VF3, QST )
C
C Now: Inertial term ( - \curl ( v \times \curl v_0 )   )
C
      CALL VCCPCC( NVI20, MAXNVI, KKA20, KKB20, KKG20, NH, MHT, MHL,
     1             MHM, NH0, MHT0, MHL0, MHM0, NH, MHTR, MHL, MHM, LH,
     2             NTHP, NPHP, MMAX, TVHI20, CVI20, GAUX, GAUW, PA, DPA,
     3             F1, F2, F3, VF1, VF2, VF3, QST )
C
C Form the linear stability matrix.
C
      CALL AVMLTA( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1             NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR,
     2             NTHP, NPHP, MMAX, LH, GAUX, GAUW, PA, DPA,
     3             F1, F2, F3, VF1, QST, PAR8, ICLS )
C
C Now add the non-linear terms to the matrix.
C We have already calculated the coefficients for this
C using VSPCC and VCCPCC ...
C
C Firstly: curl ( v0 x curl v_new ) terms
C
      FAC   = (-1.0d0)*CF
      CALL RV0CVA( NR, N1, N2, KL, KL, KL, IMF, ILN, IRN,
     1    ILNT, IRNT, INARR, MHT, MHL, MHM, MHP, MHTR, MHL, MHM,
     2    IN0, MHT0, MHL0, MHM0, MHP0, NBN, NDCS, NDRVS, NDRVM, NCFM,
     3    NBN, NDCS, NDRVS, NDRVM, NCFM, A, FAC, XARR, VEC0, SVFDC,
     4    SVFDC, NVI02, KKA02, KKB02, KKG02, TVHI02, CVI02 )
C
C Now: curl ( v_new x curl v0 ) terms
C
      CALL RVCV0A( NR, N1, N2, KL, KL, KL, IMF, ILN, IRN,
     1    ILNT, IRNT, INARR, MHT, MHL, MHM, MHP, MHTR, MHL, MHM,
     2    IN0, MHT0, MHL0, MHM0, MHP0, NBN, NDCS, NDRVS, NDRVM, NCFM,
     3    NBN, NDCS, NDRVS, NDRVM, NCFM, A, FAC, XARR, VEC0, SVFDC,
     4    SVFDC, NVI20, KKA20, KKB20, KKG20, TVHI20, CVI20 )
C
C Now: v0 . nabla ( THETA_new ) terms
C
      FAC   = (-1.0d0)*CC
      CALL RV0GTA( NR, N1, N2, KL, KL, KL, IMF, ILN, IRN, INARR,
     1    MHT, MHL, MHM, MHP, MHTR, MHL, MHM, IN0, MHT0, MHL0, MHM0,
     2    MHP0, NBN, NDCS, NDRVS, NDRVM, NCFM, NBN, NDCS, NDRVS,
     3    NDRVM, NCFM, A, FAC, XARR, VEC0, SVFDC, SVFDC,
     4    NVI01, KKA01, KKB01, KKG01, TVHI01, CVI01 )
C
C Now: v_new . nabla ( THETA_0 ) terms
C
      CALL RVGI0A( NR, N1, N2, KL, KL, KL, IMF, ILN, IRN, INARR,
     1    MHT, MHL, MHM, MHP, MHTR, MHL, MHM, IN0, MHT0, MHL0, MHM0,
     2    MHP0, NBN, NDCS, NDRVS, NDRVM, NCFM, NBN, NDCS, NDRVS,
     3    NDRVM, NCFM, A, FAC, XARR, VEC0, SVFDC, SVFDC,
     4    NVI10, KKA10, KKB10, KKG10, TVHI10, CVI10, MHI0, CAFIT )
C
C Solid body rotation solution pre-requisite
C We can use IPIV, W2, WVEC, RESID as work arrays -
C they are plenty large enough!
C
      CALL SBRRFC( NR, INARR, MHT, MHL, MHM, MHP, MHIBC, MHOBC,
     1             NCFM, XARR, SBRVEC, OSBR, W2, WVEC,
     2             IPIV, RESID )
C
C Solve the eigensystem
C
      CALL VMEPS( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1       NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR, OSBR, SBRVEC,
     2       DIAGEL, MHIBC, MHOBC, ITGN, NEV, NCV, NCVM, MXIT, NCE,
     3       DRSV, SELECT, DR, DI, D3, WORKEV, WORKD, WORKL, IPIV,
     4       ARTOL, RESID, V, W2, WVEC, CA, CE )
C
      IF ( NCE.LT.1 ) THEN
        PRINT *,' Subroutine ITSLSR.'
        PRINT *,' NCE = ', NCE
        PRINT *,' Problem with eigensolver.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RESTOL = 1.0d-5
      CALL EVALAS( NCE, IEV, DRSV, RESTOL, DR, DI, D3,
     1                   GRR, GRI, DSRSV )
C
      CALL EVECEX( N2, NCVM, IEV, V, SBRVEC )
      CALL EVECEX( N2, NCVM, IEV+1, V, W2 )
C
      RETURN
      END
C*********************************************************************
