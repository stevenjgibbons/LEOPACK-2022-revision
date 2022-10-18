C*********************************************************************
C subroutine Instability Onset Critical Rayleigh Bounds refine (Opt) *
C            -           -     -        -        -              -    *
C Steve Gibbons Tue Nov  6 12:00:09 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Seeks a critical Rayleigh number for the onset of convective       C
C instability when subjected to a flow and temperature field         C
C ( v_0 , \Theta_0 ).                                                C
C                                                                    C
C The equations are:-                                                C
C                                                                    C
C  c_a d \Theta/ dt = CD \nabla^2 \Theta                             C
C                     + v . ( CB1 r + CB2 r^{-2} , 0 , 0 )           C
C                     - CC1 v_0 . Grad ( \Theta )                    C
C                     - CC2 v . Grad ( \Theta_0 )                    C
C                                                                    C
C  c_e \curl dv/dt  = CI \nabla^2 \curl v                            C
C                     - CG \curl ( k \times v )                      C
C                     + CH \curl ( \Theta {\bm r } )                 C
C                     - CF1 \curl ( v_0 . \nabla ) v                 C
C                     - CF2 \curl ( v . \nabla ) v_0                 C
C                                                                    C
C Two bounds are taken, x1 and x2 within which the Critical          C
C Rayleigh number must lie. If the chief growth rate (GRR) at both   C
C these points has the same sign, then it is assumed that R_c        C
C does not lie in this interval and we return with IERR = -2.        C
C                                                                    C
C Otherwise, we iterate between these points.                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARR     : Int. parameter array corresponding to VI.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR      See INDFUN for details       C
C                 INARR( 3 ) = NH                                    C
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
C     KL        : Number of lower diagonals in matrix                C
C                                                                    C
C     NCFM      : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C      (NDRVM must be atleast 4 and NDRVS must be 4 for atleast      C
C       grid nodes IR = 2, NR - 1 ).                                 C
C                                                                    C
C     N1        : First dimension of matrix, A.                      C
C     N2        : Second dimension of matrix, A.                     C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     NTHP      : The number of theta points.                        C
C     NPHP      : The number of phi points.                          C
C                                                                    C
C     MMAX      : Value of highest sph. harmonic order, m.           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                                                                    C
C     IEV       : Number of the critical eigenvalue.                 C
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
C     LULOG     : Logical unit number of log file.                   C
C                LULOG = 0 if no log file is to be written.          C
C                LULOG = 44 or 144, only CH, GRR and GRI are written C
C                   out at each iteration.                           C
C                LULOG = 45 or 145, full eigensystem written out.    C
C                                                                    C
C     IERR      : Error flag on return.                              C
C                 If IERR greater than zero, all is well and         C
C                 IERR contains the number of iterations required    C
C                 for the growth rate to converge.                   C
C                                                                    C
C                 if IERR = -1, more than the maximum number of      C
C                 iterations were needed to find R_c.                C
C                                                                    C
C                 if IERR = -2, the growth rates at both limits      C
C                 have the same sign and so it is assumed that       C
C                 there is no point with zero growth rate.           C
C                                                                    C
C     IN0       : Equivalent to INARR for existing solution ( * ).   C
C     MT0       : Equivalent to MHT for existing solution ( * ).     C
C     ML0       : Equivalent to MHL for existing solution ( * ).     C
C     MM0       : Equivalent to MHM for existing solution ( * ).     C
C     MP0       : Equivalent to MHP for existing solution ( * ).     C
C     M0        : Minimum non-zero wavenumber.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 )                             C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     SBRVEC    : Work array - dim. (N2) Critical eigenvector is     C
C                 returned in SBRVEC (real part)                     C
C     RESID     : Work array - dim. (N2)                             C
C     W2        : Work array - dim. (N2) Critical eigenvector is     C
C                 returned in SBRVEC (imag part)                     C
C     WVEC      : Work array - dim. (N2)                             C
C     F1        : Work array - dim. (2*NPHP)                         C
C     F2        : Work array - dim. (2*NPHP)                         C
C     F3        : Work array - dim. (2*NPHP)                         C
C     VF1       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF2       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF3       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     QST       : Work array - dim. ( LH*(LH+2), 3)                  C
C     SF        : Work array - dim. ( NPHP, NTHP )                   C
C     SHC       : Work array - dim. ( LH*(LH+2) )                    C
C                                                                    C
C     PARAM     : Array dim ( 14 ) containing these param.s          C
C                                                                    C
C  On input:                                                         C
C                                                                    C
C            PARAM(  1 ) = CA                                        C
C            PARAM(  2 ) = CB1                                       C
C            PARAM(  3 ) = CB2                                       C
C            PARAM(  4 ) = CD                                        C
C            PARAM(  5 ) = CE                                        C
C            PARAM(  6 ) = CG                                        C
C            PARAM(  7 ) = CH1 - First value of CH to try            C
C            PARAM(  8 ) = CH2 - Second value of CH to try           C
C            PARAM(  9 ) = CI                                        C
C            PARAM( 10 ) = CTOL                                      C
C            PARAM( 11 ) = CC1                                       C
C            PARAM( 12 ) = CC2                                       C
C            PARAM( 13 ) = CF1                                       C
C            PARAM( 14 ) = CF2                                       C
C                                                                    C
C     CTOL      : Convergence tolerance for growth rate              C
C                 If abs( grr ).lt.ctol, we assume that the sol.     C
C                 is converged.                                      C
C                                                                    C
C  On output:                                                        C
C                                                                    C
C            PARAM(  7 ) = Critical Rayleigh number                  C
C            PARAM(  8 ) = CI                                        C
C            PARAM(  9 ) = GRR at Critical value.                    C
C            PARAM( 10 ) = GRI at Critical value.                    C
C            PARAM( 11 ) = CC1                                       C
C            PARAM( 12 ) = CC2                                       C
C            PARAM( 13 ) = CF1                                       C
C            PARAM( 14 ) = CF2                                       C
C                                                                    C
C     VEC0      : Vector containing ( v_0, \Theta_0 )                C
C     VEC0M     : Work array required to store modified VEC0.        C
C                                                                    C
C     REY       : Multiplying factor for the imposed vector.         C
C     REYSB1    : Initial scaling for solid body rotation.           C
C     REYSB2    : Second scaling for solid body rotation.            C
C                                                                    C
C  (Final solid body rotation is returned in REYSB1)                 C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     SELECT    : Dimension ( NCVM ). Work array.                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE IOCRBD( NR,INARR,MHT,MHL,MHM,MHP,MHTR,NBN,KL,NCFM,
     1   NDRVM,N1,N2,NDCS,NTHP,NPHP,MMAX,LH,IEV,MHIBC,MHOBC,NEV,NCV,
     2   NCVM,MXIT,IPIV,LULOG,SVFDC,A,XARR,GAUX,GAUW,PA,DPA,SBRVEC,
     3   RESID,W2,WVEC,F1,F2,F3,VF1,VF2,VF3,QST,SF,SHC,SELECT,ARTOL,
     4   DRSV,PARAM,MXATT,IERR,DR,DI,D3,WORKEV,WORKD,WORKL,V,
     5   IN0,MT0,ML0,MM0,MP0,VEC0,VEC0M,M0,REY,REYSB1,REYSB2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( 3 ), MHT( * ), MHL( * ), MHM( * ), MHP( * ),
     1        MHTR( * ), NBN, KL, NCFM, NDRVM, N1, N2, NDCS, M0,
     2        MHIBC( NDCS ), MHOBC( NDCS ), NEV, NCV, NCVM, MXIT,
     3        IPIV( * ), LULOG, MMAX, LH, IEV, NTHP, NPHP, MXATT,
     4        IERR, IN0( * ), MT0( * ), ML0( * ), MM0( * ), MP0( * )
      DOUBLE PRECISION A( N1, N2 ), SBRVEC( * ), XARR( NR ),
     1                 DR( NCVM ), DI( NCVM ), D3( NCVM ), DRSV,
     2                 WORKEV( 3*NCVM ), WORKD( 3*N2 ),
     3                 WORKL( 3*NCVM*NCVM + 6*NCVM ), ARTOL
      DOUBLE PRECISION RESID( N2 ), V( N2, NCVM ), WVEC( N2 ),
     1                 W2( N2 ), SVFDC( NCFM, NR, NDRVM+1, NDCS ),
     2                 PARAM( 14 ), SF( NPHP, NTHP ), VEC0( * ),
     3                 VEC0M( * ), REY, REYSB1, REYSB2
      DOUBLE PRECISION
     1                 QST( LH*(LH+2), 3), SHC( LH*(LH+2) ),
     2                 VF1( NPHP, NTHP, 3), VF2( NPHP, NTHP, 3),
     3                 VF3( NPHP, NTHP, 3), F1( 2*NPHP ),
     4                 F2( 2*NPHP ), F3( 2*NPHP )
      DOUBLE PRECISION
     1                 GAUX( NTHP ), GAUW( NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2, NTHP ),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2, NTHP )
      LOGICAL SELECT( NCVM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NCE, ICLS, IH, NOIT, ITGN
      DOUBLE PRECISION CA, CE, CB1, CB2, CD, CG, CH, CI, RESTOL,
     1                 GRR, GRI, DSRSV, DIAGEL, DELTA, ZERO, ROT,
     2                 CTOL, CH1, CH2, FF1, FF2, VAL1, XL, FL,
     3                 XH, FH, X, DM, DC, CC1, CC2, CF1, CF2
      DOUBLE PRECISION DIFF, DIFFM, DMVAL, DRFTNM, DRFTNP, DRFTOM,
     1                 DRFTOP, GRIO, REYSBO, REYSBR         
      PARAMETER ( ICLS = 1, DIAGEL = -200.0d0, DELTA = 1.0d-8,
     1            ZERO = 0.0d0 )
      LOGICAL OSBR
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      ITGN    = NR/2
      DMVAL   = DBLE( M0 )
C
C Unload initial values of PARAM ...
C
      CA      = PARAM(  1 )
      CB1     = PARAM(  2 )
      CB2     = PARAM(  3 )
      CD      = PARAM(  4 )
      CE      = PARAM(  5 )
      ROT     = PARAM(  6 )
      CH1     = PARAM(  7 )
      CH2     = PARAM(  8 )
      CI      = PARAM(  9 )
      CTOL    = PARAM( 10 )
      CC1     = PARAM( 11 )
      CC2     = PARAM( 12 )
      CF1     = PARAM( 13 )
      CF2     = PARAM( 14 )
C
      IF ( CTOL.LT.ZERO ) THEN
        PRINT *,' Subroutine IOCRBD.'
        PRINT *,' CTOL = ', CTOL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LULOG.NE.0 .AND. LULOG.NE.144 .AND. LULOG.NE.145 .AND.
     1         LULOG.NE.44 .AND. LULOG.NE.45 ) THEN
        PRINT *,' Subroutine IOCRBD.'
        PRINT *,' LULOG = ', LULOG
        PRINT *,' Should be 0, 44, 45, 144 or 145.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Begin iterations to find critical Rayleigh number
C
      NOIT = 0
 50   CONTINUE
      NOIT = NOIT + 1
C
      IF ( NOIT.GT.MXATT ) THEN
        IERR = -1
        GOTO 60
      ENDIF
C
C Choose value of C according to whether we are
C on 1st, 2nd or later iteration
C
      IF ( NOIT.EQ.1 ) THEN
        CH     = CH1
        REYSBO = 0.0d0
        REYSBR = REYSB1
      ENDIF
C
      IF ( NOIT.EQ.2 ) THEN
        CH     = CH2
        REYSBO = REYSB1
        REYSBR = REYSB2
        GRIO   = GRI
      ENDIF
C
      IF ( NOIT.GT.2 ) THEN
        CH     = X
        GRIO   = GRI
      ENDIF
C
C Form modified imposed vector
C
      CALL SBRVMR( IN0, MT0, ML0, MM0, VEC0, VEC0M, XARR,
     1             REY, REYSBR )
C
      CG = ROT - 2.0d0*REYSBR*CE
C
      PARAM(  1 ) = ZERO
      PARAM(  2 ) = CB1
      PARAM(  3 ) = CB2
      PARAM(  4 ) = CD
      PARAM(  5 ) = ZERO
      PARAM(  6 ) = CG
      PARAM(  7 ) = CH
      PARAM(  8 ) = CI
      PARAM(  9 ) = CC1
      PARAM( 10 ) = CC2
      PARAM( 11 ) = CF1
      PARAM( 12 ) = CF2
C
C Matrix must be zero'd on every iteration
C
C Form the linear stability matrix
C
      CALL AVMATO( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1          NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR, NTHP,
     2          NPHP, MMAX, LH, GAUX, GAUW, PA, DPA, F1, F2, F3,
     3          VF1, VF2, VF3, QST, SF, SHC, PARAM, ICLS,
     4          IN0, MT0, ML0, MM0, MP0, VEC0M, M0 )
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
C Analyse output
C
      IF ( LULOG.EQ.45 .OR. LULOG.EQ.145 ) THEN
        WRITE ( LULOG, 17 ) 
        WRITE ( LULOG, 16 ) NOIT, NCE
        WRITE ( LULOG, *  ) 'NOIT = ',NOIT, ' REYSBR = ',REYSBR
        WRITE ( LULOG, 17 ) 
        DO IH = 1, NCE
        WRITE ( LULOG, 18 ) IH, DR( IH ), DI( IH ), D3( IH )
        ENDDO
      ENDIF
 16   FORMAT('Iteration ',I4,': ',I3,' eigenvalues converged.')
 17   FORMAT('----------------------------------------------')
 18   FORMAT('Eval ',I2,' (',f16.7,',',f16.7,') Res = ',f16.7)
C
      IF ( NCE.LT.1 ) THEN
        PRINT *,' Subroutine IOCRBD.'
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
      IF ( LULOG.NE.0 ) THEN
        WRITE ( LULOG, 20 ) NOIT, CH, GRR, GRI
      ENDIF
 20   FORMAT(' Iter. ',I4,' R= ',1pd16.8,
     1       ': (',1pd16.8,',',1pd16.8,')')
C
      IF ( DSRSV.NE.DRSV ) THEN
        DRSV = DSRSV
        IF ( LULOG.EQ.45 .OR. LULOG.EQ.145 ) WRITE ( LULOG, 19 ) DSRSV
      ENDIF
 19   FORMAT('Real shift changed to ',1pd15.7)
C
C Ok - now we know GRR and GRI
C Firstly, see if we have converged ....
C
      IF ( DABS( REYSBR - REYSBO ).LT.CTOL .AND.
     1              DABS( GRR ).LT.CTOL ) THEN
        REYSB1 = REYSBR
        IERR = NOIT
        GOTO 60
      ENDIF
C
C OK - so we have not converged
C First iteration ...
C
      IF ( NOIT.EQ.1 ) THEN
        FF1 = GRR
      ENDIF
C
C Second iteration ...
C
      IF ( NOIT.EQ.2 ) THEN
        FF2  = GRR
        VAL1 = FF1*FF2
C
C val1 will be positive if growth rate at both limits
C is of the same sign - return if this is the case.
C
        IF ( VAL1.GT.ZERO ) THEN
          IERR = -2
          GOTO 60
        ENDIF
C
C OK - we appear to have Critical Rayleigh number between
C our guesses - so lets set out our low and high values
C of x and f(x) ....
C
        IF ( FF1.LT.ZERO ) THEN
          FL = FF1
          XL = CH1
          FH = FF2
          XH = CH2
        ELSE
          FL = FF2
          XL = CH2
          FH = FF1
          XH = CH1
        ENDIF
C
      ENDIF
C
C Iterations with NOIT.GT.2
C
      IF ( NOIT.GT.2 ) THEN
        IF ( GRR.LT.ZERO ) THEN
          XL = X
          FL = GRR
        ELSE
          XH = X
          FH = GRR
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
        IF ( DABS( XL - XH ).LT.(DELTA) ) THEN
          PRINT *,' Subroutine IOCRBD.'
          PRINT *,' XL = ',XL,' XH = ', XH
          PRINT *,' This will lead to an undefined DM.'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
C  In principle, the following scenario will not occur ???
C
        IF ( DABS( FL - FH ).LT.(DELTA) ) THEN
          PRINT *,' Subroutine IOCRBD.'
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
C  Now set the new value of REYSBR
C
        DRFTOP  = REYSBO + GRIO/DMVAL
        DRFTOM  = REYSBO - GRIO/DMVAL
        DRFTNP  = REYSBR + GRI/DMVAL
        DRFTNM  = REYSBR - GRI/DMVAL
C
        REYSBO  = REYSBR
        DIFF    = DABS( DRFTOP - DRFTNP )
        DIFFM   = DIFF
        REYSBR  = 0.5d0*( DRFTOP + DRFTNP )
        DIFF    = DABS( DRFTOP - DRFTNM )
        IF ( DIFF.LT.DIFFM ) THEN
          DIFFM   = DIFF
          REYSBR  = 0.5d0*( DRFTOP + DRFTNM )
        ENDIF
        DIFF    = DABS( DRFTNP - DRFTOM )
        IF ( DIFF.LT.DIFFM ) THEN
          DIFFM   = DIFF
          REYSBR  = 0.5d0*( DRFTNP + DRFTOM )
        ENDIF
        DIFF    = DABS( DRFTNM - DRFTOM )
        IF ( DIFF.LT.DIFFM ) THEN
          DIFFM   = DIFF
          REYSBR  = 0.5d0*( DRFTNM + DRFTOM )
        ENDIF
C
      ENDIF
C
      GOTO 50
C
 60   CONTINUE
      PARAM(  1 ) = CA
      PARAM(  2 ) = CB1
      PARAM(  3 ) = CB2
      PARAM(  4 ) = CD
      PARAM(  5 ) = CE
      PARAM(  6 ) = ROT
      PARAM(  7 ) = CH
      PARAM(  8 ) = CI
      PARAM(  9 ) = GRR
      PARAM( 10 ) = GRI
      PARAM( 11 ) = CC1
      PARAM( 12 ) = CC2
      PARAM( 13 ) = CF1
      PARAM( 14 ) = CF2
C
      RETURN
      END
C*********************************************************************
