C*********************************************************************
C subroutine Linear Onset Critical Rayleigh Number Find **************
C            -      -     -        -        -      -    **************
C Steve Gibbons Tue Nov 23 19:40:21 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Seeks a critical Rayleigh number for the onset of convection       C
C                                                                    C
C The equations are:-                                                C
C                                                                    C
C  c_a d \Theta/ dt = CD \nabla^2 \Theta                             C
C                     + v . ( CB1 r + CB2 r^{-2} , 0 , 0 )           C
C                                                                    C
C  c_e \curl dv/dt  = CI \nabla^2 \curl v                            C
C                     - CG \curl ( k \times v )                      C
C                     + CH \curl ( \Theta {\bm r } )                 C
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
C     MVAL      : Value of sph. harmonic order, m.                   C
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
C                 LULOG = 0 if no log file is to be written.         C
C                 LULOG = 144 or 44, only CH, GRR and GRI are        C
C                  written out at each iteration.                    C
C                 LULOG = 145 or 45, full eigensystem written out.   C
C                                                                    C
C     IERR      : Error flag on return.                              C
C                 If IERR greater than zero, all is well and         C
C                 IERR contains the number of iterations required    C
C                 for the growth rate to converge.                   C
C                                                                    C
C                 if IERR = -1, more than the maximum number of      C
C                 iterations were needed to find R_c.                C
C                                                                    C
C                 if IERR = -2, a non-positive slope of the          C
C                 growth rate was encountered with no bound on       C
C                 the R_c value.                                     C
C                 This routine assumes that the growth rate          C
C                 increases monotonically. However, if a non-pos.    C
C                 slope of the growth rate is found within bounds    C
C                 the interval is bisected for the next attempt.     C
C                                                                    C
C                 if IERR = -3, the upper limit is lower than the    C
C                 lower limit on R_c. Something is badly wrong!      C
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
C     F1      : Work array - dim. (2*NPHP)                           C
C     F2      : Work array - dim. (2*NPHP)                           C
C     F3      : Work array - dim. (2*NPHP)                           C
C     VF        : Work array - dim. ( NPHP, NTHP, 3)                 C
C     QST       : Work array - dim. ( LH*(LH+2), 3)                  C
C                                                                    C
C     PARAM     : Array dim ( 10 ) containing these param.s          C
C                                                                    C
C  On input:                                                         C
C                                                                    C
C            PARAM(  1 ) = CA                                        C
C            PARAM(  2 ) = CB1                                       C
C            PARAM(  3 ) = CB2                                       C
C            PARAM(  4 ) = CD                                        C
C            PARAM(  5 ) = CE                                        C
C            PARAM(  6 ) = CG                                        C
C            PARAM(  7 ) = CH (Initial estimate)                     C
C            PARAM(  8 ) = CI                                        C
C            PARAM(  9 ) = FACTOR (must be greater than 1.0 )        C
C                                                                    C
C  (Factor is the multiplier which sets the value of CH_{2} )        C
C  If GRR_{1} < 0.0 then CH_{2} = CH_{1} * FACTOR                    C
C  If GRR_{1} > 0.0 then CH_{2} = CH_{1} / FACTOR                    C
C                                                                    C
C            PARAM( 10 ) = CTOL                                      C
C                                                                    C
C     CTOL      : Convergence tolerance for growth rate              C
C                 If abs( grr ).lt.ctol, we assume that the sol.     C
C                 is converged.                                      C
C                                                                    C
C  On output:                                                        C
C                                                                    C
C            PARAM(  7 ) = Critical Rayleigh number                  C
C            PARAM(  9 ) = GRR at Critical value.                    C
C            PARAM( 10 ) = GRI at Critical value.                    C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     SELECT    : Dimension ( NCVM ). Work array.                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LOCRNF( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1   NCFM, NDRVM, N1, N2, NDCS, NTHP, NPHP, MVAL, LH, IEV, MHIBC,
     2   MHOBC, NEV, NCV, NCVM, MXIT, IPIV, LULOG, SVFDC, A, XARR,
     3   GAUX, GAUW, PA, DPA, SBRVEC, RESID, W2, WVEC, F1, F2, F3, VF,
     4   QST, SELECT, ARTOL, DRSV, PARAM, MXATT, IERR, DR, DI, D3,
     5   WORKEV, WORKD, WORKL, V )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( 3 ), MHT( * ), MHL( * ), MHM( * ), MHP( * ),
     1        MHTR( * ), NBN, KL, NCFM, NDRVM, N1, N2, NDCS,
     2        MHIBC( NDCS ), MHOBC( NDCS ), NEV, NCV, NCVM, MXIT,
     3        IPIV( * ), LULOG, MVAL, LH, IEV, NTHP, NPHP, MXATT,
     4        IERR
      DOUBLE PRECISION A( N1, N2 ), SBRVEC( * ), XARR( NR ),
     1                 DR( NCVM ), DI( NCVM ), D3( NCVM ), DRSV,
     2                 WORKEV( 3*NCVM ), WORKD( 3*N2 ),
     3                 WORKL( 3*NCVM*NCVM + 6*NCVM ), ARTOL
      DOUBLE PRECISION RESID( N2 ), V( N2, NCVM ), WVEC( N2 ),
     1                 W2( N2 ), SVFDC( NCFM, NR, NDRVM+1, NDCS ),
     2                 PARAM( 10 )
      DOUBLE PRECISION
     1                 QST( LH*(LH+2), 3),
     2                 VF( NPHP, NTHP, 3), F1( 2*NPHP ),
     3                 F2( 2*NPHP ), F3( 2*NPHP )
      DOUBLE PRECISION
     1                 GAUX( NTHP ), GAUW( NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
      LOGICAL SELECT( NCVM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NCE, ICLS, IH, NOIT, ITGN
      DOUBLE PRECISION CA, CE, CB1, CB2, CD, CG, CH, CI, RESTOL,
     1                 GRR, GRI, DSRSV, DIAGEL, DELTA, ONEPD, ZERO,
     2                 FACTOR, CTOL, DM, DC, CHOLD, GRROLD, DX, DY,
     3                 LOLIM, HILIM
      PARAMETER ( ICLS = 1, DIAGEL = -200.0d0, DELTA = 1.0d-6,
     1            ONEPD = 1.0d0 + DELTA, ZERO = 0.0d0 )
      LOGICAL OSBR, OLOLIM, OHILIM
C ololim and ohilim are respectively set to .true. only
C when the double precision variables dlolim and dhilim have
C been assigned values.
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      ITGN    = NR/2
C
      OLOLIM  = .FALSE.
      OHILIM  = .FALSE.
C
C Unload initial values of PARAM ...
C
      CA      = PARAM(  1 )
      CB1     = PARAM(  2 )
      CB2     = PARAM(  3 )
      CD      = PARAM(  4 )
      CE      = PARAM(  5 )
      CG      = PARAM(  6 )
      CH      = PARAM(  7 )
      CI      = PARAM(  8 )
      FACTOR  = PARAM(  9 )
      CTOL    = PARAM( 10 )
C
      IF ( FACTOR.LE.ONEPD ) THEN
        PRINT *,' Subroutine LOCRNF.'
        PRINT *,' FACTOR = ', FACTOR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( CTOL.LT.ZERO ) THEN
        PRINT *,' Subroutine LOCRNF.'
        PRINT *,' CTOL = ', CTOL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LULOG.NE.0 .AND. LULOG.NE.144 .AND. LULOG.NE.145 .AND.
     1                 LULOG.NE.44 .AND. LULOG.NE.45 ) THEN
        PRINT *,' Subroutine LOCRNF.'
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
      PARAM(  1 ) = ZERO
      PARAM(  2 ) = CB1
      PARAM(  3 ) = CB2
      PARAM(  4 ) = CD
      PARAM(  5 ) = ZERO
      PARAM(  6 ) = CG
      PARAM(  7 ) = CH
      PARAM(  8 ) = CI
C
C Matrix must be zero'd on every iteration
C
C Form the linear stability matrix
C
      CALL AVMLTA( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1             NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR,
     2             NTHP, NPHP, MVAL, LH, GAUX, GAUW, PA, DPA,
     3             F1, F2, F3, VF, QST, PARAM, ICLS )
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
      IF ( LULOG.EQ.145 .OR. LULOG.EQ.45 ) THEN
        WRITE ( LULOG, 17 ) 
        WRITE ( LULOG, 16 ) NOIT, NCE
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
        PRINT *,' Subroutine LOCRNF.'
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
        IF ( LULOG.EQ.145 .OR. LULOG.EQ.45 ) WRITE ( LULOG, 19 ) DSRSV
      ENDIF
 19   FORMAT('Real shift changed to ',1pd15.7)
C
C Ok - now we know GRR and GRI
C See if we have converged :-
C
      IF ( DABS( GRR ).LT.CTOL ) THEN
        IERR = NOIT
        GOTO 60
      ENDIF
C
C We have not converged - let's try to modify CH :-
C
      IF ( GRR.GT.ZERO ) THEN
        IF ( OHILIM .AND. CH.LT.HILIM ) THEN
          HILIM = CH
        ENDIF
        IF ( .NOT. OHILIM ) THEN
          OHILIM = .TRUE.
          HILIM  = CH
        ENDIF
      ENDIF
C
      IF ( GRR.LT.ZERO ) THEN
        IF ( OLOLIM .AND. CH.GT.LOLIM ) THEN
          LOLIM = CH
        ENDIF
        IF ( .NOT. OLOLIM ) THEN
          OLOLIM = .TRUE.
          LOLIM  = CH
        ENDIF
      ENDIF
C
      IF ( OLOLIM .AND. OHILIM .AND. HILIM.LT.LOLIM ) THEN
        IERR = -3
        GOTO 60
      ENDIF
C
      IF ( (LULOG.EQ.145 .OR. LULOG.EQ.45) .AND. OLOLIM .AND.
     1          OHILIM ) WRITE ( LULOG, 21 ) NOIT, LOLIM, HILIM
 21   FORMAT('Iter. ',I4,' R_c within [',1pd16.7,',',1pd16.7,']')
C
C Special case for NOIT = 1 ...
C
      IF ( NOIT.EQ.1 ) THEN
        IF ( GRR.LT.ZERO ) THEN
          CHOLD = CH
          CH    = CHOLD*FACTOR
        ELSE
          CHOLD = CH
          CH    = CHOLD/FACTOR
        ENDIF
        GRROLD  = GRR
        GOTO 50
      ENDIF
C
C Ok - NOIT.gt.1 and so we have values for CHOLD and GRROLD
C
      DX    = CH - CHOLD
      DY    = GRR - GRROLD
      DM    = DY/DX
C     .
C     . Gradient is zero - cannot do Newton method
C     .
      IF ( DM.LT.DELTA ) THEN
C
C If have bounds ...
C
       IF ( OHILIM .AND. OLOLIM ) THEN
         CHOLD  = CH
         GRROLD = GRR
         CH = 0.5d0*( LOLIM + HILIM )
         GOTO 50
       ELSE
C
C and if there are no bounds ...
C
        IERR = -2
        GOTO 60
C
       ENDIF
C
      ENDIF
      DC     = GRR - DM*CH
      CHOLD  = CH
      GRROLD = GRR
      CH     = (-1.0d0)*DC/DM
C
      GOTO 50
C
 60   CONTINUE
      PARAM(  1 ) = CA
      PARAM(  2 ) = CB1
      PARAM(  3 ) = CB2
      PARAM(  4 ) = CD
      PARAM(  5 ) = CE
      PARAM(  6 ) = CG
      PARAM(  7 ) = CH
      PARAM(  8 ) = CI
      PARAM(  9 ) = GRR
      PARAM( 10 ) = GRI
C
      RETURN
      END
C*********************************************************************
