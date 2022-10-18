C*********************************************************************
C subroutine Drifting Frame Convection Newton Raphson Solve **********
C            -        -     -          -      -       -     **********
C Steve Gibbons Mon Mar 13 16:12:23 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C We wish to solve a non-linear time-dependent convection problem    C
C with the only time dependence being one of a constant drift in     C
C longitude.                                                         C
C                                                                    C
C Our equations are                                                  C
C                                                                    C
C  c_a d \Theta/ dt = CD \nabla^2 \Theta                             C
C                     - CC v . \nabla ( \Theta )                     C
C                     + v . ( CB1 r + CB2 r^{-2} , 0 , 0 )           C
C                                                                    C
C  c_e \curl dv/dt  = CI \nabla^2 \curl v                            C
C                     - CG \curl ( k \times v )                      C
C                     + CH \curl ( \Theta {\bm r } )                 C
C                     + CF \curl ( v \times \curl v )                C
C                                                                    C
C  We wish to minimise the vector                                    C
C                                                                    C
C                     CD \nabla^2 \Theta                             C
C                     - CC v . \nabla ( \Theta )                     C
C                     + v . ( CB1 r + CB2 r^{-2} , 0 , 0 )           C
C                     - c_a d \Theta/ dt                             C
C                                                                    C
C                     CI \nabla^2 \curl v                            C
C                     - CG \curl ( k \times v )                      C
C                     + CH \curl ( \Theta {\bm r } )                 C
C                     + CF \curl ( v \times \curl v )                C
C                     - c_e \curl dv/dt                              C
C                                                                    C
C and we use a Newton Raphson-type iteration to find the finite      C
C amplitude solution and the drift rate, CDRIFT.                     C
C                                                                    C
C The success will depend highly upon the initial vector, VEC,       C
C supplied. VEC is supplied so that the harmonic number ITCDC        C
C is zero at grid node NR/2 (by applying an appropriate shift        C
C in longitude).                                                     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Int. parameter array corresponding to VI.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR      See INDFUN for details       C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C     MXATT     : Maximum number of iterations allowed.              C
C                                                                    C
C     IERR      : Error flag on return.                              C
C                 If IERR greater than zero, all is well and         C
C                 IERR contains the number of iterations required    C
C                 for the non-linear solution to converge.           C
C                                                                    C
C                 if IERR = -1, more than the maximum number of      C
C                 iterations were needed to converge.                C
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
C     LULOG     : Logical unit number of log file.                   C
C                 LULOG = 0 if no log file is to be written.         C
C                 Otherwise, vector norms and drift rate updates     C
C                 are written out to LULOG.                          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Initial solution vector. Dim (N2)                  C
C                                                                    C
C     PARAM     : Array dim ( 13 ) containing these param.s          C
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
C            PARAM( 11 ) = CDRIFT (guess on input, value on output)  C
C            PARAM( 12 ) = CTOL (Convergence parameter for RHS norm) C
C            PARAM( 13 ) = value to which harmonic ITCDC at grid     C
C                          node NR/2 is to be set. Preferably small. C
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
C     FDCM      : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, 1 ).                        C
C                   Array is generated by the routine fdcmbd         C
C                 See documentation for FDCMBD for details.          C
C       MUST be calculated with:                                     C
C                                                                    C
C                     NDRVM = 1                                      C
C                     NLMN  = 2                                      C
C                     NRMN  = NR - 1                                 C
C                     NLMC  = 2                                      C
C                     NRMC  = NR - 1                                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DFCNRS( VEC, INARR, NR, PARAM, MXATT, IERR, A, N1,
     1  N2, KL, MHT, MHL, MHM, MHP, MHTR, NBN, NCFM, NDRVM, NDCS,
     2  NTHP, NPHP, MMAX, LH, MHIBC, MHOBC, IPIV, LULOG, SVFDC, XARR,
     3  GAUX, GAUW, PA, DPA, FDCM, ITCDC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), NR, MXATT, IERR, N1, N2, KL, MHT( * ),
     1        MHL( * ), MHM( * ), MHP( * ), MHTR( * ), NBN, NCFM,
     2        NDRVM, NDCS, NTHP, NPHP, MMAX, LH, MHIBC( * ),
     3        MHOBC( * ), IPIV( * ), ITCDC, LULOG
      DOUBLE PRECISION VEC( * ), PARAM( * ), A( N1, N2 ), XARR( NR ),
     1        SVFDC( NCFM, NR, NDRVM+1, NDCS ), GAUX( NTHP ),
     2        GAUW( NTHP ), PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3        DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
      DOUBLE PRECISION FDCM( NCFM, NR, 1 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER INDFUN, IRZ, IELZ, NH, ICLS, NOIT, I, IOP,
     1        ILN, IRN, IHD, NDRVS, ICMP, ICMPO, NRMAX, NVI1, NVI2,
     2        LHMAX, NPHMAX, NTHMAX, M0, MAXNVI, NHMAX, ISVMAX,
     3        NTS, IH, IT10, IT11C, IT11S, IS, ILNT, IRNT, IMF
      DOUBLE PRECISION ZERO, CA, CB1, CB2, CC, CD, CE, CF, CG,
     1        CH, CI, CDRIFT, PARS( 8 ), CTOL, FAC, DLOW, DONE,
     2        SOLN, RHSN, DNRM2
      PARAMETER ( ZERO = 0.0d0, IOP = 0, NRMAX = 200,
     1            NPHMAX = 128, NTHMAX = 64, LHMAX = 62,
     2            MAXNVI = 20000, NHMAX = 100, DONE = 1.0d0, 
     3            ISVMAX = NRMAX*NHMAX, DLOW = 1.0d-8 )
      PARAMETER ( IMF = 1 )
      INTEGER KKA1( MAXNVI ), KKB1( MAXNVI ), KKG1( MAXNVI ),
     1        KKA2( MAXNVI ), KKB2( MAXNVI ), KKG2( MAXNVI ),
     2        IPIVH( 4 )
      DOUBLE PRECISION ZCFA( NRMAX ), ZCFB( NRMAX ), ZCFC( NRMAX ),
     1                 F1( 2*NPHMAX ), F2( 2*NPHMAX ), 
     2                 F3( 2*NPHMAX ), SHC( LHMAX*(LHMAX + 2) ),
     3                 DSHC( LHMAX*(LHMAX + 2) ), CVI1( MAXNVI ),
     4                 CVI2( MAXNVI )
      CHARACTER *(4) TVHI1( MAXNVI ), TVHI2( MAXNVI )
      DOUBLE PRECISION QST( LHMAX*(LHMAX + 2), 3 ),
     1                 RQST1( LHMAX*(LHMAX + 2), 3, NRMAX ),
     2                 RQST2( LHMAX*(LHMAX + 2), 3, NRMAX ),
     3                 RQST3( LHMAX*(LHMAX + 2), 3, NRMAX ),
     4                 RQSTA( LHMAX*(LHMAX + 2), 3, NRMAX )
      DOUBLE PRECISION VF1( NPHMAX, NTHMAX, 3),
     1                 VF2( NPHMAX, NTHMAX, 3),
     2                 VF3( NPHMAX, NTHMAX, 3),
     3                 SF( NPHMAX, NTHMAX ), RHS( ISVMAX )
      DOUBLE PRECISION V0( ISVMAX ), V1( ISVMAX ), V2( ISVMAX ),
     1                 V3( ISVMAX ), V4( ISVMAX ), UM( ISVMAX, 4 ),
     2                 VM( ISVMAX, 4 ), HMAT( 4, 4), VTY( 4, 1 ),
     2                 HVTY( 4, 1 ), UHVTY( ISVMAX, 1 )
      CHARACTER *(3) CHVMFF
      LOGICAL OT10, OT11C, OT11S, OTCDC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      ILN   = 2
      IRN   = NR - 1
      ILNT  = 3
      IRNT  = NR - 2
C
      OT10  = .FALSE.
      OT11C = .FALSE.
      OT11S = .FALSE.
      OTCDC = .TRUE.
C
      NH    = INARR( 3 )
C
      IF ( NR.GT.NRMAX .OR. NPHP.GT.NPHMAX .OR. LH.GT.LHMAX .OR.
     1     NTHP.GT.NTHMAX .OR. NH.GT.NHMAX .OR. N2.GT.ISVMAX ) THEN
        PRINT *,' Subroutine DFCNRS'
        PRINT *,' NR   = ', NR,  ' NRMAX  = ', NRMAX
        PRINT *,' NH   = ', NH,  ' NHMAX  = ', NHMAX
        PRINT *,' LH   = ', LH,  ' LHMAX  = ', LHMAX
        PRINT *,' N2   = ', N2,  ' ISVMAX = ', ISVMAX
        PRINT *,' NPHP = ', NPHP,' NPHMAX = ', NPHMAX
        PRINT *,' NTHP = ', NTHP,' NTHMAX = ', NTHMAX
        PRINT *,' Recompile routine with higher dimensions.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C We won't bother checking the bounds and dimensions of A
C as we are calling AVMLTA which does this all for us.
C
      CA     = PARAM(  1 )
      CB1    = PARAM(  2 )
      CB2    = PARAM(  3 )
      CC     = PARAM(  4 )
      CD     = PARAM(  5 )
      CE     = PARAM(  6 )
      CF     = PARAM(  7 )
      CG     = PARAM(  8 )
      CH     = PARAM(  9 )
      CI     = PARAM( 10 )
      CDRIFT = PARAM( 11 )
      CTOL   = PARAM( 12 )
C
      IRZ  = NR/2
      IELZ = INDFUN( IRZ, ITCDC, INARR )
C
      VEC( IELZ ) = PARAM( 13 )
C
C Calculate NTS. Now this is the number of additional
C vectors we need to use to solve the linear system.
C NTS must always be atleast 1 for a drifting frame.
C This is because an extra column  for the matrix is
C needed to modify the drift velocity, CDRIFT.
C
      NTS = 1
C
C If our boundaries are stress free, then we may
C require up to three rows to be added to the matrix
C to fix the frame of reference.
C
      DO IH = 1, NH
        IS  = MHP( IH )
        IF (  MHIBC( IS ).EQ.6 .AND. MHOBC( IS ).EQ.6 .AND.
     1        MHT( IH ).EQ.2 .AND. MHL( IH ).EQ.1  ) THEN
C         .
          IF ( MHM( IH ).EQ.0 ) THEN
            NTS   = NTS + 1
            OT10  = .TRUE.
            IT10  = IH
          ENDIF
C         .
          IF ( MHM( IH ).EQ.1 .AND. DABS( CG ).LT.DLOW ) THEN
            NTS   = NTS + 1
            OT11C = .TRUE.
            IT11C = IH
          ENDIF
C         .
          IF ( MHM( IH ).EQ.-1 .AND. DABS( CG ).LT.DLOW ) THEN
            NTS   = NTS + 1
            OT11S = .TRUE.
            IT11S = IH
          ENDIF
C         .
        ENDIF
      ENDDO
C
C We now calculate the vector interactions which
C are used to form the matrix required for Newton-Raphson
C
      NVI1 = 0
      NVI2 = 0
C
C First: Scalar product ( v . \nabla \Theta )
C
      CALL VSPCC( NVI1, MAXNVI, KKA1, KKB1, KKG1, NH, MHT, MHL, MHM,
     1            NH, MHT, MHL, MHM, NH, MHTR, MHL, MHM, LH, NTHP,
     2            NPHP, MMAX, TVHI1, CVI1, GAUX, GAUW, PA, DPA,
     3            F1, VF1, VF2, SF, SHC )
C
C Second: Inertial term ( - \curl ( v \times \curl v )   )
C
      CALL VCCPCC( NVI2, MAXNVI, KKA2, KKB2, KKG2, NH, MHT, MHL, MHM,
     1             NH, MHT, MHL, MHM, NH, MHTR, MHL, MHM, LH, NTHP,
     2             NPHP, MMAX, TVHI2, CVI2, GAUX, GAUW, PA, DPA,
     3             F1, F2, F3, VF1, VF2, VF3, QST )
C
C We will now begin the iteration of the Newton
C Raphson process.
C
      NDRVS   = 4
      NOIT    = 0
 50   CONTINUE
      NOIT = NOIT + 1
C
      IF ( NOIT.GT.MXATT ) THEN
        IERR = -1
        GOTO 151
      ENDIF
C
C Evaluate the right hand side of system
C First zero RHS vector
C
      CALL VECOP( RHS, ZERO, N2, IOP )
C
C Take derivatives of current vector, VEC.
C
      IHD     = 4
      CALL CASVDR ( VEC, ILN, IRN, NBN, IHD, NCFM, NR, NDRVS,
     1              NDRVM, INARR, NDCS, MHP, SVFDC, V0, V1,
     2              V2, V3, V4 )
C
C Add diffusive parts to the RHS
C
      FAC  = CI
      CALL SSVLC( NR, V0, V1, V2, V3, V4, INARR, MHT,
     1            MHL, MHM, RHS, INARR, MHTR, MHL, MHM,
     2            FAC, ILN, IRN, XARR )
C
      FAC  = CD
      ICMP = 3
      CALL SSVLP( NR, V0, V1, V2, INARR, MHT, MHL, MHM,
     1            ICMP, RHS, INARR, MHTR, MHL, MHM, FAC,
     2            ILN, IRN, XARR )
C
C Form drifting frame derivatives and store in
C V4 (we can do this as V4 is no longer needed in
C the current iteration).
C
      CALL VECOP( V4, ZERO, N2, IOP )
      FAC = DONE
      CALL DFTDA( INARR, MHT, MHL, MHM, MHTR, V0, V1, V2,
     1            V4, FAC, CA, CE, ZERO, XARR )
C
C Subtract drifting frame time derivative
C from RHS vector.
C
      DO I = 1, N2
        RHS( I ) = RHS( I ) - CDRIFT*V4( I )
      ENDDO
C     .
C     . Add on heat source terms
C     .
      FAC   = DONE
      CALL SSVHST( NR, V0, INARR, MHT, MHL, MHM, RHS,
     1             INARR, MHTR, MHL, MHM, FAC, ILN, IRN,
     2             XARR, CB1, CB2 )
C     .
C     . add on Buoyancy terms
C     .
      FAC   = CH
      ICMP  = 3
      ICMPO = 2
      CALL SSVTA( NR, V0, INARR, MHT, MHL, MHM, ICMP,
     1            RHS, INARR, MHTR, MHL, MHM, ICMPO,
     2            FAC, ILNT, IRNT )
C     .
      CALL VECOP( ZCFA, ZERO, NR, IOP )
      CALL VECOP( ZCFB, ZERO, NR, IOP )
      CALL VECOP( ZCFC, ZERO, NR, IOP )
C     .
C     . Now begin non-linear terms
C     . First put velocity into RQST1 array
C     .
      CHVMFF = 'VEL'
      CALL SDRQST( NR, LH, V0, V1, ILN, IRN, INARR,
     1             RQST1, XARR, MHT, MHL, MHM, CHVMFF )
C     .
C     . add on v . Grad( Theta ) terms to RHS
C     .
      FAC   = -1.0d0*CC
      CALL SDVGTA( NR, LH, MMAX, MHT, MHL, MHM, V0, V1, INARR,
     1             MHTR, MHL, MHM, RHS, INARR, RQST1, VF1, VF2,
     2             SF, F1, F2, F3, SHC, DSHC, GAUX, GAUW, PA, DPA,
     3             NTHP, NPHP, FAC, ILN, IRN, XARR, ZCFA)
C     .
C     . RQST1 contains the velocity, v, so by calling
C     . RQSTCF we can put (k x v) into RQST2.
C     .
      CALL RQSTCF( NR, LH, MMAX, ILN, IRN, NTHP, NPHP,
     1             GAUX, GAUW, PA, DPA, RQST1, ZCFA, RQST2,
     2             ZCFB, VF1, F1, F2, F3 )
C     .
C     . Taking curl of RQST2 and subtract (CG*) this
C     . amount from RQSTA
C     .
      ICLS  = 1
      FAC   = (-1.0d0)*CG
      M0    = 1
      CALL RQSTCA( LH, NR, M0, MMAX, NBN, NCFM, NDRVS, ILN, IRN,
     1             ILN, IRN, RQST2, RQSTA, XARR, FDCM, ICLS, ZERO,
     2             FAC )
C     .
C     . Take curl of velocity and store curl in RQST2
C     .
      ICLS  = 1
      FAC   = 1.0d0
      CALL RQSTCA( LH, NR, M0, MMAX, NBN, NCFM, NDRVS, ILN, IRN, ILN,
     1             IRN, RQST1, RQST2, XARR, FDCM, ICLS, ZERO, FAC )
C     .
C     . Evaluate [ v x curl v ] in RQST3
C     .
      CALL RQSTCP( NR, LH, MMAX, ILN, IRN, NTHP, NPHP,
     1             GAUX, GAUW, PA, DPA, RQST1, ZCFA, RQST2,
     2             ZCFB, RQST3, ZCFC, VF1, VF2, VF3, F1, F2, F3 )
C     .
C     . Add CF* curl of RQST3 to RQSTA
C     .
      ICLS  = 0
      FAC   = CF
      CALL RQSTCA( LH, NR, M0, MMAX, NBN, NCFM, NDRVS, ILN, IRN,
     1             ILN, IRN, RQST3, RQSTA, XARR, FDCM, ICLS, DONE, FAC )
C     .
C     . RQSTA now contains (   - CG curl (k x v)
C     .                      + CF curl[ v x (curl v) ]    )
C     . So - add this onto the vector RHS
C     .
      FAC = 1.0d0
      CHVMFF = 'VEL'
      CALL RQSTSV( NR, LH, ILN, IRN, INARR, MHTR, MHL, MHM,
     1             CHVMFF, RQSTA, RHS, XARR, FAC )
C
C Right hand side now complete. Calculate norm. (imf = 1)
C
      RHSN = DNRM2( N2, RHS, IMF )
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 190 ) RHSN
C
 190  FORMAT('dfcnrs: RHS norm = ',1pd16.7)
C
C Form the linear parts of matrix.
C If ICLS = 1, this also zeros the matrix
C
      PARS(  1 ) = ZERO
      PARS(  2 ) = (-1.0d0)*CB1
      PARS(  3 ) = (-1.0d0)*CB2
      PARS(  4 ) = (-1.0d0)*CD
      PARS(  5 ) = ZERO
      PARS(  6 ) = (-1.0d0)*CG
      PARS(  7 ) = (-1.0d0)*CH
      PARS(  8 ) = (-1.0d0)*CI
C
      ICLS = 1
      CALL AVMLTA( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1             NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR,
     2             NTHP, NPHP, MMAX, LH, GAUX, GAUW, PA, DPA,
     3             F1, F2, F3, VF1, QST, PARS, ICLS )
C
C Add time derivative terms to matrix
C First: temperature 
C
      FAC   = CA*CDRIFT
      ICMP  = 3
      CALL AMDFTD( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1             MHTR, MHL, MHM, ICMP, FAC, NBN, NDRVS,
     2             NDRVM, ILN, IRN, XARR, NCFM, SVFDC,
     3             A, N1, N2, IMF, KL, KL, KL, NDCS )
C
C Now: poloidal and toroidal velocities
C
      FAC = CE*CDRIFT
      ICMP  = 1
      ICMPO = 2
      CALL AMDFTD( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1             MHTR, MHL, MHM, ICMPO, FAC, NBN, NDRVS,
     2             NDRVM, ILNT, IRNT, XARR, NCFM, SVFDC,
     3             A, N1, N2, IMF, KL, KL, KL, NDCS )
C
      ICMP  = 2
      ICMPO = 1
      CALL AMDFTD( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1             MHTR, MHL, MHM, ICMPO, FAC, NBN, NDRVS,
     2             NDRVM, ILN, IRN, XARR, NCFM, SVFDC,
     3             A, N1, N2, IMF, KL, KL, KL, NDCS )
C
C Now add non-linear terms to matrix.
C We have already calculated the coefficients for this
C using VSPCC and VCCPCC ...
C
C Firstly: curl ( v0 x curl v_new ) terms
C
      FAC   = CF
      CALL RV0CVA( NR, N1, N2, KL, KL, KL, IMF, ILN, IRN,
     1    ILNT, IRNT, INARR, MHT, MHL, MHM, MHP, MHTR, MHL, MHM,
     2    INARR, MHT, MHL, MHM, MHP, NBN, NDCS, NDRVS, NDRVM, NCFM,
     3    NBN, NDCS, NDRVS, NDRVM, NCFM, A, FAC, XARR, VEC, SVFDC,
     4    SVFDC, NVI2, KKA2, KKB2, KKG2, TVHI2, CVI2 )
C
C Now: curl ( v_new x curl v0 ) terms
C
      CALL RVCV0A( NR, N1, N2, KL, KL, KL, IMF, ILN, IRN,
     1    ILNT, IRNT, INARR, MHT, MHL, MHM, MHP, MHTR, MHL, MHM,
     2    INARR, MHT, MHL, MHM, MHP, NBN, NDCS, NDRVS, NDRVM, NCFM,
     3    NBN, NDCS, NDRVS, NDRVM, NCFM, A, FAC, XARR, VEC, SVFDC,
     4    SVFDC, NVI2, KKA2, KKB2, KKG2, TVHI2, CVI2 )
C
C Now: v0 . nabla ( THETA_new ) terms
C
      FAC   = CC
      CALL RV0GTA( NR, N1, N2, KL, KL, KL, IMF, ILN, IRN, INARR,
     1    MHT, MHL, MHM, MHP, MHTR, MHL, MHM, INARR, MHT, MHL, MHM,
     2    MHP, NBN, NDCS, NDRVS, NDRVM, NCFM, NBN, NDCS, NDRVS,
     3    NDRVM, NCFM, A, FAC, XARR, VEC, SVFDC, SVFDC,
     4    NVI1, KKA1, KKB1, KKG1, TVHI1, CVI1 )
C
C Now: v_new . nabla ( THETA_0 ) terms
C
      CALL RVGT0A( NR, N1, N2, KL, KL, KL, IMF, ILN, IRN, INARR,
     1    MHT, MHL, MHM, MHP, MHTR, MHL, MHM, INARR, MHT, MHL, MHM,
     2    MHP, NBN, NDCS, NDRVS, NDRVM, NCFM, NBN, NDCS, NDRVS,
     3    NDRVM, NCFM, A, FAC, XARR, VEC, SVFDC, SVFDC,
     4    NVI1, KKA1, KKB1, KKG1, TVHI1, CVI1 )
C
C Now prevent singularity of matrix through checking
C of boundary points ...
C
      CALL AMSDEA( A, N1, N2, KL, KL, KL, IMF, INARR,
     1             MHP, MHIBC, 'Inner', DONE, NDCS )
      CALL AMSDEA( A, N1, N2, KL, KL, KL, IMF, INARR,
     1             MHP, MHOBC, 'Outer', DONE, NDCS )
C
C Treat matrix for free rotations and additional
C columns with appropriate adjustments to RHS and
C formation of the correct Woodbury formula matrices
C 
      CALL NRCWMF( N1, N2, KL, KL, KL, NTS, ITCDC, IT10, IT11C,
     1             IT11S, NR, INARR, OTCDC, OT10, OT11C, OT11S,
     2             A, RHS, UM, VM, XARR, V4 )
C
C We have now fully prepared the matrix for Woodbury
C formula solution with the routine BMWDFS. We need to
C LU decompose so call ILUDF as IMF as this is 1.
C
      CALL BMWDFS( N1, N2, NTS, KL, KL, KL, IPIV, IPIVH, A,
     1             RHS, UM, VM, HMAT, VTY, HVTY, UHVTY, IMF )
C
      CALL ASVCPL( RHS, NR, NDCS, INARR, MHP, MHIBC, MHOBC,
     1             NCFM, NDRVS, NDRVM, NBN, SVFDC )
C
C Solution must have been successful. So extract value
C of new CDRIFT
C
      CDRIFT = CDRIFT + RHS( IELZ )
      RHS( IELZ ) = 0.0d0
C
C Copy RHS back into VEC for next iteration (imf = 1)
C
      RHSN = DNRM2( N2, RHS, IMF )
      DO I = 1, N2
        VEC( I ) = VEC( I ) + RHS( I )
      ENDDO
      SOLN = DNRM2( N2, VEC, IMF )
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 191 ) NOIT, CDRIFT
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 192 ) RHSN, SOLN
C
C Judge whether convergence has been achieved
C
      IF ( RHSN.LT.CTOL ) THEN
        IERR = NOIT
        GOTO 151
      ENDIF
C
 191  FORMAT('dfcnrs: Iteration ',I3,' cdrift = ',1pd16.7)
 192  FORMAT('Sol. norms: change: ',1pd16.7,' total: ',1pd16.7)
C
C Return to the beginning of loop for next iteration
C
      GOTO 50
C
 151  CONTINUE
      PARAM( 11 ) = CDRIFT
      RETURN
      END
C*********************************************************************
