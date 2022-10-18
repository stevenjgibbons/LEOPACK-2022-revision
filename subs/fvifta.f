C*********************************************************************
C subroutine Fixed Velocity Instability Forcing Term Add *************
C            -     -        -           -       -    -   *************
C Steve Gibbons Mon May 15 19:09:56 BST 2000                         *
C____________________________________________________________________C
C                                                                    C
C Adds to a vector RHS the following terms in the MHD equations:     C
C                                                                    C
C F_Theta =  u . ( CB1 r + CB2 r^{-2} , 0 , 0 ) :                    C
C            - CC u_0 . Grad ( Theta )          :                    C
C                                                                    C
C F_vort. =     CH curl ( Theta )               :                    C
C             - CG curl ( K x v )               :                    C
C             - CF curl ( v_0. Grad) v          :                    C
C             - CF curl ( v. Grad) v_0          :                    C
C                                                                    C
C Here, the arrays V0RVF and CV0RVF store v_0 and curl v_0 in        C
C RVF format.                                                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARR     : Format array for solution vectors. Dim (3).        C
C                                                                    C
C     MHT       : Dim ( * ). Harmonic type.                          C
C     MHL       : Dim ( * ). Degree, l.                              C
C     MHM       : Dim ( * ). m for cos m phi, -m for sin m phi.      C
C                                                                    C
C     MHTR      : Dim ( * ). MHT after calling CINDSW.               C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NCFM      : Leading coefficient of FDCM array.                 C
C     NBN       : Number of bounding nodes. See ASVDR.               C
C     MMAX      : Maximum value of m.                                C
C     M0        : Minimum value of non-zero m such that              C
C                  MOD( m, M0 ) = 0 for all m.                       C
C     NPHP      : Number of points in phi.                           C
C     NTHP      : Number of points in theta.                         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     CB1       : Coefficient of v_r.r term in heat equation.        C
C     CB2       : Coefficient of v_r/(rr) term in heat equation.     C
C     CC        : Coefficient of v.Grad(Theta) term in heat eqn.     C
C     CH        : Coefficient of buoyancy term in vorticity eqn.     C
C     CG        : Coefficient of Coriolis term in vorticity eqn.     C
C     CF        : Coefficient of Inertial term in vorticity eqn.     C
C                                                                    C
C     V0        : Vector at time step i, zero^th derivatives.        C
C     V1        : Vector at time step i, first   derivatives.        C
C  (Pre-calculate these with CASVDR)                                 C
C                                                                    C
C     V0RVF     : Fixed velocity, v_0. Dim( NR, NPHP, NTHP, 3 )      C
C     CV0RVF    : Curl of v_0. Dim( NR, NPHP, NTHP, 3 )              C
C                                                                    C
C     RVF1      : Work arr. dim. ( NR, NPHP, NTHP, 3 )               C
C     RVF2      : Work arr. dim. ( NR, NPHP, NTHP, 3 )               C
C     RVF3      : Work arr. dim. ( NR, NPHP, NTHP, 3 )               C
C     SF        : Work arr. dim. ( NPHP, NTHP )                      C
C     VF        : Work arr. dim. ( NPHP, NTHP, 3 )                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Gauss weights computed bu GAUWTS. Dim ( NTHP )     C
C     SHC       : Work array. Dim ( LH*(LH+2) ).                     C
C     DSHC      : Work array. Dim ( LH*(LH+2) ).                     C
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
C     RHS       : Right hand side vector.                            C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     RQST1     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C     RQST2     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C     RQST3     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C     RQSTA     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FVIFTA( NR, INARR, MHT, MHL, MHM, MHTR, LH, NCFM,
     1                 NBN, MMAX, M0, NTHP, NPHP, CB1, CB2, CC, CH,
     2                 CG, CF, V0, V1, RVF1, RVF2, RVF3, SF, VF, PA,
     3                 DPA, GAUX, GAUW, SHC, DSHC, FDCM, RHS, XARR,
     4                 RQST1, RQST2, RQST3, RQSTA, V0RVF, CV0RVF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * ), MHTR( * ),
     1        LH, NCFM, NBN, MMAX, M0, NTHP, NPHP, NR
      DOUBLE PRECISION FDCM( NCFM, NR, 1 ), SF( NPHP, NTHP ),
     1                 RVF1( NR, NPHP, NTHP, 3 ), GAUX( NTHP ),
     2                 RVF2( NR, NPHP, NTHP, 3 ), GAUW( NTHP ),
     3                 RVF3( NR, NPHP, NTHP, 3 ), VF( NPHP, NTHP, 3 )
      DOUBLE PRECISION V0RVF( NR, NPHP, NTHP, 3 ),
     1                 CV0RVF( NR, NPHP, NTHP, 3 )
      DOUBLE PRECISION SHC( LH*(LH+2) ), DSHC( LH*(LH+2) ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
      DOUBLE PRECISION CB1, CB2, CC, CH, CG, CF, V0( * ),
     1                 V1( * ), RHS( * ), XARR( * ),
     2                 RQSTA( LH*(LH+2), 3, NR )
      DOUBLE PRECISION RQST1( LH*(LH+2), 3, NR ),
     1                 RQST2( LH*(LH+2), 3, NR ),
     2                 RQST3( LH*(LH+2), 3, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER ILNR, IRNR, ILNT, IRNT, ICMPI, ICMPO,
     1        NPHMAX, NRMAX, NDRVS, ICLS, IOP, IRES
      PARAMETER ( NPHMAX = 128, NRMAX = 300 )
      DOUBLE PRECISION FAC, A, B, FTF1( 2*NPHMAX ), FTF2( 2*NPHMAX ),
     1                 FTF3( 2*NPHMAX ), ZERO, ZCFA( NRMAX ), DLOW
      PARAMETER ( IOP = 0, ZERO = 0.0d0, IRES = 2, DLOW = 1.0d-10 )
      CHARACTER *(3) CHVMFF
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      CALL VECOP( ZCFA, ZERO, NRMAX, IOP )
C
      IF ( NR.NE.INARR( 2 ) ) THEN
        PRINT *,' Subroutine FVIFTA.'
        PRINT *,' NR = ', NR,', INARR(2) = ', INARR( 2 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check our work arrays are sufficiently large.
C     .
      IF ( NPHP.GT.NPHMAX .OR. NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine FVIFTA.'
        PRINT *,' NR   = ', NR,  ' NRMAX  = ', NRMAX
        PRINT *,' NPHP = ', NPHP,' NPHMAX = ', NPHMAX
        PRINT *,' Redimension and recompile routine.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . First add on heat source terms
C     .
      ILNR  = 2
      IRNR  = NR - 1
      FAC   = 1.0d0
      CALL SSVHST( NR, V0, INARR, MHT, MHL, MHM, RHS,
     1             INARR, MHTR, MHL, MHM, FAC, ILNR, IRNR,
     2             XARR, CB1, CB2 )
C     .
C     . add on Buoyancy terms
C     .
      ILNT  = 3
      IRNT  = NR - 2
      FAC   = CH
      ICMPI = 3
      ICMPO = 2
      CALL SSVTA( NR, V0, INARR, MHT, MHL, MHM, ICMPI,
     1            RHS, INARR, MHTR, MHL, MHM, ICMPO,
     2            FAC, ILNT, IRNT )
C     .
C     . Now begin non-linear terms
C     . First put velocity into RQST1 array
C     . (We will need it in this form in order to
C     . take the curl later).
C     .
      ILNR  = 2
      IRNR  = NR - 1
      CHVMFF = 'VEL'
      CALL SDRQST( NR, LH, V0, V1, ILNR, IRNR, INARR,
     1             RQST1, XARR, MHT, MHL, MHM, CHVMFF )
C     .
C     . Now store the velocity in RVF format - do this
C     . directly from the solution vector (SDVRVF).
C     .
      CALL SDVRVF( NTHP, NPHP, LH, IRES, ILNR, IRNR, NR,
     1             INARR, MHT, MHL, MHM, MMAX, CHVMFF, V0, V1,
     2             GAUX, PA, DPA, RVF1, FTF1, FTF2, FTF3, XARR )
C     .
C     . V0RVF contains the fixed velocity.
C     . add on v_0 . Grad( Theta ) terms to RHS
C     .
      FAC   = -1.0d0*CC
      CALL RVVGTA( NR, LH, MMAX, MHT, MHL, MHM, V0, V1,
     1             INARR, MHTR, MHL, MHM, RHS, INARR, V0RVF, VF,
     2             SF, FTF1, FTF2, FTF3, SHC, DSHC, GAUX, GAUW,
     3             PA, DPA, NTHP, NPHP, FAC, ILNR, IRNR, XARR )
C     .
C     . Now deal with the Coriolis force.
C     . RVF1 contains v, so call RVFCF2 to put (k x v)
C     . into RVF2.
C     .
      CALL RVFCF2( NTHP, NPHP, NR, ILNR, IRNR, RVF1,
     1                    RVF2, GAUX )
C     .
C     . RVF2 now contains (k x v) so we put this
C     . back into RQST format. (the array RQST2).
C     .
      CALL RVFRQS( RQST2, RVF2, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1        FTF3, ZCFA, LH, NTHP, NPHP, MMAX, NR, ILNR, IRNR )
C     .
C     . Taking curl of RQST2 and subtract (CG*) this
C     . amount from RQSTA
C     .
      NDRVS = 1
      ICLS  = 1
      A     = 0.0d0
      B     = (-1.0d0)*CG
      CALL RQSTCA( LH, NR, M0, MMAX, NBN, NCFM, NDRVS, ILNR, IRNR,
     1             ILNR, IRNR, RQST2, RQSTA, XARR, FDCM, ICLS, A, B)
C     .
C     . Take curl of velocity and store curl in RQST2
C     .
      NDRVS = 1
      ICLS  = 1
      A     = 0.0d0     
      B     = 1.0d0     
      CALL RQSTCA( LH, NR, M0, MMAX, NBN, NCFM, NDRVS, ILNR,
     1     IRNR, ILNR, IRNR, RQST1, RQST2, XARR, FDCM, ICLS, A, B )
C     .
C     . Evaluate (curl v) in space using RQSRVF (from RQST2)
C     . Store curl v in RVF2: (k x v) is no longer needed.
C     .
      CALL VECOP( ZCFA, ZERO, NR, IOP )
      CALL RQSRVF( RQST2, RVF2, GAUX, PA, DPA, FTF1, FTF2, FTF3,
     1              LH, NTHP, NPHP, MMAX, ZCFA, NR, ILNR, IRNR )
C     .
C     . The remaining terms are dependent upon the
C     . multiplier CF. If this number is sufficiently
C     . small, we should move to line 50 straight away.
C     .
      IF ( DABS( CF ).LT.DLOW ) GOTO 50
C     .
C     . Evaluate [ v x curl v_0 ] in RVF3
C     .
      CALL RVFCP ( RVF1, CV0RVF, RVF3, NPHP, NTHP, NR, ILNR, IRNR )
C     .
C     . RVF3 now contains [ v x curl v_0 ].
C     . Use RVFRQS to convert this function to RQST format,
C     . storing output in RQST3.
C     .
      CALL RVFRQS( RQST3, RVF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1        FTF3, ZCFA, LH, NTHP, NPHP, MMAX, NR, ILNR, IRNR )
C     .
C     . Add CF* curl of RQST3 to RQSTA
C     .
      NDRVS = 1
      ICLS  = 0
      A     = 1.0d0
      B     = CF
      CALL RQSTCA( LH, NR, M0, MMAX, NBN, NCFM, NDRVS, ILNR,
     1     IRNR, ILNR, IRNR, RQST3, RQSTA, XARR, FDCM, ICLS, A, B )
C     .
C     . Evaluate [ v_0 x curl v ] in RVF3
C     .
      CALL RVFCP ( V0RVF, RVF2, RVF3, NPHP, NTHP, NR, ILNR, IRNR )
C     .
C     . RVF3 now contains [ v_0 x curl v ].
C     . Use RVFRQS to convert this function to RQST format,
C     . storing output in RQST3.
C     .
      CALL RVFRQS( RQST3, RVF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1        FTF3, ZCFA, LH, NTHP, NPHP, MMAX, NR, ILNR, IRNR )
C     .
C     . Add CF* curl of RQST3 to RQSTA
C     .
      NDRVS = 1
      ICLS  = 0
      A     = 1.0d0
      B     = CF
      CALL RQSTCA( LH, NR, M0, MMAX, NBN, NCFM, NDRVS, ILNR,
     1     IRNR, ILNR, IRNR, RQST3, RQSTA, XARR, FDCM, ICLS, A, B )
C     .
 50   CONTINUE
C     .
C     . Now all of the non-magnetic terms have been
C     . evaluated. 
C     .
C     . RQSTA now contains (   - CG curl (k x v)
C     .                      + CF curl[ v_0 x (curl v) ]
C     .                      + CF curl[ v x (curl v_0) ]
C     . So - add this onto the vector RHS
C     .
      FAC = 1.0d0
      CHVMFF = 'VEL'
      CALL RQSTSV( NR, LH, ILNR, IRNR, INARR, MHTR, MHL, MHM,
     1             CHVMFF, RQSTA, RHS, XARR, FAC )
C     .
      RETURN
      END
C*********************************************************************
