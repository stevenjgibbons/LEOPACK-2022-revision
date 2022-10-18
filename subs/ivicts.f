C*********************************************************************
C subroutine Imposed Vector conducting Inner Core Time Step routine  *
C            -       -                 -     -    -    -             *
C Steve Gibbons Tue Apr  4 09:33:16 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Advances the solution of the dynamo problem by a single time step. C
C                                                                    C
C The equations are integrated subject to the addition of a          C
C constant vector contained in VEC0, in the same radial format as    C
C the velocity and temperature solution.                             C
C                                                                    C
C for instance U = u + FAC*u_0     etc.                              C
C                                                                    C
C VEC0 is indexed by IN0, MHT0, MHL0, MHM0 and MHP0                  C
C All harmonics indicated by these arrays must also be contained     C
C somewhere in the standard solution vector in INARRV, MHTV, MHLV,   C
C MHMV and MHPV.                                                     C
C                                                                    C
C FAC is input as the 18th element of DPARS.                         C
C____________________________________________________________________C
C                                                                    C
C  Let F_Theta, F_Vort and R_B be the forcing terms in the           C
C  heat, vorticity and induction equations respectively with         C
C                                                                    C
C  F_Theta = v ( CB1 r + CB2 r^{-2} , 0, 0 ) - CC v . Grad (Theta)   C
C                                                                    C
C  F_Vort  =   CH curl ( Theta r ) - CG curl ( k x v )               C
C             -CF curl ( v . Grad ) v + CJ curl ( B . Grad ) B       C
C                                                                    C
C  F_B     =  CM = curl ( v x B )                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Let v^{i}, Theta^{i} and B^{i} be velocity, temperature and       C
C magnetic field respectively at time step i.                        C
C                                                                    C
C  Using MCDRHF and SSVFTA, we form the forcing terms                C
C                                                                    C
C      F1 =  [CA+DELTAT.CD.CFAC Lap] Theta^{i} +                     C
C                                DELTAT. F_Theta( v,The,B^{i} ),     C
C                                                                    C
C      F2 =  [CE+DELTAT.CI.CFAC Lap] curl v^{i} +                    C
C                                DELTAT.F_Vort( v,The,B^{i} ),       C
C and                                                                C
C      F3 =  [CK+DELTAT.CL.CFAC Lap] B^{i} +                         C
C                                DELTAT.F_B( v,The,B^{i} ).          C
C                                                                    C
C A predictor ( Theta_(0)^{i+1}, v_(0)^{i+1}, B_(0)^{i+1} )          C
C is then formed by solving                                          C
C                                                                    C
C [ CA + DELTAT.CD.(CFAC-1) Lap] Theta_(0)^{i+1} = F1                C
C [ CE + DELTAT.CI.(CFAC-1) Lap] v_(0)^{i+1} = F2                    C
C [ CK + DELTAT.CL.(CFAC-1) Lap] B_(0)^{i+1} = F3                    C
C                                                                    C
C and then functional iteration is used to semi-implicitly solve     C
C the equations                                                      C
C                                                                    C
C [ CA + DELTAT.CD.(CFAC-1) Lap] Theta_(j)^{i+1} = F1                C
C [ CE + DELTAT.CI.(CFAC-1) Lap] v_(j)^{i+1} = F2                    C
C [ CK + DELTAT.CL.(CFAC-1) Lap] B_(j)^{i+1} = F3                    C
C                                                                    C
C  where                                                             C
C                                                                    C
C F1 = [CA+DELTAT.CD.CFAC Lap] Theta^{i}                             C
C      + CFAC.DELTAT.F_Theta( v,The,B^{i} )                          C
C      + CFAC.DELTAT.F_Theta( v,The,B_(j-1)^{i+1} )                  C
C                                                                    C
C F2 = [CE+DELTAT.CI.CFAC Lap] Vort^{i}                              C
C      + CFAC.DELTAT.F_Vort( v,The,B^{i} )                           C
C      + CFAC.DELTAT.F_Vort( v,The,B_(j-1)^{i+1} )                   C
C                                                                    C
C F3 = [CK+DELTAT.CL.CFAC Lap] B^{i}                                 C
C      + CFAC.DELTAT.F_B( v,The,B^{i} )                              C
C      + CFAC.DELTAT.F_B( v,The,B_(j-1)^{i+1} )                      C
C                                                                    C
C This functional iteration continues until                          C
C                                                                    C
C DABS( |V_(j)^{i+1}| - |V_{j-1}^{i+1}| ) < TOL                      C
C                                                                    C
C where V represents the full solution vector.                       C
C                                                                    C
C INARRV( 1 ) for VECV = 4.                                          C
C INARRM( 1 ) for VECM = 4.                                          C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     KLV       : Number of diagonals in DMATV matrix.               C
C     KLMF      : Number of diagonals in DMATM matrix.               C
C     NBN       : Number of bounding nodes.                          C
C     NCFM      : Leading dim. of SVFDCV, SVFDCM and FDCM arrays.    C
C     NR        : Number of radial grid nodes (velocity).            C
C     NRMF      : Number of radial grid nodes (magnetic field).      C
C     NHV       : Number of spherical harmonics in velocity vec.     C
C     NHM       : Number of spherical harmonics in magnetic vec.     C
C     NDCSV     : Number of finite difference schemes (velocity).    C
C     NDCSM     : Number of finite difference schemes (mag. field).  C
C                                                                    C
C     MHTV      : Dim ( * ). Harmonic type.  (Velocity)              C
C     MHLV      : Dim ( * ). Harmonic degree, l. (Velocity)          C
C     MHMV      : Dim ( * ). Harmonic order, m. (Velocity)           C
C     MHPV      : Pointer array for harmonics. If HMP( ih ) = is     C
C                  then 'is' is the finite difference scheme used    C
C                   to take derivatives of that harm. radial func.   C
C                   If MHP is negative, the harmonic is avoided and  C
C                   CASVDR moves on to the next harmonic.            C
C                   (Velocity function).                             C
C                                                                    C
C     MHTR      : Dim ( * ). Harmonic type for curl.                 C
C                                                                    C
C     MHTM      : Dim ( * ). Harmonic type.  (Magnetic field)        C
C     MHLM      : Dim ( * ). Harmonic degree, l. (Magnetic field)    C
C     MHMM      : Dim ( * ). Harmonic order, m. (Magnetic field)     C
C     MHPM      : Pointer array for harmonics. If HMP( ih ) = is     C
C                  then 'is' is the finite difference scheme used    C
C                   to take derivatives of that harm. radial func.   C
C                   If MHP is negative, the harmonic is avoided and  C
C                   CASVDR moves on to the next harmonic.            C
C                   (Magnetic field ).                               C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum value of m.                                C
C     M0        : Minimum value of non-zero m such that              C
C                  MOD( m, M0 ) = 0 for all m.                       C
C                                                                    C
C     NPHP      : Number of points in phi.                           C
C     NTHP      : Number of points in theta.                         C
C                                                                    C
C     N2V       : Dimension of vectors. (2nd dim of DMATV).          C
C     N2M       : Dimension of vectors. (2nd dim of DMATM).          C
C                                                                    C
C     NOIT      : On input. NOIT = Maximum number of iterations      C
C                 permitted in the functional iteration part of the  C
C                 process in the routine.                            C
C                 On output. NOIT is number of iterations taken      C
C                 UNLESS - maximum number is exceeded (NOIT = -1)    C
C                 OR the iteration norm,                             C
C                 DABS( |V_(j)^{i+1}| - |V_{j-1}^{i+1}| )            C
C                 is actually growing (NOIT = -2).                   C
C                                                                    C
C     LULOG     : Set to zero to suppress output of information.     C
C                 Otherwise, set to output chanel wth this number.   C
C                                                                    C
C     IPIVV     : Dim (N2V). Output from SFDDMF.                     C
C     IPIVM     : Dim (N2M). Output from MCVDMF.                     C
C                                                                    C
C     MHIBCV    : MHIBCV( is ) describes the inner boundary          C
C                  condition for scheme IS (velocity harmonics).     C
C                                                                    C
C  MHIBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( is ) = 7 --> insulating magnetic field.                    C
C                                                                    C
C     MHOBCV    : MHOBCV( is ) describes the outer boundary          C
C                  condition for scheme IS. (See above for key.)     C
C                  (velocity harmonics).                             C
C                                                                    C
C     MHIBCM    : MHIBCM( is ) describes the inner boundary          C
C                  condition for scheme IS. (See above for key.)     C
C                  (magnetic field harmonics).                       C
C                                                                    C
C     MHOBCM    : MHOBCM( is ) describes the outer boundary          C
C                  condition for scheme IS. (See above for key.)     C
C                  (magnetic field harmonics).                       C
C                                                                    C
C     NTS       : Number of toroidal singularities.                  C
C                 (See SFDDMF, DMWMF)                                C
C                                                                    C
C     NNDS      : Number of nodes (between 2 and NR)                 C
C                                                                    C
C     IN0       : Index arr. for constant vector. (See INDFUN).      C
C     MHT0      : Dim ( * ). Harmonic type. [Constant vector]        C
C     MHL0      : Dim ( * ). Harmonic degree, l. [Constant vector]   C
C     MHM0      : Dim ( * ). Harmonic order, m. [Constant vector]    C
C     MHP0      : Pointer array for harmonics. [Constant vector]     C
C                                                                    C
C                 Note that the pointer must indicate a finite diff. C
C                 scheme from the velocity and not the magnetic      C
C                 field set of schemes.                              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VECV      : Dim(*). Velocity vector at time step (i)           C
C     VECV1     : Dim(*). Velocity vector at time step (i+1)         C
C     VECM      : Dim(*). Mag. field vector at time step (i)         C
C     VECM1     : Dim(*). Mag. field vector at time step (i+1)       C
C                                                                    C
C     SVFDCV    : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, 5, NDCS ).                 C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     SVFDCM    : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, 3, NDCS ).                 C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
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
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Gauss weights computed bu GAUWTS. Dim ( NTHP )     C
C                                                                    C
C     XARRV     : Dim ( NR ). i^{th} element gives x_i.              C
C     XARRM     : Dim ( NRMF ). i^{th} element gives x_i.            C
C     DMATV     : Vel. diffusion matrix. Dim ( 3*KLV + 1, N2V ).     C
C     DMATM     : Mag. diffusion matrix. Dim ( 3*KLMF + 1, N2M ).    C
C     DPARS     : Dim( * ). Supplies double precision parameters.    C
C                                                                    C
C                  DPARS(  1 ) = DELTAT                              C
C                  DPARS(  2 ) = TOL                                 C
C                  DPARS(  3 ) = CFAC                                C
C                  DPARS(  4 ) = CA                                  C
C                  DPARS(  5 ) = CB1                                 C
C                  DPARS(  6 ) = CB2                                 C
C                  DPARS(  7 ) = CC                                  C
C                  DPARS(  8 ) = CD                                  C
C                  DPARS(  9 ) = CE                                  C
C                  DPARS( 10 ) = CF                                  C
C                  DPARS( 11 ) = CG                                  C
C                  DPARS( 12 ) = CH                                  C
C                  DPARS( 13 ) = CI                                  C
C                  DPARS( 14 ) = CJ                                  C
C                  DPARS( 15 ) = CK                                  C
C                  DPARS( 16 ) = CL                                  C
C                  DPARS( 17 ) = CM                                  C
C                  DPARS( 18 ) = SCAL - multiplier for the           C
C                                 constant term.                     C
C                                                                    C
C     VEC0      : Constant vector (see IN0, MHT0 etc.)               C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE IVICTS( KLV, KLMF, NBN, NCFM, NR, NRMF, NHV, NHM, 
     1 NDCSV, NDCSM, MHTV, MHLV, MHMV, MHPV, MHTR, MHTM, MHLM, MHMM,
     2 MHPM, LH, MMAX, M0, NPHP, NTHP, N2V, N2M, NOIT, LULOG, IPIVV,
     3 IPIVM, MHIBCV, MHOBCV, MHIBCM, MHOBCM, NTS, NNDS, VECV, VECV1,
     4 VECM, VECM1, SVFDCV, SVFDCM, FDCM, PA, DPA, GAUX, GAUW, XARRV,
     5 XARRM, DMATV, DMATM, DPARS, IN0, MHT0, MHL0, MHM0, MHP0, VEC0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KLV, KLMF, NBN, NCFM, NR, NRMF, NHV, NHM, NDCSV, NDCSM,
     1        MHTV( * ), MHLV( * ), MHMV( * ), MHPV( * ), MHTR( * ),
     2        MHTM( * ), MHLM( * ), MHMM( * ), MHPM( * ), LH, MMAX,
     3        M0, NPHP, NTHP, N2V, N2M, NOIT, LULOG
C
      INTEGER IPIVV( * ), IPIVM( * ), MHIBCV( * ), MHOBCV( * ),
     1        MHIBCM( * ), MHOBCM( * ), NTS, NNDS, IN0( * ),
     2        MHT0( * ), MHL0( * ), MHM0( * ), MHP0( * )
C
      DOUBLE PRECISION VECV( * ), VECV1( * ), VEC0( * ),
     1                 VECM( * ), VECM1( * )
C
      DOUBLE PRECISION SVFDCV( NCFM, NR, 5, NDCSV ), XARRV( * ),
     1                 SVFDCM( NCFM, NRMF, 3, NDCSM ), XARRM( * ),
     2                 DMATV( 3*KLV + 1, N2V ),
     3                 DMATM( 3*KLMF + 1, N2M )
C
      DOUBLE PRECISION FDCM( NCFM, NR, 1 ), DPARS( * ),
     1                 GAUX( NTHP ), GAUW( NTHP )
C
      DOUBLE PRECISION PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     1                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - work arrays -                              C
C The variables listed in this section are merely designed           C
C to reduce the number of arrays paseed into this routine.           C
C This has the disadvantage that some parameters, NPHP, LH           C
C etc. will not correspond to those set inside this routine.         C
C If they are significantly smaller, this may mean declaring         C
C more memory than necessary.                                        C
C If the arrays declared inside this routine are too small,          C
C execution will abort and you will need to recompile with           C
C new values. This should not present too much misery.               C
C____________________________________________________________________C
C
      INTEGER LHMAX, NPHMAX, NTHMAX, NRVMAX, NRMMAX, NNDMX, NHVMAX,
     1        NHMMAX, N2VMAX, N2MMAX, NTSMAX
      PARAMETER ( LHMAX = 62, NPHMAX = 128, NTHMAX = 64,
     1            NRVMAX = 60, NRMMAX = 80, NNDMX = 6,
     2            NHVMAX = 500, NHMMAX = 500, 
     3            N2VMAX = NRVMAX*NHVMAX, N2MMAX = NRMMAX*NHMMAX,
     4            NTSMAX = 3 )
C____________________________________________________________________C
C  Work array declarations                                           C
C____________________________________________________________________C
C
      INTEGER          IW1( NNDMX )
      DOUBLE PRECISION RQST1( LHMAX*(LHMAX+2), 3, NRVMAX ),
     1                 RQST2( LHMAX*(LHMAX+2), 3, NRVMAX ),
     2                 RQST3( LHMAX*(LHMAX+2), 3, NRVMAX ),
     3                 RQSTA( LHMAX*(LHMAX+2), 3, NRVMAX )
      DOUBLE PRECISION W1( NNDMX ), W2( NNDMX ), W3( NNDMX ),
     1                 U( N2VMAX, NTSMAX ), V( N2VMAX, NTSMAX )
      DOUBLE PRECISION SHC( LHMAX*(LHMAX+2) ),
     1                 DSHC( LHMAX*(LHMAX+2) ),
     2                 SF( NPHMAX, NTHMAX )
      DOUBLE PRECISION VF1( NPHMAX, NTHMAX, 3 ),
     1                 VF2( NPHMAX, NTHMAX, 3 ),
     2                 VF3( NPHMAX, NTHMAX, 3 )
      DOUBLE PRECISION DV0( N2VMAX ), DV1( N2VMAX ),
     1                 DV2( N2VMAX ), DV3( N2VMAX ),
     2                 DV4( N2VMAX ), DTVV( N2VMAX ),
     3                 RFIV( N2VMAX ), RFIV1( N2VMAX )
      DOUBLE PRECISION DM0( N2MMAX ), DM1( N2MMAX ),
     1                 DM2( N2MMAX ), DTVM( N2MMAX ),
     2                 RFIM( N2MMAX ), RFIM1( N2MMAX )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITMX, ILN, IRN, IHD, NDRVS, INARRV( 3 ), INARRM( 3 ),
     1        IOP, INC, NBANDV, NBANDM, IPIVH( 3 ), ILUDF, NRHS,
     2        INFO, ICLS, ICMP
      DOUBLE PRECISION DELTAT, TOL, CFAC, CA, CB1, CB2, CC, CD, CE,
     1                 CF, CG, CH, CI, CJ, CK, CL, CM, ZERO, AOLD,
     2                 AVEC1, DNRM2, FAC, DIFF, ODIFF, HMAT( 3, 3),
     2                 VTY( 3, 1 ), HVTY( 3, 1 ), SCAL
      CHARACTER *(1) TRANS
      PARAMETER ( IOP = 0, ZERO = 0.0d0, INC = 1, ILUDF = 2,
     1            TRANS = 'N', NRHS = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check that array-size-controlling parameters
C are in the correct range: i.e. that the above parameters
C are not exceeded.
C
      IF (     NR.GT.NRVMAX     .OR.     NRMF.GT.NRMMAX     .OR.
     1         NNDS.GT.NNDMX    .OR.     LH.GT.LHMAX        .OR.
     2         NTS.GT.NTSMAX    .OR.     NHV.GT.NHVMAX      .OR.
     3         NHM.GT.NHMMAX    .OR.     NTHP.GT.NTHMAX     .OR.
     4         NPHP.GT.NPHMAX                          ) THEN
        PRINT *,' Subroutine IVICTS.'
        PRINT *,' Some array size is wrong.'
        PRINT *,' Check from the following list and '
        PRINT *,' recompile routine. Steve Gibbons 3.4.2000 '
        PRINT *,' ------------------------------------------ '
        PRINT *,' NR    = ',NR    ,' NRVMAX = ', NRVMAX
        PRINT *,' NRMF  = ',NRMF  ,' NRMMAX = ', NRMMAX
        PRINT *,' NNDS  = ',NNDS  ,' NNDMX  = ', NNDMX
        PRINT *,' LH    = ',LH    ,' LHMAX  = ', LHMAX
        PRINT *,' NTS   = ',NTS   ,' NTSMAX = ', NTSMAX
        PRINT *,' NHV   = ',NHV   ,' NHVMAX = ', NHVMAX
        PRINT *,' NHM   = ',NHM   ,' NHMMAX = ', NHMMAX
        PRINT *,' NTHP  = ',NTHP  ,' NTHP   = ', NTHMAX
        PRINT *,' NPHP  = ',NPHP  ,' NPHP   = ', NPHMAX
        PRINT *,' N2V   = ',N2V   ,' N2VMAX = ', N2VMAX
        PRINT *,' N2M   = ',N2M   ,' N2MMAX = ', N2MMAX
        PRINT *,' ------------------------------------------ '
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NBANDV = 3*KLV + 1
      NBANDM = 3*KLMF + 1
      ODIFF = 1.0d8
C
      DELTAT  = DPARS(  1 )
      TOL     = DPARS(  2 )
      CFAC    = DPARS(  3 )
      CA      = DPARS(  4 )
      CB1     = DPARS(  5 )
      CB2     = DPARS(  6 )
      CC      = DPARS(  7 )
      CD      = DPARS(  8 )
      CE      = DPARS(  9 )
      CF      = DPARS( 10 )
      CG      = DPARS( 11 )
      CH      = DPARS( 12 )
      CI      = DPARS( 13 )
      CJ      = DPARS( 14 )
      CK      = DPARS( 15 )
      CL      = DPARS( 16 )
      CM      = DPARS( 17 )
      SCAL    = DPARS( 18 )
C
      INARRV( 1 ) = 4
      INARRV( 2 ) = NR
      INARRV( 3 ) = NHV
C
      INARRM( 1 ) = 4
      INARRM( 2 ) = NRMF
      INARRM( 3 ) = NHM
C
C Store the diffusion parts from the constant
C vector to the velocity forcing term in DTVV
C
      CALL VECOP( DTVV, ZERO, N2V, IOP )
C
C First the diffusion of vorticity
C
      ILN     = 2
      IRN     = NR - 1
      FAC     = CI*DELTAT*SCAL
      NDRVS   = 4
      CALL ASVLC( NR, NDCSV, VEC0, IN0, MHT0, MHL0, MHM0, MHP0,
     1            DTVV, INARRV, MHTR, MHLV, MHMV, FAC, NBN, NDRVS,
     2            NDRVS, ILN, IRN, SVFDCV, XARRV, NCFM )
C
C Now diffusion of temperature
C
      ILN     = 2
      IRN     = NR - 1
      FAC     = CD*DELTAT*SCAL
      NDRVS   = 4
      ICMP    = 3
      CALL ASVLP( NR, NDCSV, VEC0, IN0, MHT0, MHL0, MHM0, ICMP,
     1            MHP0, DTVV, INARRV, MHTR, MHLV, MHMV, FAC, NBN,
     2            NDRVS, NDRVS, ILN, IRN, SVFDCV, XARRV, NCFM )
C
C Get all derivatives of the step 'i' velocity vector
C
      ILN     = 2
      IRN     = NR - 1
      IHD     = 4
      NDRVS   = 4
      CALL CASVDR( VECV, ILN, IRN, NBN, IHD, NCFM, NR, NDRVS,
     1             NDRVS, INARRV, NDCSV, MHPV, SVFDCV, DV0, DV1,
     2             DV2, DV3, DV4 )
C
C Get all derivatives of the step 'i' magnetic field vector
C
      ILN     = 2
      IRN     = NRMF - 1
      IHD     = 2
      NDRVS   = 2
      CALL CASVD2( VECM, ILN, IRN, NBN, IHD, NCFM, NRMF, NDRVS,
     1             NDRVS, INARRM, NDCSM, MHPM, SVFDCM, DM0, DM1,
     2             DM2 )
C
C Add the diffusion parts of the velocity forcing term to DTVV
C
      CALL MCDRHF( NR, NHV, MHTV, MHLV, MHMV, MHTR, CA, CE,
     1             ZERO, CD, CI, ZERO, CFAC, DELTAT, DV0, DV1, DV2,
     2             DV3, DV4, DTVV, XARRV )
C
C Store the diffusion parts of the magnetic forcing term in DTVM
C
      CALL VECOP( DTVM, ZERO, N2M, IOP )
      CALL MCDRHF( NRMF, NHM, MHTM, MHLM, MHMM, MHTM, ZERO, ZERO,
     1             CK, ZERO, ZERO, CL, CFAC, DELTAT, DM0, DM1, DM2,
     2             DV3, DV4, DTVM, XARRM )
C
C Add to V0, V1, V2, V3 and V4 the contributions
C from constant vector.
C
      ICLS    = 0
      ILN     = 2
      IRN     = NR - 1
      IHD     = 4
      NDRVS   = 4
      CALL ACASVD ( ICLS, VEC0, ILN, IRN, NBN, IHD, NCFM, NR,
     1              NDRVS, NDRVS, NDCSV, IN0, MHT0, MHL0, MHM0,
     2              MHP0, INARRV, MHTV, MHLV, MHMV, SVFDCV, DV0,
     3              DV1, DV2, DV3, DV4, SCAL )
C
C Calculate the non-linear forcing terms for 'step i' terms
C Store this in the RFIV and RFIM arrays.
C
      CALL VECOP( RFIV, ZERO, N2V, IOP )
      CALL VECOP( RFIM, ZERO, N2M, IOP )
C
      CALL STCIFT( NR, INARRV, NRMF, INARRM, MHTV, MHLV, MHMV,
     1  MHTR, MHTM, MHLM, MHMM, LH, NCFM, NBN, MMAX, M0, NTHP, NPHP,
     2  CB1, CB2, CC, CH, CG, CF, CJ, CM, DV0, DV1, DM0, DM1, VF1,
     3  VF2, VF3, SF, PA, DPA, GAUX, GAUW, SHC, DSHC, FDCM, RFIV,
     4  RFIM, XARRV, XARRM, RQST1, RQST2, RQST3, RQSTA )
C
C Construct the right hand side vectors for finding predictor
C Store these in VECV1 and VECM1.
C Copy DTVV into VECV1 and DTVM into VECM1.
C
      CALL VECCP( DTVV, VECV1, N2V )
      CALL DAXPY( N2V, DELTAT, RFIV, INC, VECV1, INC )
C
      CALL VECCP( DTVM, VECM1, N2M )
      CALL DAXPY( N2M, DELTAT, RFIM, INC, VECM1, INC )
C
C Solve system to form predictor.
C First form U and V matrices for the Woodbury formula
C
      CALL DMWMF( N2V, NR, INARRV, MHTV, MHLV, MHPV, MHIBCV, MHOBCV,
     1            NNDS, IW1, NTS, XARRV, U, V, W1, W2, W3, VECV1 )
C
C Now solve for velocity (using Woodbury formula)
C
      CALL BMWDFS( NBANDV, N2V, NTS, KLV, KLV, KLV, IPIVV, IPIVH,
     1             DMATV, VECV1, U, V, HMAT, VTY, HVTY, DV4, ILUDF )
C
C Complete the velocity solution vector
C
      NDRVS = 4
      CALL ASVCPL( VECV1, NR, NDCSV, INARRV, MHPV, MHIBCV, MHOBCV,
     1             NCFM, NDRVS, NDRVS, NBN, SVFDCV )
      AVEC1 = DNRM2( N2V, VECV1, INC )
C
C AVEC1 is now the norm of the (predictor) velocity estimate.
C Now solve for magnetic field (we do not need the Woodbury
C formula for this as there are no singularities)
C
      CALL DGBTRS( TRANS, N2M, KLMF, KLMF, NRHS, DMATM, NBANDM,
     1             IPIVM, VECM1, N2M, INFO )
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine IVICTS.'
        PRINT *,' DGBTRS called and INFO returned value ',INFO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now complete the magnetic field solution vector
C
      NDRVS = 2
      CALL ASVCPL( VECM1, NRMF, NDCSM, INARRM, MHPM, MHIBCM, MHOBCM,
     1             NCFM, NDRVS, NDRVS, NBN, SVFDCM )
      AVEC1 = AVEC1 + DNRM2( N2M, VECM1, INC )
C
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 81 ) AVEC1
 81   FORMAT('ivicts: Predictor norm = ',1PD16.7)
C
      ITMX = NOIT
      NOIT = 0
 50   CONTINUE
      NOIT = NOIT + 1
      IF ( NOIT.GT.ITMX ) THEN
        IF ( LULOG.NE.0 ) WRITE ( LULOG, 83 )
 83     FORMAT('ivicts: Max. iterations exceeded.')
        NOIT = -1
        RETURN
      ENDIF
C
C 'Predictor' is now stored in VECV1 and VECM1
C Calculate derivatives of the step 'i+1' velocity vector
C
      ILN     = 2
      IRN     = NR - 1
      IHD     = 4
      NDRVS   = 4
      CALL CASVDR( VECV1, ILN, IRN, NBN, IHD, NCFM, NR, NDRVS,
     1             NDRVS, INARRV, NDCSV, MHPV, SVFDCV, DV0, DV1,
     2             DV2, DV3, DV4 )
C
C Get all derivatives of the step 'i+1' magnetic field vector
C
      ILN     = 2
      IRN     = NRMF - 1
      IHD     = 2
      NDRVS   = 2
      CALL CASVD2( VECM1, ILN, IRN, NBN, IHD, NCFM, NRMF, NDRVS,
     1             NDRVS, INARRM, NDCSM, MHPM, SVFDCM, DM0, DM1,
     2             DM2 )
C
C Add to V0, V1, V2, V3 and V4 the contributions
C from constant vector.
C
      ICLS    = 0
      ILN     = 2
      IRN     = NR - 1
      IHD     = 4
      NDRVS   = 4
      CALL ACASVD ( ICLS, VEC0, ILN, IRN, NBN, IHD, NCFM, NR,
     1              NDRVS, NDRVS, NDCSV, IN0, MHT0, MHL0, MHM0,
     2              MHP0, INARRV, MHTV, MHLV, MHMV, SVFDCV, DV0,
     3              DV1, DV2, DV3, DV4, SCAL )
C
C Calculate the non-linear forcing terms for 'step i+1' terms
C Store this in the RFIV1 and RFIM1 arrays.
C
      CALL VECOP( RFIV1, ZERO, N2V, IOP )
      CALL VECOP( RFIM1, ZERO, N2M, IOP )
C
      CALL STCIFT( NR, INARRV, NRMF, INARRM, MHTV, MHLV, MHMV,
     1  MHTR, MHTM, MHLM, MHMM, LH, NCFM, NBN, MMAX, M0, NTHP, NPHP,
     2  CB1, CB2, CC, CH, CG, CF, CJ, CM, DV0, DV1, DM0, DM1, VF1,
     3  VF2, VF3, SF, PA, DPA, GAUX, GAUW, SHC, DSHC, FDCM, RFIV1,
     4  RFIM1, XARRV, XARRM, RQST1, RQST2, RQST3, RQSTA )
C
C Form the new right hand side vectors
C
      CALL VECCP( DTVV, VECV1, N2V )
      CALL VECCP( DTVM, VECM1, N2M )
      FAC = CFAC*DELTAT
      CALL DAXPY( N2V, FAC, RFIV, INC, VECV1, INC )
      CALL DAXPY( N2M, FAC, RFIM, INC, VECM1, INC )
      FAC = (1.0d0 - CFAC)*DELTAT
      CALL DAXPY( N2V, FAC, RFIV1, INC, VECV1, INC )
      CALL DAXPY( N2M, FAC, RFIM1, INC, VECM1, INC )
C
C Solve the system for new iteration of solution.
C First form U and V matrices for the Woodbury formula
C
      CALL DMWMF( N2V, NR, INARRV, MHTV, MHLV, MHPV, MHIBCV, MHOBCV,
     1            NNDS, IW1, NTS, XARRV, U, V, W1, W2, W3, VECV1 )
C
C Now solve for velocity using Woodbury formula
C
      CALL BMWDFS( NBANDV, N2V, NTS, KLV, KLV, KLV, IPIVV, IPIVH,
     1             DMATV, VECV1, U, V, HMAT, VTY, HVTY, DV4, ILUDF )
C
C Complete the velocity vector
C
      NDRVS = 4
      CALL ASVCPL( VECV1, NR, NDCSV, INARRV, MHPV, MHIBCV, MHOBCV,
     1             NCFM, NDRVS, NDRVS, NBN, SVFDCV )
      AOLD  = AVEC1
      AVEC1 = DNRM2( N2V, VECV1, INC )
C
C Now solve directly for the magnetic field
C
      CALL DGBTRS( TRANS, N2M, KLMF, KLMF, NRHS, DMATM, NBANDM,
     1             IPIVM, VECM1, N2M, INFO )
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine IVICTS.'
        PRINT *,' DGBTRS called and INFO returned value ',INFO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now complete the magnetic field solution vector
C
      NDRVS = 2
      CALL ASVCPL( VECM1, NRMF, NDCSM, INARRM, MHPM, MHIBCM, MHOBCM,
     1             NCFM, NDRVS, NDRVS, NBN, SVFDCM )
      AVEC1 = AVEC1 + DNRM2( N2M, VECM1, INC )
C
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 84 ) NOIT, AVEC1
 84   FORMAT('ivicts: Iteration ',I4,' Solution norm = ',1PD16.7)
      DIFF = DABS( AVEC1 - AOLD )
      IF ( DIFF.LT.TOL ) THEN
        IF ( LULOG.NE.0 ) WRITE ( LULOG, 85 )
 85     FORMAT('Solution converged.')
        RETURN
      ENDIF
C
C ok - so our reiteration wasn't sufficiently
C close to return - check to see if we are getting worse
C
      IF ( DIFF.GT.ODIFF .AND. NOIT.GT.3 ) THEN
        NOIT = -2
        RETURN
      ENDIF
      ODIFF = DIFF
C
      GOTO 50
C
      END
C*********************************************************************
