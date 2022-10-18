C*********************************************************************
C subroutine New Imposed Vector Time Step Routine ********************
C            -   -       -      -    -    -       ********************
C Steve Gibbons Mon May 15 14:12:51 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C We supply a solution vector VEC at time step i, along with the     C
C diffusion matrix DMAT which is pre-calculated by the routine       C
C SFDDMF. NIVTSR will advance the solution by one time step (size    C
C DELTAT).                                                           C
C Differs from TSIDTS in that the diffusion matrix must be formed    C
C by SFDDMF as opposed to MCVDMF and the linear system must be       C
C carefully solved by the Woodbury formula (routine BMWDFS) after    C
C the matrices U and V are formed by DFWMF.                          C
C                                                                    C
C The equations are integrated subject to the addition of a          C
C constant vector contained in VEC0.                                 C
C                                                                    C
C for instance U = u + FAC*u_0     etc.                              C
C                                                                    C
C VEC0 is indexed by IN0, MHT0, MHL0, MHM0 and MHP0                  C
C All harmonics indicated by these arrays must also be contained     C
C somewhere in the standard solution vector in INARR, MHT, MHL, MHM  C
C and MHP.                                                           C
C                                                                    C
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
C INARR( 1 ) for VEC = 4.                                            C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     KL        : Number of diagonals in matrix.                     C
C     NBN       : Number of bounding nodes.                          C
C     NCFM      : Leading dimension of SVFDC and FDCM arrays.        C
C     NR        : Number of radial grid nodes.                       C
C     NDRVM     : Max deriv.s stored in SVFDC and FDCM arrays.       C
C     NH        : Number of spherical harmonics in solution vec.     C
C     NDCS      : Number of finite difference schemes.               C
C                                                                    C
C     MHT       : Dim ( * ). Harmonic type.                          C
C     MHL       : Dim ( * ). Harmonic degree, l.                     C
C     MHM       : Dim ( * ). Harmonic order, m.                      C
C     MHP       : Pointer array for harmonics. If HMP( ih ) = is     C
C                  then 'is' is the finite difference scheme used    C
C                   to take derivatives of that harm. radial func.   C
C                   If MHP is negative, the harmonic is avoided and  C
C                   CASVDR moves on to the next harmonic.            C
C     MHTR      : Dim ( * ). Harmonic type for curl.                 C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum value of m.                                C
C     M0        : Minimum value of non-zero m such that              C
C                  MOD( m, M0 ) = 0 for all m.                       C
C     N2        : Dimension of vectors. (2nd dim of DMAT).           C
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
C     IPIV      : Dim (N2). Output from MCVDMF.                      C
C                                                                    C
C     MHIBC     : MHIBC( is ) describes the inner boundary           C
C                  condition for scheme IS.                          C
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
C     MHOBC     : MHOBC( is ) describes the outer boundary           C
C                  condition for scheme IS. (See above for key.)     C
C                                                                    C
C     NTS       : Number of toroidal singularities.                  C
C                 (See SFDDMF, DMWMF )                               C
C                                                                    C
C     NNDS      : Number of nodes (between 2 and NR)                 C
C     IW1       : Dim (NNDS). Work array.                            C
C                                                                    C
C     IN0       : Index arr. for constant vector. (See INDFUN).      C
C     MHT0      : Dim ( * ). Harmonic type. [Constant vector]        C
C     MHL0      : Dim ( * ). Harmonic degree, l. [Constant vector]   C
C     MHM0      : Dim ( * ). Harmonic order, m. [Constant vector]    C
C     MHP0      : Pointer array for harmonics. [Constant vector]     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Dim(*). Solution vector at time step (i)           C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
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
C     SHC       : Work array. Dim ( LH*(LH+2) ).                     C
C     DSHC      : Work array. Dim ( LH*(LH+2) ).                     C
C     RVF1      : Work arr. dim. ( NR, NPHP, NTHP, 3 )               C
C     RVF2      : Work arr. dim. ( NR, NPHP, NTHP, 3 )               C
C     RVF3      : Work arr. dim. ( NR, NPHP, NTHP, 3 )               C
C     SF        : Work arr. dim. ( NPHP, NTHP )                      C
C     VF        : Work arr. dim. ( NPHP, NTHP, 3 )                   C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     DMAT      : Diffusion matrix. Dim ( 3*KL + 1, N2 ).            C
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
C     V0        : Dim ( * ). Work array.                             C
C     V1        : Dim ( * ). Work array.                             C
C     V2        : Dim ( * ). Work array.                             C
C     V3        : Dim ( * ). Work array.                             C
C     V4        : Dim ( * ). Work array.                             C
C     DTV       : Dim ( * ). Work array.                             C
C     RFI       : Dim ( * ). Work array.                             C
C     RFI1      : Dim ( * ). Work array.                             C
C     VEC1      : Dim ( * ). Work array.                             C
C                                                                    C
C     RQST1     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C     RQST2     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C     RQST3     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C     RQSTA     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C                                                                    C
C     W1        : Work array. Dim ( NNDS )                           C
C     W2        : Work array. Dim ( NNDS )                           C
C     W3        : Work array. Dim ( NNDS, NNDS )                     C
C     U         : U matrix :  Dim ( N2, NTS )                        C
C     V         : V matrix :  Dim ( N2, NTS )                        C
C     VEC0      : Constant vector (see IN0, MHT0 etc.)               C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NIVTSR( KL, NBN, NCFM, NR, NDRVM, NH, NDCS, MHT, MHL,
     1   MHM, MHP, MHTR, LH, MMAX, M0, NPHP, NTHP, VEC, SVFDC, FDCM,
     2   PA, DPA, GAUX, GAUW,SHC,DSHC,RVF1,RVF2,RVF3,SF, VF, XARR, N2,
     3   NOIT, LULOG, IPIV, MHIBC, MHOBC, DMAT, DPARS, V0, V1, V2, V3,
     4   V4, DTV, RFI, RFI1, VEC1, RQST1, RQST2, RQST3, RQSTA, NTS,
     5   NNDS, IW1, W1, W2, W3, U, V, IN0,MHT0,MHL0,MHM0,MHP0,VEC0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KL, NBN, NCFM, NR, NDRVM, NH, NDCS, MHT( * ), MHL( * ),
     1        MHM( * ), MHP( * ), MHTR( * ), LH, MMAX, M0, NPHP, NTHP,
     2        N2, NOIT, LULOG, IPIV( * ), MHIBC( * ), MHOBC( * ), NTS,
     3        NNDS, IW1( NNDS ), IN0( * ), MHT0( * ), MHL0( * ),
     4        MHM0( * ), MHP0( * )
      DOUBLE PRECISION VEC( * ), SVFDC( NCFM, NR, NDRVM+1, NDCS ),
     1                 XARR( * ), DMAT( 3*KL + 1, N2 ),
     2                 DPARS( * ), V0( * ), V1( * ), V2( * ),
     3                 V3( * ), V4( * ), DTV( * ), RFI( * ),
     2                 RFI1( * ), VEC1( * ), VEC0( * )
      DOUBLE PRECISION FDCM( NCFM, NR, 1 ), SF( NPHP, NTHP ),
     1                 RVF1( NR, NPHP, NTHP, 3 ), GAUX( NTHP ),
     2                 RVF2( NR, NPHP, NTHP, 3 ), GAUW( NTHP ),
     3                 RVF3( NR, NPHP, NTHP, 3 ), VF( NPHP, NTHP, 3 )
      DOUBLE PRECISION SHC( LH*(LH+2) ), DSHC( LH*(LH+2) ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
      DOUBLE PRECISION RQST1( LH*(LH+2), 3, NR ),
     1                 RQST2( LH*(LH+2), 3, NR ),
     2                 RQST3( LH*(LH+2), 3, NR ),
     3                 RQSTA( LH*(LH+2), 3, NR )
      DOUBLE PRECISION W1( NNDS ), W2( NNDS ), W3( NNDS ),
     1                 U( N2, NTS ), V( N2, NTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITMX, ILN, IRN, IHD, NDRVS, INARR( 3 ), IOP, INC,
     1        NBAND, IPIVH( 3 ), ILUDF, ICLS, ICMP
      DOUBLE PRECISION DELTAT, TOL, CFAC, CA, CB1, CB2, CC, CD, CE,
     1                 CF, CG, CH, CI, CJ, CK, CL, CM, ZERO, AOLD,
     2                 AVEC1, DNRM2, FAC, DIFF, ODIFF, HMAT( 3, 3),
     2                 VTY( 3, 1 ), HVTY( 3, 1 ), SCAL
      PARAMETER ( IOP = 0, ZERO = 0.0d0, INC = 1, ILUDF = 2 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NBAND = 3*KL + 1
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
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C
C Store the diffusion parts from the constant vector
C to the forcing term in DTV
C
      CALL VECOP( DTV, ZERO, N2, IOP )
C
C First the diffusion of vorticity
C
      ILN     = 2
      IRN     = NR - 1
      FAC     = CI*DELTAT*SCAL
      NDRVS   = 4
      CALL ASVLC( NR, NDCS, VEC0, IN0, MHT0, MHL0, MHM0, MHP0,
     1            DTV, INARR, MHTR, MHL, MHM, FAC, NBN, NDRVS,
     2            NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C
C Now diffusion of temperature
C
      ILN     = 2
      IRN     = NR - 1
      FAC     = CD*DELTAT*SCAL
      NDRVS   = 4
      ICMP    = 3
      CALL ASVLP( NR, NDCS, VEC0, IN0, MHT0, MHL0, MHM0, ICMP,
     1            MHP0, DTV, INARR, MHTR, MHL, MHM, FAC, NBN,
     2            NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C
C Now diffusion of poloidal magnetic field
C
      ILN     = 2
      IRN     = NR - 1
      FAC     = CL*DELTAT*SCAL
      NDRVS   = 4
      ICMP    = 4
      CALL ASVLP( NR, NDCS, VEC0, IN0, MHT0, MHL0, MHM0, ICMP,
     1            MHP0, DTV, INARR, MHTR, MHL, MHM, FAC, NBN,
     2            NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C
C Now diffusion of toroidal magnetic field
C
      ILN     = 2
      IRN     = NR - 1
      FAC     = CL*DELTAT*SCAL
      NDRVS   = 4
      ICMP    = 5
      CALL ASVLP( NR, NDCS, VEC0, IN0, MHT0, MHL0, MHM0, ICMP,
     1            MHP0, DTV, INARR, MHTR, MHL, MHM, FAC, NBN,
     2            NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C
C Get all derivatives of the step 'i' solution vector
C
      ILN     = 2
      IRN     = NR - 1
      IHD     = 4
      NDRVS   = 4
      CALL CASVDR( VEC, ILN, IRN, NBN, IHD, NCFM, NR, NDRVS,
     1             NDRVM, INARR, NDCS, MHP, SVFDC, V0, V1,
     2             V2, V3, V4 )
C
C Add the diffusion parts from the changing vector
C to the forcing term to DTV
C
      CALL MCDRHF( NR, NH, MHT, MHL, MHM, MHTR, CA, CE, CK, CD,
     1             CI, CL, CFAC, DELTAT, V0, V1, V2, V3, V4,
     2             DTV, XARR )
C
C Add to V0, V1, V2, V3 and V4 the contributions
C from constant vector.
C
      ICLS = 0
      CALL ACASVD ( ICLS, VEC0, ILN, IRN, NBN, IHD, NCFM, NR,
     1              NDRVS, NDRVM, NDCS, IN0, MHT0, MHL0, MHM0,
     2              MHP0, INARR, MHT, MHL, MHM, SVFDC, V0, V1,
     3              V2, V3, V4, SCAL )
C
C Calculate the non-linear forcing terms for 'step i' terms
C Store this in the RFI array.
C
      CALL VECOP( RFI, ZERO, N2, IOP )
      CALL NTSVFT( NR, INARR, MHT, MHL, MHM, MHTR, LH, NCFM,
     1             NBN, MMAX, M0, NTHP, NPHP, CB1, CB2, CC, CH,
     2             CG, CF, CJ, CM, V0, V1, RVF1, RVF2, RVF3, SF,
     3             VF, PA, DPA, GAUX, GAUW, SHC, DSHC, FDCM, RFI,
     4             XARR, RQST1, RQST2, RQST3, RQSTA )
C
C Construct the right hand side vector for finding predictor
C Store this in VEC1.  Copy DTV into VEC1.
C
      CALL VECCP( DTV, VEC1, N2 )
      CALL DAXPY( N2, DELTAT, RFI, INC, VEC1, INC )
C
C Solve system to form predictor.
C First form U and V matrices for the Woodbury formula
C
      CALL DMWMF( N2, NR, INARR, MHT, MHL, MHP, MHIBC, MHOBC,
     1            NNDS, IW1, NTS, XARR, U, V, W1, W2, W3, VEC1 )
C
C Now solve using Woodbury formula
C
      CALL BMWDFS( NBAND, N2, NTS, KL, KL, KL, IPIV, IPIVH, DMAT,
     1             VEC1, U, V, HMAT, VTY, HVTY, V4, ILUDF )
C
      CALL ASVCPL( VEC1, NR, NDCS, INARR, MHP, MHIBC, MHOBC,
     1             NCFM, NDRVS, NDRVM, NBN, SVFDC )
      AVEC1 = DNRM2( N2, VEC1, INC )
C
C AVEC1 is now the norm of the (predictor) estimate.
C
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 81 ) AVEC1
 81   FORMAT('ivsits: Predictor norm = ',1PD16.7)
C
      ITMX = NOIT
      NOIT = 0
 50   CONTINUE
      NOIT = NOIT + 1
      IF ( NOIT.GT.ITMX ) THEN
        IF ( LULOG.NE.0 ) WRITE ( LULOG, 83 )
 83     FORMAT('ivsits: Max. iterations exceeded.')
        NOIT = -1
        RETURN
      ENDIF
C
C 'Predictor' is now stored in VEC1
C Calculate derivatives
C
      ILN     = 2
      IRN     = NR - 1
      IHD     = 4
      NDRVS   = 4
      CALL CASVDR( VEC1, ILN, IRN, NBN, IHD, NCFM, NR, NDRVS,
     1             NDRVM, INARR, NDCS, MHP, SVFDC, V0, V1,
     2             V2, V3, V4 )     
C
C Add to V0, V1, V2, V3 and V4 the contributions
C from constant vector.
C
      ICLS = 0
      CALL ACASVD ( ICLS, VEC0, ILN, IRN, NBN, IHD, NCFM, NR,
     1              NDRVS, NDRVM, NDCS, IN0, MHT0, MHL0, MHM0,
     2              MHP0, INARR, MHT, MHL, MHM, SVFDC, V0, V1,
     3              V2, V3, V4, SCAL )
C
C Calculate the non-linear forcing terms for 'step i+1' terms
C Store this in the RFI1 array.
C
      CALL VECOP( RFI1, ZERO, N2, IOP )
      CALL NTSVFT( NR, INARR, MHT, MHL, MHM, MHTR, LH, NCFM,
     1             NBN, MMAX, M0, NTHP, NPHP, CB1, CB2, CC, CH,
     2             CG, CF, CJ, CM, V0, V1, RVF1, RVF2, RVF3, SF,
     3             VF, PA, DPA, GAUX, GAUW, SHC, DSHC, FDCM, RFI1,
     4             XARR, RQST1, RQST2, RQST3, RQSTA )
C
C Form the new right hand side vector
C
      CALL VECCP( DTV, VEC1, N2 )
      FAC = CFAC*DELTAT
      CALL DAXPY( N2, FAC, RFI, INC, VEC1, INC )
      FAC = (1.0d0 - CFAC)*DELTAT
      CALL DAXPY( N2, FAC, RFI1, INC, VEC1, INC )
C
C Solve the system for new iteration of solution.
C First form U and V matrices for the Woodbury formula
C
      CALL DMWMF( N2, NR, INARR, MHT, MHL, MHP, MHIBC, MHOBC,
     1            NNDS, IW1, NTS, XARR, U, V, W1, W2, W3, VEC1 )
C
C Now solve using Woodbury formula
C
      CALL BMWDFS( NBAND, N2, NTS, KL, KL, KL, IPIV, IPIVH, DMAT,
     1             VEC1, U, V, HMAT, VTY, HVTY, V4, ILUDF )
C
      CALL ASVCPL( VEC1, NR, NDCS, INARR, MHP, MHIBC, MHOBC,
     1             NCFM, NDRVS, NDRVM, NBN, SVFDC )
      AOLD  = AVEC1
      AVEC1 = DNRM2( N2, VEC1, INC )
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 84 ) NOIT, AVEC1
 84   FORMAT('ivsits: Iteration ',I4,' Solution norm = ',1PD16.7)
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
