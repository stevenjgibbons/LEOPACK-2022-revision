C*********************************************************************
C subroutine Semi-Implicit Magnetic and Convection Time Step routine *
C            -    -        -            -          -    -            *
C Steve Gibbons Fri Feb 11 09:32:14 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C We supply a solution vector VEC at time step i, along with the     C
C diffusion matrix DMAT which is pre-calculated by the routine       C
C MCVDMF. SIMCTS will advance the solution by one time step (size    C
C DELTAT).                                                           C
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
C     NCFM      : Leading dimension of SVFDC array.                  C
C     NR        : Number of radial grid nodes.                       C
C     NDRVM     : Max deriv.s stored in SVFDC array.                 C
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
C                                                                    C
C     ISLARR    : Starting locations from MCNLIC. Dim ( 5 ).         C
C     NVIARR    : Number of coefficients from MCNLIC. Dim ( 5 ).     C
C     IHNA      : Number of alpha harmonics. Dim ( * ) [MCNLIC]      C
C     IHNB      : Number of beta  harmonics. Dim ( * ) [MCNLIC]      C
C     IHNG      : Number of gamma harmonics. Dim ( * ) [MCNLIC]      C
C                                                                    C
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
C     LULOG     : 0 to suppress output. File number otherwise.       C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Dim(*). Solution vector at time step (i)           C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     CVI       : Dim(*). Interaction coefficients. [MCNLIC]         C
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
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. ( * ).       C
C                 = 'CQSS', 'CQST' etc. according to the corresp.    C
C                 vector interaction.                                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SIMCTS( KL, NBN, NCFM, NR, NDRVM, NH, NDCS, MHT, MHL,
     1             MHM, MHP, MHTR, VEC, SVFDC, XARR, ISLARR, NVIARR,
     2             IHNA, IHNB, IHNG, N2, NOIT, LULOG, IPIV, MHIBC,
     3             MHOBC, TVHI, CVI, DMAT, DPARS, V0, V1, V2, V3, V4,
     4             DTV, RFI, RFI1, VEC1 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KL, NBN, NCFM, NR, NDRVM, NH, NDCS, MHT( * ), MHL( * ),
     1        MHM( * ), MHP( * ), MHTR( * ), ISLARR( 5 ), NVIARR( 5 ),
     2        IHNA( * ), IHNB( * ), IHNG( * ), N2, NOIT, LULOG,
     3        IPIV( * ), MHIBC( * ), MHOBC( * )
      DOUBLE PRECISION VEC( * ), SVFDC( NCFM, NR, NDRVM+1, NDCS ),
     1                 XARR( * ), CVI( * ), DMAT( 3*KL + 1, N2 ),
     2                 DPARS( * ), V0( * ), V1( * ), V2( * ),
     3                 V3( * ), V4( * ), DTV( * ), RFI( * ),
     2                 RFI1( * ), VEC1( * )
      CHARACTER *(4) TVHI( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITMX, ILN, IRN, IHD, NDRVS, INARR( 3 ), IOP, INC,
     1        INFO, NBAND, NRHS
      DOUBLE PRECISION DELTAT, TOL, CFAC, CA, CB1, CB2, CC, CD, CE,
     1                 CF, CG, CH, CI, CJ, CK, CL, CM, ZERO, AOLD,
     2                 AVEC1, DNRM2, FAC, DIFF, ODIFF
      CHARACTER *(1) TRANS
      PARAMETER ( IOP = 0, ZERO = 0.0d0, INC = 1, TRANS = 'N',
     1            NRHS = 1 )
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
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH
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
C Store the diffusion parts of the forcing term in DTV
C
      CALL VECOP( DTV, ZERO, N2, IOP )
      CALL MCDRHF( NR, NH, MHT, MHL, MHM, MHTR, CA, CE, CK, CD,
     1             CI, CL, CFAC, DELTAT, V0, V1, V2, V3, V4,
     2             DTV, XARR )
C
C Calculate the non-linear forcing terms for 'step i' terms
C Store this in the RFI array.
C
      CALL VECOP( RFI, ZERO, N2, IOP )
      CALL SSVFTA( INARR, MHT, MHL, MHM, MHTR, ISLARR, NVIARR,
     1             IHNA, IHNB, IHNG, CB1, CB2, CC, CH, CG, CF,
     2             CJ, CM, V0, V1, V2, V3, RFI, XARR, CVI, TVHI)
C
C Construct the right hand side vector for finding predictor
C Store this in VEC1.  Copy DTV into VEC1.
C
      CALL VECCP( DTV, VEC1, N2 )
      CALL DAXPY( N2, DELTAT, RFI, INC, VEC1, INC )
C
C Solve system to form predictor.
C
      CALL DGBTRS( TRANS, N2, KL, KL, NRHS, DMAT, NBAND,
     1             IPIV, VEC1, N2, INFO )
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine SIMCTS.'
        PRINT *,' DGBTRS called and INFO returned value ',INFO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      CALL ASVCPL( VEC1, NR, NDCS, INARR, MHP, MHIBC, MHOBC,
     1             NCFM, NDRVS, NDRVM, NBN, SVFDC )
      AVEC1 = DNRM2( N2, VEC1, INC )
C
C AVEC1 is now the norm of the (predictor) estimate.
C
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 81 ) AVEC1
 81   FORMAT('simcts: Predictor norm = ',1PD16.7)
C
      ITMX = NOIT
      NOIT = 0
 50   CONTINUE
      NOIT = NOIT + 1
      IF ( NOIT.GT.ITMX ) THEN
        IF ( LULOG.NE.0 ) WRITE ( LULOG, 83 )
 83     FORMAT('simcts: Max. iterations exceeded.')
        NOIT = -1
        RETURN
      ENDIF
C
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 82 ) NOIT
 82   FORMAT('simcts: Iteration number = ',I4)
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
C Calculate the non-linear forcing terms for 'step i+1' terms
C Store this in the RFI1 array.
C
      CALL VECOP( RFI1, ZERO, N2, IOP )
      CALL SSVFTA( INARR, MHT, MHL, MHM, MHTR, ISLARR, NVIARR,
     1             IHNA, IHNB, IHNG, CB1, CB2, CC, CH, CG, CF,
     2             CJ, CM, V0, V1, V2, V3, RFI1, XARR, CVI, TVHI)
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
C
      CALL DGBTRS( TRANS, N2, KL, KL, NRHS, DMAT, NBAND,
     1             IPIV, VEC1, N2, INFO )
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine SIMCTS.'
        PRINT *,' DGBTRS called and INFO returned value ',INFO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      CALL ASVCPL( VEC1, NR, NDCS, INARR, MHP, MHIBC, MHOBC,
     1             NCFM, NDRVS, NDRVM, NBN, SVFDC )
      AOLD  = AVEC1
      AVEC1 = DNRM2( N2, VEC1, INC )
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 84 ) AVEC1
 84   FORMAT('simcts: Solution norm = ',1PD16.7)
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
