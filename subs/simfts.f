C*********************************************************************
C subroutine Semi-Implicit Magnetic Field Time Step routine **********
C            -    -        -        -     -    -            **********
C Steve Gibbons Sun Nov  7 11:43:03 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C We supply the magnetic field, BI, and the velocity, VI, evaluated  C
C at the time step i ( TIME = t say ). We also supply the velocity   C
C VIP1 evaluated at time step i+1 ( TIME = t + deltat ).             C
C                                                                    C
C SIMFTS will solve for BIP1, the magnetic field at TIME t + deltat  C
C by solving the equation                                            C
C                                                                    C
C B^{i+1} - B^{i}                                                    C
C --------------- =  c \nabla^{2} B^{i} + c RM \curl (V^{i} x B^{i}) C
C    dt               + (1-c) \nabla^{2} B^{i+1}                     C
C                       + (1-c) RM \curl ( V^{i+1} x B^{i+1} )       C
C                                                                    C
C The diffusive part is trivially implemented by forming the DMAT    C
C matrix with the routine MFDMFM with NBAND = 3*KL + 1, DC = 1 and   C
C DELTAT = dt. CFAC is the factor c: c = 1 --> fully explicit,       C
C c = 0 --> fully implicit and c = 0.5, Crack-Nicolson.              C
C                                                                    C
C Note DMAT MUST be treated with MFDMFM BEFORE calling SIMFTS.       C
C                                                                    C
C The non-linear implicit term is solved for by repeated iteration.  C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     NMH       : Number of magnetic field harmonics.                C
C     KL        : Number of lower diagonals in matrix.               C
C                 (The leading order of the matrix is 3*KL+1. )      C
C     ILN       : Left-most node for evaluating RHS.                 C
C     IRN       : Right-most node for evaluating RHS.                C
C     INVA      : Integer array dimension ( * ). Details of Velocity C
C                  Elements may be arbitrary except for              C
C                  INVA ( 1 ) = IFORMV - flag for vector format.     C
C                                        See INDFUN                  C
C                  INVA ( 2 ) = NRR. Must be consistent with NR.     C
C                  INVA ( 3 ) = NHV = total number of radial func.s  C
C                                                                    C
C     MT        : Array length ( * ) - atleast length NMH            C
C                                                                    C
C         MT/VT( IH ) = 1 --> harmonic is poloidal velocity          C
C         MT/VT( IH ) = 2 --> harmonic is toroidal velocity          C
C         MT/VT( IH ) = 3 --> harmonic is temperature.               C
C         MT/VT( IH ) = 4 --> harmonic is poloidal magnetic field    C
C         MT/VT( IH ) = 5 --> harmonic is toroidal magnetic field    C
C                                                                    C
C     ML        : Array length ( * ) - atleast length NMH            C
C                  Sph. harm. degree, l.                             C
C     MM        : Array length ( * ) - atleast length NMH            C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MP        : Array length ( * ) - atleast length NMH            C
C                  Pointer array to finite difference coefficients.  C
C                  MP ( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C   MT, ML, MM and MP all relate to the magnetic field:              C
C   VT, VL, VM and VP are the equivalent dimensions of the velocity. C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                                                                    C
C     M0        : Lowest non-zero mode if simple multiplicity        C
C                  applies. (Put 1 if this is safest).               C
C                                                                    C
C     MMAX      : Maximum mode number to be used.                    C
C                                                                    C
C     NTHP      : Number of points in theta                          C
C     NPHP      : Number of points in phi                            C
C                                                                    C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dimension of FDCM. At least (2*NBN+1)      C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     IRES      : The resolution flag.                               C
C                 If IRES = 1, harmonics with l greater than LH      C
C                 or m greater than MMAX are simply ignored.         C
C                 If IRES = 2, the program aborts if LH or MMAX      C
C                 is exceeded.                                       C
C                                                                    C
C     NBN       : Number of bounding nodes.                          C
C                                                                    C
C     IPIV      : Dim ( * ) aleast NMH*NR. Output from MFDMFM.       C
C                                                                    C
C     MHBCS     : Dim ( NDCS ). MHBCS( is ) describes the boundary   C
C                 condition for finite difference scheme is.         C
C                 See SVFDCF for details. Note that inner and        C
C                 outer boundary conditions are assumed to be        C
C                 identical.                                         C
C                                                                    C
C     ITMX      : On input: Max. number of iterations permitted.     C
C                 On output: No. of iterations taken.                C
C                 If ITMX is returned = -1 then successive           C
C                 iterations on BIP1 are not converging to any one   C
C                 vector - the time step is probably too large.      C
C                 If ITMX is returned = -2 then the maximum number   C
C                 of iterations were exceeded.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DMAT      : Diffusion matrix created by MFDMFM (see above)     C
C                 Dimensions are ( 3*KL+1, NR*NMH )                  C
C     VI        : Velocity at time step I. (time t). Dim ( * )       C
C                 aleast NVH*NR. (Input vector)                      C
C     VIP1      : Velocity at time step I+1. (time t + deltat).      C
C                 Dim ( * )  aleast NVH*NR. (Input vector)           C
C     BI        : Magnetic field at time step I. Dim ( * )           C
C                 aleast NMH*NR. (Input vector)                      C
C     BIP1      : Magnetic field at time step I+1. (time t + deltat) C
C                 Dim ( * )  aleast NMH*NR. (Output vector)          C
C     ATI       : Advection term at TIME t (time step i)             C
C     ATIP1     : Advection term at TIME t+deltat (time step i+1)    C
C                 Both ATI and ATIP1 must have Dim aleast NMH*NR.    C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C                                                                    C
C     GAUW      : Gauss weights computed bu GAUWTS. Dim ( NTHP )     C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C                                                                    C
C     DPA       : Derivatives of the above.                          C
C       PA and DPA must be formed in advance by a call to SCHNLA.    C
C                                                                    C
C     RQST1     : Dim. ( LH*(LH+2) ,3, NR ). Work array.             C
C     RQST2     : Dim. ( LH*(LH+2) ,3, NR ). Work array.             C
C     ZCF       : Dim ( NR ). Work array.                            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     FDCM      : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM ).                   C
C                  Formed by a call to FDCM with                     C
C                  NLMN = ILN                                        C
C                  NRMN = IRN                                        C
C                  NLMC = ILN                                        C
C                  NRMC = IRN                                        C
C                                                                    C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     V1        : Dim ( NPHP, NTHP, 3 ). Work array.                 C
C     V2        : Dim ( NPHP, NTHP, 3 ). Work array.                 C
C     V3        : Dim ( NPHP, NTHP, 3 ). Work array.                 C
C                                                                    C
C     FTF1      : Dim ( 2*NPHP ). Work array.                        C
C     FTF2      : Dim ( 2*NPHP ). Work array.                        C
C     FTF3      : Dim ( 2*NPHP ). Work array.                        C
C                                                                    C
C     PARS      : Dim ( * ). An array containing the double          C
C                 precision parameters of our time step.             C
C                 Array elements are as follows:-                    C
C                                                                    C
C          PARS( 1 ) = DELTAT    : Size of time step.                C
C          PARS( 2 ) = CFAC      : See above. Weights current        C
C                                  versus subsequent time step when  C
C                                  evaluating derivatives.           C
C          PARS( 3 ) = RM        : Magnetic Reynolds number.         C
C                                  Multiplies the advective term.    C
C          PARS( 4 ) = RTOL      : The convergence tolerance         C
C                                 for the reiteration of BIP1.       C
C                                 If the difference between the      C
C                                 norms of the j^{th} and (j+1)^{th} C
C                                 iterations of BIP1 is less than    C
C                                 RTOL, we judge our iteration to be C
C                                 successful.                        C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     OATI      : Set to .TRUE. if ATI is already evaluated.         C
C                 Set to .FALSE. to zero and fill ATI.               C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SIMFTS( NR, NMH, KL, ILN, IRN, INVA, MT, ML, MM, MP,
     1     VT, VL, VM, VP, LH, M0, MMAX, NTHP, NPHP, NDRVS, NDRVM,
     2     NFDCM, NDCS, IRES, NBN, IPIV, MHBCS, ITMX,
     3     DMAT, VI, VIP1, BI, BIP1, ATI, ATIP1, GAUX, GAUW, PA, DPA,
     4     RQST1, RQST2, ZCF, SVFDC, FDCM, XARR, V1, V2, V3, FTF1,
     5     FTF2, FTF3, OATI, PARS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NMH, KL, ILN, IRN, INVA( * ), MT( * ), ML( * ),
     1        MM( * ), MP( * ), VT( * ), VL( * ), VM( * ), VP( * ),
     2        LH, M0, MMAX, NTHP, NPHP, NDRVS, NDRVM,
     3        NFDCM, NDCS, IRES, NBN, IPIV( * ), MHBCS( NDCS ), ITMX
      DOUBLE PRECISION DMAT( 3*KL + 1, NR*NMH ), VI( * ), VIP1( * ),
     1                 BI( * ), BIP1( * ), ATI( * ), ATIP1( * )
      DOUBLE PRECISION GAUX( NTHP ), GAUW( NTHP ),
     1                 XARR( NR ), ZCF( NR ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
      DOUBLE PRECISION RQST1( LH*(LH+2) ,3 ,NR ),
     1                 RQST2( LH*(LH+2) ,3 ,NR ),
     2                 SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     3                 FDCM( NFDCM, NR, NDRVM ), PARS( * )
      DOUBLE PRECISION V1( NPHP, NTHP, 3), FTF1( 2*NPHP ),
     1                 V2( NPHP, NTHP, 3), FTF2( 2*NPHP ),
     2                 V3( NPHP, NTHP, 3), FTF3( 2*NPHP )
      LOGICAL OATI
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER N2, IOP, INC, NBAND, INFO, NRHS, NOIT, INMA( 3 )
      DOUBLE PRECISION ZERO, FAC, DELTAT, CFAC, RM, DC, ABIP1,
     1                 RTOL, AOLD, DIFF, ODIFF, DNRM2
      CHARACTER *(1) TRANS
      PARAMETER ( ZERO = 0.0d0, DC = 1.0d0, INC = 1, TRANS = 'N',
     1            NRHS = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      INMA( 1 ) = 4
      INMA( 2 ) = NR
      INMA( 3 ) = NMH
C
C Fill in double precision parameters.
C
      DELTAT  = PARS( 1 )
      CFAC    = PARS( 2 )
      RM      = PARS( 3 )
      RTOL    = PARS( 4 )
C
      NBAND   = 3*KL + 1
C
      ODIFF = 1.0d8
C
C VI and BI are vectors which respectively contain
C the velocity and magnetic field at time step 'i'.
C Let ATI be the vector which contains the Advective Term
C at time step 'i' - oati tells us whether we need to evaluate
C it or not.
C
      IF ( OATI ) GOTO 49
C
C OK - we need to evaluate ATI. First zero it.
C
      IOP = 0
      N2  = NMH*NR
      CALL VECOP ( ATI, ZERO, N2, IOP )
C
C Now add R_m multiplied by \curl ( v^{i} x B^{i} ) to ATI
C
      FAC = RM
      CALL CVCBTA( NR, NMH, BI, INVA, VI, ATI, MT, ML, MM, MP,
     1          VT, VL, VM, VP, LH, M0, MMAX, ILN, IRN, NTHP, NPHP,
     2          NDRVS, NDRVM, NFDCM, NDCS, IRES, NBN, GAUX, GAUW,
     3          PA, DPA, RQST1, RQST2, ZCF, SVFDC, FDCM, XARR,
     4          V1, V2, V3, FTF1, FTF2, FTF3, FAC )
C
 49   CONTINUE
C
C the vector ATI now contains current advective term
C We will now solve for BIP1 as a 'predictor' solution
C using the guess from this rough estimate.
C First zero BIP1 ...
C
      IOP = 0
      N2  = NMH*NR
      CALL VECOP ( BIP1, ZERO, N2, IOP )
C
C Now add on the non-linear part of the forcing term
C
      FAC = DELTAT
      CALL DAXPY( N2, FAC, ATI, INC, BIP1, INC )
C
C Now add on diffusion parts ...
C
      CALL MFDRHF( NR, NMH, NDCS, MT, ML, MM, MP, NBN,
     1             NDRVS, NDRVM, NFDCM, DC, CFAC, DELTAT, BI,
     2             BIP1, SVFDC, XARR )
C
C Ok. our right hand side is now complete.
C Now let's solve the system so that BIP1 is returned
C containing the magnetic field predictor.
C
      CALL DGBTRS( TRANS, N2, KL, KL, NRHS, DMAT, NBAND,
     1             IPIV, BIP1, N2, INFO )
C
C Check solution has passed without incident.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine SIMFTS.'
        PRINT *,' The LAPACK subroutine DGBTRS'
        PRINT *,' has returned ',INFO,' as a value of'
        PRINT *,' INFO solving DMAT.x = BIP1.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Solved for BIP1 without problem.
C We 'complete' the vector here in order to
C calculate a norm.
C
      CALL ASVCPL( BIP1, NR, NDCS, INMA, MP, MHBCS, MHBCS,
     1             NFDCM, NDRVS, NDRVM, NBN, SVFDC )
C
      ABIP1 = DNRM2( N2, BIP1, INC )
C
C ABIP1 now contains the norm of the estimate
C Now, begin to iterate to a semi-implicit solution for BIP1
C
      NOIT = 0
 50   CONTINUE
      NOIT = NOIT + 1
      IF ( NOIT.GT.ITMX ) THEN
        ITMX = -2
        RETURN
      ENDIF
C
C We now calculate the advective term at time = t + deltat
C using our new guess for BIP1
C Zero ATIP1
C
      IOP = 0
      CALL VECOP ( ATIP1, ZERO, N2, IOP )
C
C Now add R_m multiplied by \curl (v^{i+1} x B^{i+1}) to ATIP1
C
      FAC = RM
      CALL CVCBTA( NR, NMH, BIP1, INVA, VIP1, ATIP1, MT, ML, MM, MP,
     1          VT, VL, VM, VP, LH, M0, MMAX, ILN, IRN, NTHP, NPHP,
     2          NDRVS, NDRVM, NFDCM, NDCS, IRES, NBN, GAUX, GAUW,
     3          PA, DPA, RQST1, RQST2, ZCF, SVFDC, FDCM, XARR,
     4          V1, V2, V3, FTF1, FTF2, FTF3, FAC )
C
C the vector ATIP1 now contains advective term
C for step i+1. We will now solve for BIP1 as a
C 'corrector' solution using our refined estimate
C First zero BIP1 ...
C
      IOP = 0
      CALL VECOP ( BIP1, ZERO, N2, IOP )
C
C Now add on advective part of forcing term from
C time step i
C
      FAC = CFAC*DELTAT
      CALL DAXPY( N2, FAC, ATI, INC, BIP1, INC )
C
C Now add on advective part of forcing term from
C time step i+1
C
      FAC = (1.0d0 - CFAC)*DELTAT
      CALL DAXPY( N2, FAC, ATIP1, INC, BIP1, INC )
C
C Now add on diffusion parts ...
C
      CALL MFDRHF( NR, NMH, NDCS, MT, ML, MM, MP, NBN,
     1             NDRVS, NDRVM, NFDCM, DC, CFAC, DELTAT, BI,
     2             BIP1, SVFDC, XARR )
C
C Ok. our right hand side is now complete.
C Now let's solve the system so that BIP1 is returned
C containing the magnetic field corrector.
C
      CALL DGBTRS( TRANS, N2, KL, KL, NRHS, DMAT, NBAND,
     1             IPIV, BIP1, N2, INFO )
C
C Check solution has passed without incident.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine SIMFTS.'
        PRINT *,' The LAPACK subroutine DGBTRS'
        PRINT *,' has returned ',INFO,' as a value of'
        PRINT *,' INFO solving DMAT.x = BIP1.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Solved for BIP1 without problem.
C We 'complete' the vector here in order to
C calculate a norm.
C
      CALL ASVCPL( BIP1, NR, NDCS, INMA, MP, MHBCS, MHBCS,
     1             NFDCM, NDRVS, NDRVM, NBN, SVFDC )
C
      AOLD = ABIP1
      ABIP1 = DNRM2( N2, BIP1, INC )
      DIFF = ABS( ABIP1 - AOLD )
      IF ( DIFF.LT.RTOL ) THEN
        ITMX = NOIT
        RETURN
      ENDIF
C
C ok - so our reiteration wasn't sufficiently
C close to return - check to see if we are getting worse
C
      IF ( DIFF.GT.ODIFF ) THEN
        ITMX = -1
        RETURN
      ENDIF
      ODIFF = DIFF
C
      GOTO 50
      END
C*********************************************************************
