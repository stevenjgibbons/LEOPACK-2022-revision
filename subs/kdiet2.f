C*********************************************************************
C subroutine Kinematic Dynamo Induction Equation Time step (optim 2) *
C            -         -      -         -        -                -  *
C Steve Gibbons Wed Feb 21 13:39:35 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Progress the induction equation by one time-step. The velocity is  C
C already stored in XSV form for both time-step (i) and time-step    C
C (i+1). The indices for the location of the velocity is stored in   C
C the COMMON block integers ICRAD0, ICTHE0, ICPHI0 (for time-step i) C
C and ICRAD1, ICTHE1, ICPHI1 (for time-step i+1).                    C
C                                                                    C
C  The induction equation is                                         C
C                                                                    C
C    CK d B / dt      =   CL Lap B                                   C
C                       + CM curl ( v x B )                          C
C                                                                    C
C This routine only needs to know the value of CM as the matrices    C
C used to calculate the implicit part of the time-integration are    C
C precalculated using the sub-routines PFTSMF and TFTSMF.            C
C                                                                    C
C Several matrices must be pre-calculated in order to take curls     C
C etc. of the magnetic field functions. They are as follows          C
C                                                                    C
C  SV4Q:  Dimension ( 1, NH4*NR )                                    C
C  ====                                                              C
C         Will multiply the solution vector SV4 to give a vector     C
C         of scaloidal radial functions. Form it using the following C
C         piece of code:-                                            C
C                                                                    C
c     N1     = 1
c     N2     = NH4*NR
c     IOPT   = 1
c     K      = 0
c     ILN    = 2
c     IRN    = NR-1
c     CALL VOBMBR( N1, N2, K, NR, NH4, NBN, NCFM, ILN, IRN,
c    1             ML4, IOPT, XARR, SV4Q, FDCM )
c     WRITE ( LULOG, * ) 'Formed SV4Q.'
C
C  SV4S:  Dimension ( 1+2*NBN, NH4*NR )                              C
C  ====                                                              C
C         Will multiply the solution vector SV4 to give a vector of  C
C         spheroidal radial functions. Form it using the following   C
C         piece of code:-                                            C
c     N1     = 1+2*NBN
c     N2     = NH4*NR
c     IOPT   = 2
c     ILN    = 2
c     IRN    = NR-1
c     K      = NBN
c     CALL AVBMBR( N1, N2, K, NR, NH4, NBN, NCFM, NDRVM, ILN,
c    1             IRN, ML4, MP4, IOPT, NDCS, XARR, SV4S, SVFDC )
c     WRITE ( LULOG, * ) 'Formed SV4S.'
C
C  SV5T:  Dimension ( 1, NH5*NR )                                    C
C  ====                                                              C
C         Will multiply the solution vector SV5 to give a vector of  C
C         toroidal radial functions. Form it using the following     C
C         piece of code:-                                            C
c     N1     = 1
c     N2     = NH5*NR
c     IOPT   = 3
c     K      = 0
c     ILN    = 2
c     IRN    = NR-1
c     CALL VOBMBR( N1, N2, K, NR, NH5, NBN, NCFM, ILN, IRN,
c    1             ML5, IOPT, XARR, SV5T, FDCM )
c     WRITE ( LULOG, * ) 'Formed SV5T.'
C                                                                    C
C  CQ5T:  Dimension ( 1, NH5*NR )                                    C
C  ====                                                              C
C         Will multiply the vector VQ5 to give toroidal forcing term C
C         R5. Form it using the following piece of code:-            C
c     N1     = 1
c     N2     = NH5*NR
c     IOPT   = 11
c     K      = 0
c     ILN    = 2
c     IRN    = NR-1
c     CALL VOBMBR( N1, N2, K, NR, NH5, NBN, NCFM, ILN, IRN,
c    1             ML5, IOPT, XARR, CQ5T, FDCM )
c     WRITE ( LULOG, * ) 'Formed CQ5T.'
C                                                                    C
C  CS5T:  Dimension ( 1+2*NBN, NH5*NR )                              C
C  ====                                                              C
C         Will multiply the vector VS5 to give toroidal forcing term C
C         R5. Form it using the following piece of code:-            C
c     N1     = 1+2*NBN
c     N2     = NH5*NR
c     IOPT   = 12
c     K      = NBN
c     ILN    = 2
c     IRN    = NR-1
c     CALL VOBMBR( N1, N2, K, NR, NH5, NBN, NCFM, ILN, IRN,
c    1             ML5, IOPT, XARR, CS5T, FDCM )
c     WRITE ( LULOG, * ) 'Formed CS5T.'
C                                                                    C
C  CT4P:  Dimension ( 1, NH4*NR )                                    C
C  ====                                                              C
C         Will multiply the vector VT4 to give poloidal forcing term C
C         R4. Form it using the following piece of code:-            C
c     N1     = 1
c     N2     = NH4*NR
c     IOPT   = 13
c     K      = 0
c     ILN    = 2
c     IRN    = NR-1
c     CALL VOBMBR( N1, N2, K, NR, NH4, NBN, NCFM, ILN, IRN,
c    1             ML4, IOPT, XARR, CT4P, FDCM )
c     WRITE ( LULOG, * ) 'Formed CT4P.'
C                                                                    C
C The array FDCM is not called by KDIET2 directly, but must be       C
C for the above auxiliary routines. Use the following code to form   C
C FDCM:-       [dimensions are (NCFM, NR, 1 )]                       C
C               for NCFM = 2*NBN+1                                   C
c     NDRVS = 1
c     ILN   = 2
c     IRN   = NR - 1
c     CALL FDCMBD( NR, NBN, ILN, IRN, ILN, IRN, NCFM,
c    1             NCFM, NDRVS, NDRVS, IWORK, XARR, FDCM,
c    2             COEFM1, W1, W2 )
C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Common Block                                                      C
C  ------------                                                      C
C                                                                    C
C   /NDIMPARS/                                                       C
C                                                                    C
C  Contains the integer numbers NR, NH4, NH5, NBN,                   C
C                    M0, IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX     C
C                                                                    C
C   /DPHYSPARS/                                                      C
C                                                                    C
C  Physical parameters for problem. All double precision scalars     C
C                                                                    C
C   CFAC, DELTAT, DTOL, CM                                           C
C                                                                    C
C     CFAC      : Determines how implicit/explicit integration is.   C
C                 Must be strictly greater than 0.0 and strictly     C
C                 less than 1.0  The higher CFAC is, the more        C
C                 explicit the integration. CFAC = 0.5 --> Crank-    C
C                 Nicolson scheme.                                   C
C     DELTAT    : Time-step size.                                    C
C     DTOL      : How close the solution norms on consecutive        C
C                 iterations need to be.                             C
C                                                                    C
C  where                                                             C
C                                                                    C
C     NR        : Total number of radial grid nodes.                 C
C     NH4       : Number of poloidal field harmonics in soln.        C
C     NH5       : Number of toroidal field harmonics in soln.        C
C     NBN       : Number of bounding nodes for finite differences.   C
C     M0        : Lowest non-zero wavenumber in solution.            C
C     IOUTF     : Logical Unit number of output info file.           C
C                 (If IOUTF = 0 then no output is written).          C
C     LH        : Highest spherical harmonic degree, l.              C
C     ITMX      : Maximum number of iterations allowed per           C
C                  time-step                                         C
C     NOIT      : Actual number of iterations taken.                 C
C                  If number of iterations required exceeds ITMX     C
C                  then NOIT is returned as -1.                      C
C                  If soln. norm appears to be increasing after      C
C                  3 iterations, NOIT is returned as -2.             C
C     NTHP      : Number of theta points                             C
C     NPHP      : Number of phi points                               C
C     NCMX      : Leading dimension of XSV array.                    C
C                                                                    C
C   /ICOMPS/                                                         C
C                                                                    C
C     ICRAD0, ICTHE0, ICPHI0, ICRAD1, ICTHE1, ICPHI1                 C
C                                                                    C
C   ICRAD0 is the location in XSV for v_rad at time-step i.          C
C   ICTHE0 is the location in XSV for v_the at time-step i.          C
C   ICPHI0 is the location in XSV for v_phi at time-step i.          C
C                                                                    C
C   ICRAD1 is the location in XSV for v_rad at time-step i+1.        C
C   ICTHE1 is the location in XSV for v_the at time-step i+1.        C
C   ICPHI1 is the location in XSV for v_phi at time-step i+1.        C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IP4       : Dim (NH4*NR). Pivotting information for AM4        C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C     IP5       : Dim (NH5*NR). Pivotting information for AM5        C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     AM4       : Matrix for solution of f^{i+1}. (poloidal m.f.)    C
C                 Has dimensions ( 3*NBN + 1, NH4*NR )               C
C                                                                    C
C     BM4       : Matrix for multiplication of f^i. (poloidal m.f.)  C
C                 Has dimensions ( 2*NBN + 1, NH4*NR )               C
C                                                                    C
C     AM5       : Matrix for solution of f^{i+1}. (toroidal m.f.)    C
C                 Has dimensions ( 3*NBN + 1, NH5*NR )               C
C                                                                    C
C     BM5       : Matrix for multiplication of f^i. (toroidal m.f.)  C
C                 Has dimensions ( 2*NBN + 1, NH5*NR )               C
C                                                                    C
C     SV4       : Dim ( NH4*NR ). Poloidal field initial soln.       C
C     SV5       : Dim ( NH5*NR ). Toroidal field initial soln.       C
C                                                                    C
C     E4        : Dim ( NH4*NR ). Poloidal field final soln.         C
C     E5        : Dim ( NH5*NR ). Toroidal field final soln.         C
C                                                                    C
C     D4        : Dim ( NH4*NR ). Work array.                        C
C     D5        : Dim ( NH5*NR ). Work array.                        C
C                                                                    C
C     R4        : Dim ( NH4*NR ). Work array.                        C
C     R5        : Dim ( NH5*NR ). Work array.                        C
C                                                                    C
C     SV4Q      : Dim (         1, NH4*NR ) turns p. mag. to Q( r )  C
C     SV4S      : Dim ( 2*NBN + 1, NH4*NR ) turns p. mag. to S( r )  C
C     SV5T      : Dim (         1, NH5*NR ) turns t. mag. to T( r )  C
C                                                                    C
C     CQ5T      : Dim (         1, NH5*NR ) takes curl of Q( r )     C
C     CS5T      : Dim ( 2*NBN + 1, NH5*NR ) takes curl of S( r )     C
C     CT4P      : Dim (         1, NH4*NR ) takes curl of T( r )     C
C                                                                    C
C     F1        : Dim (2*NPHP). Work array for fourier transforming. C
C     F2        : Dim (2*NPHP). Work array for fourier transforming. C
C     F3        : Dim (2*NPHP). Work array for fourier transforming. C
C                                                                    C
C     VQ4       : Dim ( NH4*NR ). Work array.                        C
C     VS4       : Dim ( NH4*NR ). Work array.                        C
C     VT5       : Dim ( NH5*NR ). Work array.                        C
C     VT4       : Dim ( NH4*NR ). Work array.                        C
C     VQ5       : Dim ( NH5*NR ). Work array.                        C
C     VS5       : Dim ( NH5*NR ). Work array.                        C
C                                                                    C
C     XSV       : Work array. Dim ( NCMX, NPHP, NTHP, NR )           C
C                                                                    C
C     FTFPM     : Dim (3,NH4,NTHP). Coeff.s from RSDV2C              C
C     FTFTM     : Dim (2,NH5,NTHP). Coeff.s from RSDV2C              C
C     INFPM     : Integer! Dim (2, NH4 ) Indices from RSDV2C         C
C     INFTM     : Integer! Dim (2, NH5 ) Indices from RSDV2C         C
C                                                                    C
C To create the above 4 arrays, make the following call:             C
C                                                                    C
C     CALL RSDV2C( NTHP, M0, LH, NH4, ML4, MM4, NH5, ML5, MM5,       C
C    1             GAUX, PA, DPA, FTFPM, FTFTM, INFPM, INFTM )       C
C                                                                    C
C     PFAM      : Dim ( 3, NH5, NTHP ). Coeffs from XSVSDC.          C
C     TFAM      : Dim ( 2, NH4, NTHP ). Coeffs from XSVSDC.          C
C     JPFAM     : Dim(2,NH5). (Integers!) Locations from XSVSDC      C
C     JTFAM     : Dim(2,NH4). (Integers!) Locations from XSVSDC      C
C                                                                    C
C To create the above 4 arrays, make the following call:             C
C                                                                    C
C     CALL XSVSDC( NTHP, M0, LH, NH5, ML5, MM5, NH4, ML4, MM4,       C
C                  GAUX, GAUW, PA, DPA, PFAM, TFAM, JPFAM, JTFAM )   C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE KDIET2(     IP4,    IP5,AM4,BM4,AM5,BM5,
     1             SV4,SV5,D4,D5,E4,E5,R4,R5,F1,F2,F3,XSV,SV4Q,SV4S,
     2             SV5T,VQ4,VS4,VT5,VT4,VQ5,VS5,CT4P,CQ5T,CS5T,
     3             FTFPM,FTFTM,INFPM,INFTM,PFAM,TFAM,JPFAM,JTFAM)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Common block contents .....................C
      INTEGER          NR, NH4, NH5, NBN, M0,
     1                 IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX
      COMMON  /NDIMPARS/ NR, NH4, NH5, NBN, M0,
     1                   IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX
      INTEGER          ICRAD0, ICTHE0, ICPHI0, ICRAD1, ICTHE1, ICPHI1
      COMMON  /ICOMPS/ ICRAD0, ICTHE0, ICPHI0, ICRAD1, ICTHE1, ICPHI1
      DOUBLE PRECISION CFAC, DELTAT, DTOL, CM
      COMMON  /DPHYSPARS/ CFAC, DELTAT, DTOL, CM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
C Basic solution vector indices
C
      INTEGER                      IP4( NH4*NR )
      INTEGER                      IP5( NH5*NR )
C
C Diffusion/time-step matrices
C
      DOUBLE PRECISION AM4( 3*NBN+1, NH4*NR ), BM4( 2*NBN+1, NH4*NR )
      DOUBLE PRECISION AM5( 3*NBN+1, NH5*NR ), BM5( 2*NBN+1, NH5*NR )
C
C Auxiliary matrices
C
      DOUBLE PRECISION SV4Q(         1, NH4*NR )
      DOUBLE PRECISION SV4S( 2*NBN + 1, NH4*NR )
      DOUBLE PRECISION SV5T(         1, NH5*NR )
C
      DOUBLE PRECISION CQ5T(       1, NH5*NR )
      DOUBLE PRECISION CS5T( 2*NBN+1, NH5*NR )
      DOUBLE PRECISION CT4P(       1, NH4*NR )
C
C Initial solution vectors
C
      DOUBLE PRECISION SV4( NH4*NR ), SV5( NH5*NR )
C
C Final solution vectors
C
      DOUBLE PRECISION E4( NH4*NR ), E5( NH5*NR )
C
C Arrays necessary for transform
C
      DOUBLE PRECISION F1(2*NPHP), F2(2*NPHP), F3(2*NPHP)
C
C (pre-calculated by RSDV2C)
C
      INTEGER          INFPM( 2, NH4 ), INFTM( 2, NH5 )
C
      DOUBLE PRECISION FTFPM( 3, NH4, NTHP ),
     1                 FTFTM( 2, NH5, NTHP )
C
C (pre-calculated by XSVSDC)
C
      INTEGER          JPFAM( 2, NH5 ), JTFAM( 2, NH4 )
C
      DOUBLE PRECISION PFAM( 3, NH5, NTHP ),
     1                 TFAM( 2, NH4, NTHP )
C
C Work arrays
C
      DOUBLE PRECISION D4( NH4*NR ), D5( NH5*NR )
      DOUBLE PRECISION R4( NH4*NR ), R5( NH5*NR )
      DOUBLE PRECISION VQ4( NH4*NR ), VS4( NH4*NR ), VT5( NH5*NR ),
     1                 VT4( NH4*NR ), VQ5( NH5*NR ), VS5( NH5*NR ),
     2                 XSV( NCMX, NPHP, NTHP, NR)
C
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          KL, KU, LDA, M, N, INCX, INCY, NRHS, INFO,
     1                 LLU, ICMRAD, ICMTHE, ICMPHI
      LOGICAL          LW
      CHARACTER *(1)   TRANS
      DOUBLE PRECISION ALPHA, BETA, ODIFF, DIFF, AVEC1, AOLD,
     1                 DNRM2, FAC
      EXTERNAL         DNRM2
C
      PARAMETER ( TRANS = 'N', INCX   = 1, INCY   = 1, NRHS = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Odiff should begin very large
C
      ODIFF = 1.0d8
C
      KL = NBN
      KU = NBN
C
      IF ( IOUTF.EQ.0 ) THEN
        LW  = .FALSE.
        LLU = 6
      ELSE
        LW  = .TRUE.
        LLU = IOUTF
      ENDIF
C
      IF ( LW ) WRITE ( IOUTF, * ) 'Entered KDIET2.'
C
C Calculate the contribution from the diffusive parts.
C
C Poloidal magnetic field:
C   Multiply banded matrix BM4 by vector SV4 to give D4
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = KL + KU + 1
      M      = NR*NH4
      N      = M
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, BM4, LDA, SV4, INCX,
     1             BETA, D4, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'KDIET2: Called DGBMV: dv4:= bm4 . sv4.'
C
C Toroidal magnetic field:
C   Multiply banded matrix BM5 by vector SV5 to give D5
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = KL + KU + 1
      M      = NR*NH5
      N      = M
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, BM5, LDA, SV5, INCX,
     1             BETA, D5, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'KDIET2: Called DGBMV: dv5:= bm5 . sv5.'
C
C Copy the diffusion forcing terms into E4 and E5.
C
      N = NR*NH4
      CALL DCOPY( N, D4, INCX, E4, INCX )
C
      N = NR*NH5
      CALL DCOPY( N, D5, INCX, E5, INCX )
C
C We now need to form non-linear forcing terms for step i
C These are put into the vectors R4 and R5.
C
      ICMRAD = ICRAD0
      ICMTHE = ICTHE0
      ICMPHI = ICPHI0
C
      CALL IENLT2(                 SV4,SV5,R4,R5,       F1,F2,F3,
     1             XSV,          SV4Q,SV4S,SV5T,VQ4,VS4,VT5,VT4,VQ5,VS5,
     2             CT4P,CQ5T,CS5T,ICMRAD,ICMTHE,ICMPHI,
     3             FTFPM,FTFTM,INFPM,INFTM,PFAM,TFAM,JPFAM,JTFAM)
C
C We now add the non-linear terms (stored in R4 and R5) to E4 and E5
C
      N = NR*NH4
      CALL DAXPY( N, DELTAT, R4, INCX, E4, INCX )
C
      N = NR*NH5
      CALL DAXPY( N, DELTAT, R5, INCX, E5, INCX )
C
C Solve the system of equations for to form predictor
C We use the LAPACK routine DGBTRS
C Now, solve for poloidal magnetic field:
C
      N      = NR*NH4
      LDA    = 3*NBN + 1
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM4, LDA, IP4, E4, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine KDIET2.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for predictor solution of E4.'
      ENDIF
C
C Now, solve for toroidal magnetic field:
C
      N      = NR*NH5
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM5, LDA, IP5, E5, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine KDIET2.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for predictor solution of E5.'
      ENDIF
C
C The predictor is now stored in E4 and E5
C Calculate the norm of this solution.
C
      N     = NR*NH4
      AVEC1 = DNRM2( N, E4, INCX )
      N     = NR*NH5
      AVEC1 = AVEC1 + DNRM2( N, E5, INCX )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'KDIET2: Predictor norm = ', AVEC1
C
C We then add the step^i forcing terms (non-linear)
C to the diffusion forcing terms. This will reduce
C the number of calculations done if more than one
C iteration is required.
C
      FAC = CFAC*DELTAT
C
      N   = NR*NH4
      CALL DAXPY( N, FAC, R4, INCX, D4, INCX )
C
      N   = NR*NH5
      CALL DAXPY( N, FAC, R5, INCX, D5, INCX )
C
C Now begin the loop around corrector iterations
C
      NOIT  = 0
 50   CONTINUE
      NOIT = NOIT + 1
C
      IF ( NOIT.GT.ITMX ) THEN
        WRITE ( LLU, * ) 'Subroutine KDIET2.'
        WRITE ( LLU, * ) 'Max iterations exceeded.'
        NOIT = -1
        RETURN
      ENDIF
C
C Calculate the non-linear forcing terms in vectors R4 and R5
C
      ICMRAD = ICRAD1
      ICMTHE = ICTHE1
      ICMPHI = ICPHI1
C
      CALL IENLT2(                 E4,E5,R4,R5,       F1,F2,F3,
     1             XSV,          SV4Q,SV4S,SV5T,VQ4,VS4,VT5,VT4,VQ5,VS5,
     2             CT4P,CQ5T,CS5T,ICMRAD,ICMTHE,ICMPHI,
     3             FTFPM,FTFTM,INFPM,INFTM,PFAM,TFAM,JPFAM,JTFAM)
C
      N = NR*NH4
      CALL DCOPY( N, D4, INCX, E4, INCX )
C
      N = NR*NH5
      CALL DCOPY( N, D5, INCX, E5, INCX )
C
C Now add the step^{i+1} non-lin. terms to E4 and E5
C
      FAC = (1.0d0-CFAC)*DELTAT
C
      N   = NR*NH4
      CALL DAXPY( N, FAC, R4, INCX, E4, INCX )
C
      N   = NR*NH5
      CALL DAXPY( N, FAC, R5, INCX, E5, INCX )
C
C Solve the system of equations for iteration NOIT
C We use the LAPACK routine DGBTRS
C Now, solve for poloidal magnetic field:
C
      N      = NR*NH4
      LDA    = 3*NBN + 1
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM4, LDA, IP4, E4, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine KDIET2.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for solution ',NOIT,' of E4.'
      ENDIF
C
C Now, solve for toroidal magnetic field:
C
      N      = NR*NH5
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM5, LDA, IP5, E5, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine KDIET2.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for solution ',NOIT,' of E5.'
      ENDIF
C
      AOLD  = AVEC1
      N     = NR*NH4
      AVEC1 = DNRM2( N, E4, INCX )
      N     = NR*NH5
      AVEC1 = AVEC1 + DNRM2( N, E5, INCX )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'KDIET2: Iteration ',NOIT,' norm = ', AVEC1
C
      DIFF = DABS( AVEC1 - AOLD )
      IF ( DIFF.LT.DTOL ) THEN
        IF ( LW ) WRITE ( IOUTF, * )
     1    'KDIET2: Soln. converged iteration ', NOIT
          RETURN
      ENDIF
C
C Check to see if our norm appears to be getting bigger
C
      IF ( DIFF.GT.ODIFF .AND. NOIT.GT.3 ) THEN
        IF ( LW ) WRITE ( IOUTF, * )
     1    'KDIET2: Soln. norm increasing: iteration ', NOIT
        NOIT = -2
        RETURN
      ENDIF
C
      ODIFF = DIFF
      GOTO 50
C
      END
C*********************************************************************
