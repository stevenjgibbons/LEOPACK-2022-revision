C*********************************************************************
C subroutine Induction Equation Non Linear Terms form  optimised 2 ***
C            -         -        -   -      -                     - ***
C Steve Gibbons Wed Feb 21 11:26:16 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Fills vectors R4 and R5 with the following terms in the induction  C
C equation:    CM curl ( v x B ) ...                                 C
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
C     NTHP      : Number of theta points                             C
C     NPHP      : Number of phi points                               C
C     NCMX      : Leading dimension of XSV array.                    C
C                                                                    C
C     ICMRAD    : Location in XSV of v_rad.                          C
C     ICMTHE    : Location in XSV of v_the.                          C
C     ICMPHI    : Location in XSV of v_phi.                          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV4       : Dim ( NH4*NR ). Poloidal field initial soln.       C
C     SV5       : Dim ( NH5*NR ). Toroidal field initial soln.       C
C                                                                    C
C     R4        : Dim ( NH4*NR ). Work array.                        C
C     R5        : Dim ( NH5*NR ). Work array.                        C
C                                                                    C
C     SV4Q      : Dim (         1, NH4*NR ) turns p. mag. to Q( r )  C
C     SV4S      : Dim ( 2*NBN + 1, NH4*NR ) turns p. mag. to S( r )  C
C     SV5T      : Dim (         1, NH5*NR ) turns t. mag. to T( r )  C
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
      SUBROUTINE IENLT2( SV4,SV5,R4,R5,F1,F2,F3,XSV,SV4Q,
     1                   SV4S,SV5T,VQ4,VS4,VT5,VT4,VQ5,VS5,CT4P,
     2                   CQ5T,CS5T,ICMRAD,ICMTHE,ICMPHI,FTFPM,
     3                   FTFTM,INFPM,INFTM,PFAM,TFAM,JPFAM,JTFAM)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Common block contents .....................C
      INTEGER          NR, NH4, NH5, NBN, M0,
     1                 IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX
      COMMON  /NDIMPARS/ NR, NH4, NH5, NBN, M0,
     1                   IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX
      DOUBLE PRECISION CFAC, DELTAT, DTOL, CM
      COMMON  /DPHYSPARS/ CFAC, DELTAT, DTOL, CM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
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
C Forcing term vectors
C
      DOUBLE PRECISION R4( NH4*NR ), R5( NH5*NR )
C
C Work arrays
C
      DOUBLE PRECISION XSV( NCMX, NPHP, NTHP, NR)
      DOUBLE PRECISION VQ4( NH4*NR ), VS4( NH4*NR ), VT5( NH5*NR ),
     1                 VT4( NH4*NR ), VQ5( NH5*NR ), VS5( NH5*NR )
C
C Arrays necessary for transform
C
      DOUBLE PRECISION F1(2*NPHP), F2(2*NPHP), F3(2*NPHP)
C
C Transform coefficients
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
      INTEGER          ICMRAD, ICMTHE, ICMPHI
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          K, N, INCX, INCY, LDA, M, ILNR, IRNR,
     1                 ICMR, ICMT, ICMP, IFORMF
      CHARACTER *(1)   TRANS
      LOGICAL          LW
      DOUBLE PRECISION ALPHA, BETA, ZERO
C
      PARAMETER ( TRANS = 'N', INCX   = 1, INCY   = 1,
     1            ZERO = 0.0d0, IFORMF = 4 )
C____________________________________________________________________C
C
      IF ( IOUTF.EQ.0 ) THEN
        LW  = .FALSE.
      ELSE
        LW  = .TRUE.
      ENDIF
C
C Now calculate scaloidal part of magnetic field.
C Use DGBMV to multiply SV4 by matrix SV4Q to give VQ4
C
      ALPHA  = 1.0d0
      BETA   = ZERO
      LDA    = 1
      M      = NR*NH4
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, SV4Q, LDA, SV4, INCX,
     1             BETA, VQ4, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'IENLT2: Called DGBMV: vq4:= sv4q . sv4.'
C
C Now calculate spheroidal part of magnetic field.
C Use DGBMV to multiply SV4 by matrix SV4S to give VS4
C
      ALPHA  = 1.0d0
      BETA   = ZERO
      LDA    = 1+2*NBN
      M      = NR*NH4
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, SV4S, LDA, SV4, INCX,
     1             BETA, VS4, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'IENLT2: Called DGBMV: vs4:= sv4s . sv4.'
C
C Now calculate toroidal part of magnetic field.
C Use DGBMV to multiply SV5 by matrix SV5T to give VT5
C
      ALPHA  = 1.0d0
      BETA   = ZERO
      LDA    = 1
      M      = NR*NH5
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, SV5T, LDA, SV5, INCX,
     1             BETA, VT5, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'IENLT2: Called DGBMV: vt5:= sv5t . sv5.'
C
C we now have the mag. field in QST format and so can transform
C into real space.
C
      ILNR = 2
      IRNR = NR - 1
C     .  put radial component in icmr = 4 
      ICMR = 4
C     .  put theta component in icmt = 5 
      ICMT = 5
C     .  put phi component in icmp = 6 
      ICMP = 6
      CALL RSDV2D( NTHP, NPHP,     NR, ILNR, IRNR, NH4,
     1                  NH5,           NCMX, ICMR, ICMT, ICMP,
     2             VQ4, VS4, VT5, XSV, F1, F2, F3,
     3             FTFPM, FTFTM, INFPM, INFTM, IFORMF )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'IENLT2: Called RSDV2D: XSV now contains mag. field.'
C
C Magnetic field is now stored in XSV components 4-6.
C Now evaluate all of the non-linear operations
C within XSV. This is now a single subroutine call!
C
      CALL XSVCBE( NTHP, NPHP, NR, ILNR, IRNR, NCMX, ICMRAD,
     1             ICMTHE, ICMPHI, XSV )
C
C Transform the non-linear terms for the
C induction equation back into QST space.
C Note that we will use VQ5, VS5 and VT4
C
      ILNR = 2
      IRNR = NR - 1
C     .  radial component is stored in icmr = 1
      ICMR = 1
C     .  theta component is stored in icmt = 2
      ICMT = 2
C     .  phi component is stored in icmp = 3
      ICMP = 3
      CALL XSVSDD( NTHP, NPHP,     NR, ILNR, IRNR, NH5,
     1                  NH4,           NCMX, ICMR, ICMT, ICMP,
     2             VQ5, VS5, VT4, XSV, F1, F2, F3,
     3             PFAM, TFAM, JPFAM, JTFAM, IFORMF )
C
C Now add the curl of the scaloidal function to the
C toroidal forcing term of the magnetic field (R5)
C
      ALPHA  = CM
      BETA   = ZERO
      LDA    = 1
      M      = NR*NH5
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CQ5T, LDA, VQ5, INCX,
     1             BETA, R5, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'IENLT2: Called DGBMV: r5:= CM*cq5t . vq5.'
C
C Now add the curl of the spheroidal function to the
C toroidal forcing term of the magnetic field (R5)
C
      ALPHA  = CM
      BETA   = 1.0d0
      LDA    = 1+2*NBN
      M      = NR*NH5
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CS5T, LDA, VS5, INCX,
     1             BETA, R5, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'IENLT2: Called DGBMV: r5:= r5 + CM*cs5t . vs5.'
C
C Now add the curl of the toroidal function to the
C poloidal forcing term of the magnetic field (R4)
C
      ALPHA  = CM
      BETA   = ZERO
      LDA    = 1
      M      = NR*NH4
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CT4P, LDA, VT4, INCX,
     1             BETA, R4, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'IENLT2: Called DGBMV: r4:= CM*ct4p . vt4.'
C
      RETURN
      END
C*********************************************************************
