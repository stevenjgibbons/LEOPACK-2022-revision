C*********************************************************************
C subroutine Curl of Velocity Cross B Term Add ***********************
C            -       -        -     - -    -   ***********************
C Steve Gibbons Tue Nov  2 14:24:25 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C If B is a magnetic field vector containing NMH harmonics each      C
C represented at NR grid nodes, with IFORMF = 4 (see INDFUN) and     C
C V is a velocity whose details are given by INVA, then CVCBTA       C
C will add FAC (possibly a magnetic Reynolds number?) multiplied by  C
C curl ( V cross B ) to RHS (which has identical format to B.        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C     NMH       : Number of magnetic field harmonics.                C
C                                                                    C
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
C     ILN       : Leftmost node to be used.                          C
C     IRN       : Rightmost node to be used.                         C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     B         : Magnetic field. Dim ( * ) aleast NMH*NR            C
C                                                                    C
C     V         : Velocity. Dim ( * ) aleast NVH*NR                  C
C                                                                    C
C     RHS       : Returned with addition of curl ( v x B ) * FAC     C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C                                                                    C
C     GAUW      : Gauss weights computed bu GAUWTS. Dim ( NTHPTS )   C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
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
C     FAC       : Multiplicative constant.                           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CVCBTA( NR, NMH, B, INVA, V, RHS, MT, ML, MM, MP,
     1          VT, VL, VM, VP, LH, M0, MMAX, ILN, IRN, NTHP, NPHP,
     2          NDRVS, NDRVM, NFDCM, NDCS, IRES, NBN, GAUX, GAUW,
     3          PA, DPA, RQST1, RQST2, ZCF, SVFDC, FDCM, XARR,
     4          V1, V2, V3, FTF1, FTF2, FTF3, FAC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NMH, INVA( * ), MT( * ), ML( * ), MM( * ), MP( * ),
     1        VT( * ), VL( * ), VM( * ), VP( * ), LH, M0, MMAX, ILN,
     2        IRN, NTHP, NPHP, NDRVS, NDRVM, NFDCM, NDCS, IRES, NBN
      DOUBLE PRECISION B( * ), V( * ), RHS( * ), GAUX( NTHP ),
     1                 GAUW( NTHP ), XARR( NR ), ZCF( NR ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
      DOUBLE PRECISION RQST1( LH*(LH+2) ,3 ,NR ),
     1                 RQST2( LH*(LH+2) ,3 ,NR ),
     2                 SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     3                 FDCM( NFDCM, NR, NDRVM ), FAC
      DOUBLE PRECISION V1( NPHP, NTHP, 3), FTF1( 2*NPHP ),
     1                 V2( NPHP, NTHP, 3), FTF2( 2*NPHP ),
     2                 V3( NPHP, NTHP, 3), FTF3( 2*NPHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION DLOW
      INTEGER INMA( 3 )
      PARAMETER ( DLOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Quick return for vanishingly small R_m
C
      IF ( ABS( FAC ).LT.DLOW ) RETURN
C
      IF ( INVA( 2 ).NE.NR ) THEN
        PRINT *,' Subroutine CVCBTA.'
        PRINT *,' NR        = ', NR
        PRINT *,' INVA( 2 ) = ', INVA( 2 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INMA( 1 ) = 4
      INMA( 2 ) = NR
      INMA( 3 ) = NMH
C
      CALL ASVCP( NDCS, NBN, IRES, INVA, VT, VL, VM, VP, INMA, MT,
     1            ML, MM, MP, NR, LH, MMAX, ILN, IRN, NTHP, NPHP,
     2            NDRVS, NDRVM, NFDCM, V, B, GAUX, GAUW, PA, DPA,
     3            RQST1, ZCF, SVFDC, XARR, 'VEL', 'MAG', V1, V2,
     4            V3, FTF1, FTF2, FTF3 )
C
C V cross B is now contained in RQST1 - now take curl
C
      CALL RQSTCL( LH, NR, M0, MMAX, NBN, NFDCM, NDRVS, ILN,
     1             IRN, ILN, IRN, RQST1, RQST2, XARR, FDCM )
C
C Curl ( V cross B ) is now in RQST2
C So add it to the appropriate terms in RHS
C
      CALL RQSTSV( NR, LH, ILN, IRN, INMA, MT, ML, MM,
     1             'MAG', RQST2, RHS, XARR, FAC )
C
      RETURN
      END
C*********************************************************************
