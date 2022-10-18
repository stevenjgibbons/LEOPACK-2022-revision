C*********************************************************************
C subroutine Adapted Solution Vector Curl of Coriolis Force **********
C            -       -        -      -       -        -     **********
C Steve Gibbons Sun Nov 14 10:06:55 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     NTHP      : Number of theta points.                            C
C                                                                    C
C     NPHP      : Number of phi points.                              C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                 LH is the maximum l which can 'safely' be          C
C                 computed on this theta/phi grid. If a harmonic in  C
C                 the solution vector has l exceeding LH then the    C
C                 action is determined by IRES.                      C
C                                                                    C
C     IRES      : The resolution flag.                               C
C                 If IRES = 1, harmonics with l greater than LH      C
C                 or m greater than MMAX are simply ignored.         C
C                 If IRES = 2, the program aborts if LH or MMAX      C
C                 is exceeded.                                       C
C                                                                    C
C     MMAX      : Maximum spherical harmonic order, m.               C
C                                                                    C
C     M0        : Lowest non-zero mode if simple multiplicity        C
C                  applies. (Put 1 if this is safest).               C
C                                                                    C
C     NBN       : Number of bounding nodes.                          C
C                                                                    C
C     NFDCM     : Leading dimension of FDCM. At least (2*NBN+1)      C
C                                                                    C
C     NR        : Number of radial grid nodes in each function.      C
C                                                                    C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     ILN       : Leftmost node to be used.                          C
C     IRN       : Rightmost node to be used.                         C
C                                                                    C
C     INIA      : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INIA( 1 ) = IFORMI - flag for vector format.      C
C                                        See INDFUN                  C
C                  INIA( 2 ) = NRI. Must be consistent with NR.      C
C                  INIA( 3 ) = NIH = total number of radial func.s   C
C                                                                    C
C     MIT       : Array length ( * ) - atleast length NIH            C
C                                                                    C
C         MIT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MIT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MIT( IH ) = 3 --> harmonic is temperature.                 C
C         MIT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MIT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MIL       : Array length ( * ) - atleast length NIH            C
C                  Sph. harm. degree, l.                             C
C                                                                    C
C     MIM       : Array length ( * ) - atleast length NIH            C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MIP       : Array length ( * ) - atleast length NIH            C
C                  Pointer array to finite difference coefficients.  C
C                  MIP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C (Note that INIA, MIT, MIL and MIP all correspond to the input      C
C vector. INOA, MOT and  MOL all contain identical information       C
C for the output vector.)                                            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VI        : Input solution vector. Dim ( * ) but length        C
C                 atleast  NR*NIH.                                   C
C                                                                    C
C     VO        : Output solution vector. Dim ( * ) but length       C
C                 atleast  NR*NOH.                                   C
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
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     V1        : Dim ( NPHP, NTHP, 3 ). Work array.                 C
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
      SUBROUTINE ASVCCF( NDCS, NTHP, NPHP, LH, IRES, MMAX, M0, NBN, 
     1          NFDCM, NR, NDRVS, NDRVM, ILN, IRN, INIA, MIT, MIL,
     2          MIM, MIP, INOA, MOT, MOL, MOM, VI, VO, SVFDC, FDCM,
     3          GAUX, GAUW, PA, DPA, RQST1, RQST2, ZCF, XARR, V1,
     4          FTF1, FTF2, FTF3, FAC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDCS, NTHP, NPHP, LH, IRES, MMAX, M0, NBN, NFDCM, NR,
     1        NDRVS, NDRVM, ILN, IRN, INIA( * ), MIT( * ), MIL( * ),
     2        MIM( * ), MIP( * ), INOA( * ), MOT( * ), MOL( * ),
     3        MOM( * )
      DOUBLE PRECISION VI( * ), VO( * ), FTF1( 2*NPHP ),
     1                 FTF2( 2*NPHP ), FTF3( 2*NPHP ),
     2                 V1( NPHP, NTHP, 3)
      DOUBLE PRECISION GAUX( NTHP ),
     1                 GAUW( NTHP ), XARR( NR ), ZCF( NR ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
      DOUBLE PRECISION RQST1( LH*(LH+2) ,3 ,NR ),
     1                 RQST2( LH*(LH+2) ,3 ,NR ),
     2                 SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     3                 FDCM( NFDCM, NR, NDRVM ), FAC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR
      DOUBLE PRECISION DLOW
      PARAMETER ( DLOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Quick return for vanishingly small FAC
C
      IF ( ABS( FAC ).LT.DLOW ) RETURN
C
C Loop around grid nodes from ILN to IRN and at each
C evaluate velocity, VI, as a function of (r, theta, phi)
C then calculate k x V and put results (in RQST format into
C the array RQST1.
C
      DO IR = ILN, IRN
C       .
C       . evaluate VI in space ...
C       .
        CALL ASV2VF( NDCS, NTHP, NPHP, LH, IRES, NBN, NFDCM,
     1               IR, NR, NDRVS, NDRVM, INIA, MIT, MIL, MIM,
     2               MIP, MMAX, 'Velocity', VI, SVFDC, GAUX, PA,
     3               DPA, V1, FTF1, FTF2, FTF3, XARR )
C       .
C       . Now evaluate the action of the Coriolis force
C       .
        CALL VFCOR ( NTHP, NPHP, V1, GAUX )
C       .
C       . now put vector back into RQST1 in appropriate space
C       .
        CALL VFRQST ( RQST1, V1, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1                FTF3, ZCF, LH, NTHP, NPHP, MMAX, NR, IR )
C       .
      ENDDO
C
C RQST1 now contains ( k times V ) -
C so now evaluate curl and put into RQST2
C
      CALL RQSTCL( LH, NR, M0, MMAX, NBN, NFDCM, NDRVS, ILN,
     1             IRN, ILN, IRN, RQST1, RQST2, XARR, FDCM )
C
C Now, just need to add FAC*the poloidal and toroidal components
C to the correct elements of VO
C
      CALL RQSTSV( NR, LH, ILN, IRN, INOA, MOT, MOL, MOM,
     1             'Velocity', RQST2, VO, XARR, FAC )
C
      RETURN
      END
C*********************************************************************
