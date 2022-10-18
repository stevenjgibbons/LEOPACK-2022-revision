C*********************************************************************
C subroutine Adapted Solution Vector Cross Product *******************
C            -       -        -      -     -       *******************
C Steve Gibbons Mon Nov  1 19:49:35 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Given a vector V1 (defined by MT1, ML1, MM1 and MP1 ) and a vector C
C V2 (defined by MT2, ML2, MM2 and MP2 ), ASVCP will calculate       C
C which is V1 \times V2 and store the result in the array RQST.      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     NBN       : Number of bounding nodes.                          C
C     IRES      : The resolution flag.                               C
C                 If IRES = 1, harmonics with l greater than LH      C
C                 or m greater than MMAX are simply ignored.         C
C                 If IRES = 2, the program aborts if LH or MMAX      C
C                 is exceeded.                                       C
C                                                                    C
C     INAR1     : Integer array dimension ( * ). Details of V1.      C
C                  Elements may be arbitrary except for              C
C                  INAR1( 1 ) = IFORM1 - flag for vector format.     C
C                                        See INDFUN                  C
C                  INAR1( 2 ) = NRR. Must be consistent with NR.     C
C                  INAR1( 3 ) = NH1 = total number of radial func.s  C
C                                                                    C
C     MT1       : Array length ( * ) - atleast length NH1            C
C                                                                    C
C         MT1( IH ) = 1 --> harmonic is poloidal velocity            C
C         MT1( IH ) = 2 --> harmonic is toroidal velocity            C
C         MT1( IH ) = 3 --> harmonic is temperature.                 C
C         MT1( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MT1( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     ML1       : Array length ( * ) - atleast length NH1            C
C                  Sph. harm. degree, l.                             C
C     MM1       : Array length ( * ) - atleast length NH1            C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MP1       : Array length ( * ) - atleast length NH1            C
C                  Pointer array to finite difference coefficients.  C
C                  MP1( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C Note that INAR2, MT2, ML2, MM2 and MP2 are entirely analogous      C
C to INAR1, MT1, ML1, MM1 and MP1 but for vector 2.                  C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum mode number to be used.                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C     NTHPTS    : Number of points in theta.                         C
C     NPHPTS    : Number of points in phi.                           C
C                                                                    C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dimension of FDCM. At least (2*NBN+1)      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V1        : First vector. Dim ( * ) aleast NH1*NR              C
C     V2        : Second vector. Dim ( * ) aleast NH2*NR             C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Gauss weights computed bu GAUWTS. Dim ( NTHPTS )   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C       PA and DPA must be formed in advance by a call to SCHNLA.    C
C                                                                    C
C  RQST and ZCF store the qst decomposition of the output vector.    C
C                                                                    C
C   ZCF. ( j ) contains the coeff. to Q_0^0( r_j )                   C
C                                                                    C
C     RQST      : Dim. (  LH*(LH+2) ,3, NR )                         C
C     ZCF       : Dim. (  NR )                                       C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHVM1     : Velocity/magnetic field select for V1.             C
C     CHVM2     : Velocity/magnetic field select for V2.             C
C                 Set to either 'VEL' or 'MAG'                       C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VF1       : Storage for vector function. Dim.(NPHPTS,NTHPTS,3) C
C     VF2       : Storage for vector function. Dim.(NPHPTS,NTHPTS,3) C
C     VF3       : Storage for vector function. Dim.(NPHPTS,NTHPTS,3) C
C     FTF1      : Array for fourier transforming. Dim. ( 2*NPHPTS )  C
C     FTF2      : Array for fourier transforming. Dim. ( 2*NPHPTS )  C
C     FTF3      : Array for fourier transforming. Dim. ( 2*NPHPTS )  C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ASVCP( NDCS, NBN, IRES, INAR1, MT1, ML1, MM1, MP1,
     1           INAR2, MT2, ML2, MM2, MP2, NR, LH, MMAX, ILNR,
     2           IRNR, NTHPTS, NPHPTS, NDRVS, NDRVM, NFDCM, V1, V2,
     3           GAUX, GAUW, PA, DPA, RQST, ZCF, SVFDC, XARR, CHVM1,
     4           CHVM2, VF1, VF2, VF3, FTF1, FTF2, FTF3 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDCS, NBN, IRES, INAR1( * ), MT1( * ), ML1( * ),
     1        MM1( * ), MP1( * ), INAR2( * ), MT2( * ), ML2( * ),
     2        MM2( * ), MP2( * ), NR, LH, MMAX, ILNR, IRNR,
     3        NTHPTS, NPHPTS, NFDCM, NDRVM, NDRVS
      DOUBLE PRECISION RQST( LH*(LH+2) ,3 ,NR ), ZCF( NR ),
     1                 V1( * ), V2( * ), XARR( NR ),
     2                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
      DOUBLE PRECISION FTF1( 2*NPHPTS ),
     1                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
      DOUBLE PRECISION VF1( NPHPTS, NTHPTS, 3), GAUX( NTHPTS ),
     1                 VF2( NPHPTS, NTHPTS, 3), GAUW( NTHPTS ),
     2                 VF3( NPHPTS, NTHPTS, 3)
      CHARACTER *(3) CHVM1, CHVM2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Loop around points IR = ILNR, IRNR
C
      DO IR = ILNR, IRNR
C       .
C       . Put the vector function of V1 into VF1 for node IR
C       .
        CALL ASV2VF( NDCS, NTHPTS, NPHPTS, LH, IRES, NBN, NFDCM,
     1               IR, NR, NDRVS, NDRVM, INAR1, MT1, ML1, MM1,
     2               MP1, MMAX, CHVM1, V1, SVFDC, GAUX, PA, DPA,
     3               VF1, FTF1, FTF2, FTF3, XARR )
C       .
C       . Put the vector function of V2 into VF2 for node IR
C       .
        CALL ASV2VF( NDCS, NTHPTS, NPHPTS, LH, IRES, NBN, NFDCM,
     1               IR, NR, NDRVS, NDRVM, INAR2, MT2, ML2, MM2,
     2               MP2, MMAX, CHVM2, V2, SVFDC, GAUX, PA, DPA,
     3               VF2, FTF1, FTF2, FTF3, XARR )
C       .
C       . Take cross product of vectors 1 and 2
C       .
        CALL VFCP( VF1, VF2, VF3, NPHPTS, NTHPTS )
C       .
C       . Transform back into RQST
C       .
        CALL VFRQST( RQST, VF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1               FTF3, ZCF, LH, NTHPTS, NPHPTS, MMAX, NR, IR )
C       .
C
      ENDDO
C
      RETURN
      END
C*********************************************************************
