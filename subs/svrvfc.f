C*********************************************************************
C subroutine Solution Vector Radial Vector Function with Curl ********
C            -        -      -      -      -             -    ********
C Steve Gibbons Sat May 13 13:35:12 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given a solution vector, SV, with the usual indexing arrays:       C
C MHT, MHL, MHM, MHP and INARR; SVRVFC will fill an RVF type array   C
C VRVF with the velocity, and CRVF with curl of the velocity.        C
C                                                                    C
C Acts from radial grid nodes 2 to NR - 1.                           C
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
C     NTHPTS    : Number of theta points.                            C
C                                                                    C
C     NPHPTS    : Number of phi points.                              C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                 LH is the maximum l which can 'safely' be          C
C                 computed on this theta/phi grid.                   C
C                                                                    C
C     NBN       : Number of bounding nodes.                          C
C                                                                    C
C     NCFM      : Leading dimension of FDCM and SVFDC.               C
C                   At least (2*NBN+1).                              C
C                                                                    C
C     NR        : Number of radial grid nodes in each function.      C
C                                                                    C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NRR. Must be consistent with NR.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     MMAX      : Maximum spherical harmonic order, m.               C
C                                                                    C
C     M0        : Increment in wavenumbers, m. Must be non-zero.     C
C                  M0 = 1 is safest - valid for all cases.           C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     VRVF      : Output. Dim ( NR, NPHPTS, NTHPTS, 3). Vel. Vector  C
C     CRVF      : Output. Dim ( NR, NPHPTS, NTHPTS, 3). Curl of vel. C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C     FTF1      : Work array. Dim (2*NPHPTS).                        C
C     FTF2      : Work array. Dim (2*NPHPTS).                        C
C     FTF3      : Work array. Dim (2*NPHPTS).                        C
C                                                                    C
C     FDCM      : Dim ( NCFM, NR, 1 ). Finite diff. coeff.s          C
C                Must be formed in advance by a call to FDCMBD       C
C                with NDRVS = NDRVM = 1, NLMN = 2, NRMN = NR - 1,    C
C                NLMC = 2, NRMC = NR - 1.                            C
C                                                                    C
C     RQST1     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C     RQST2     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVRVFC( NDCS, NTHPTS, NPHPTS, LH, NBN, NCFM, NR,
     1               NDRVS, NDRVM, INARR, MHT, MHL, MHM, MHP, MMAX,
     2               M0, SV, SVFDC, GAUX, PA, DPA, VRVF, CRVF, XARR,
     3               FTF1, FTF2, FTF3, FDCM, RQST1, RQST2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDCS, NTHPTS, NPHPTS, LH, NBN, NCFM, NR, NDRVS, NDRVM,
     1        INARR( * ), MHT( * ), MHL( * ), MHM( * ), MHP( * ),
     2        MMAX, M0
      DOUBLE PRECISION SV( * ), SVFDC( NCFM, NR, NDRVM+1, NDCS ),
     1                 GAUX( NTHPTS ), XARR( NR ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
      DOUBLE PRECISION VRVF( NR, NPHPTS, NTHPTS, 3),
     1                 CRVF( NR, NPHPTS, NTHPTS, 3),
     2                 FDCM( NCFM, NR, 1 )
      DOUBLE PRECISION FTF1( 2*NPHPTS ), FTF2( 2*NPHPTS ), 
     1                 FTF3( 2*NPHPTS ), RQST1( LH*(LH+2), 3, NR ),
     2                 RQST2( LH*(LH+2), 3, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRES, ICLS, NRMAX, ILNR, IRNR, IOP
      CHARACTER *(3) CHVMFF
      DOUBLE PRECISION A, B
      PARAMETER ( IRES = 1, CHVMFF = 'VEL', NRMAX = 200 )
      DOUBLE PRECISION ZCOEFA( NRMAX )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IOP = 0
      B   = 0.0d0
      CALL VECOP( ZCOEFA, B, NR, IOP )
C
      ILNR = 2
      IRNR = NR - 1
C
      IF ( NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine SVRVFC.'
        PRINT *,' NR = ', NR,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C First evaluate velocity directly into VRVF.
C
      CALL ASVRVF( NDCS, NTHPTS, NPHPTS, LH, IRES, NBN, NCFM,
     1             ILNR, IRNR, NR, NDRVS, NDRVM, INARR, MHT,
     2             MHL, MHM, MHP, MMAX, CHVMFF, SV, SVFDC, GAUX,
     3             PA, DPA, VRVF, FTF1, FTF2, FTF3, XARR )
C
C Now evaluate velocity in RQST format in the array RQST1
C
      CALL ASRQST( NR, NDCS, LH, SV, ILNR, IRNR, INARR,
     1             NCFM, NDRVS, NDRVM, RQST1, NBN, XARR,
     2             MHT, MHL, MHM, MHP, CHVMFF, SVFDC )
C
C Now take curl of velocity and store in RQST2
C
      ICLS = 1
      A    = 0.0d0
      B    = 1.0d0
C
      CALL RQSTCA( LH, NR, M0, MMAX, NBN, NCFM, 1, ILNR, IRNR,
     1             ILNR, IRNR, RQST1, RQST2, XARR, FDCM,
     2             ICLS, A, B )
C
C Now evaluate curl of velocity in space into CRVF.
C
      CALL RQSRVF( RQST2, CRVF, GAUX, PA, DPA, FTF1, FTF2, FTF3,
     1            LH, NTHPTS, NPHPTS, MMAX, ZCOEFA, NR, ILNR, IRNR )
C
      RETURN
      END
C*********************************************************************
