C*********************************************************************
C subroutine Adapted Solution Vector Inhomog. Temp. Laplacian ********
C            -       -        -      -        -     -         ********
C Steve Gibbons Wed Jan 26 11:02:58 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Adds a multiple of FAC* the laplacian of all harmonics IHI in      C
C solution vector VI which satisfy MHTI( IHI ) = 3 to the            C
C corresponding harmonics IHO in vector VO with MHTO( IHO ) = 3.     C
C                                                                    C
C In addition to Theta, the inhomogeneous temperature function is    C
C added by the use of the routine ITFA.                              C
C                                                                    C
C The inhomog. temp. coef.s are passed in through the array CAFIT    C
C (see ITHCAR) and the indices are passed in by the array MHII.      C
C                                                                    C
C Now MHTI defines what each scalar function in a solution vector    C
C represents.                                                        C
C                                                                    C
C         MHTI( IH ) = 1 --> harmonic is poloidal velocity           C
C         MHTI( IH ) = 2 --> harmonic is toroidal velocity           C
C         MHTI( IH ) = 3 --> harmonic is temperature.                C
C         MHTI( IH ) = 4 --> harmonic is poloidal magnetic field     C
C         MHTI( IH ) = 5 --> harmonic is toroidal magnetic field     C
C                                                                    C
C The lapacian of a scalar with radial function f(r)                 C
C is a scalar with radial function D_l f(r) - l is degree of S.H.    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     INARRI    : Int. parameter array corresponding to VI.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARRI( 1 ) = IFORMF                               C
C                 INARRI( 2 ) = NRI     See INDFUN for details       C
C                 INARRI( 3 ) = NHI                                  C
C                                                                    C
C     MHTI      : Array length ( * ) - atleast length NHI            C
C                  See above for key. (corresponds to input vec.)    C
C     MHLI      : Array length ( * ) - atleast length NHI            C
C                  Sph. harm. degree, l.                             C
C     MHMI      : Array length ( * ) - atleast length NHI            C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHPI      : Array length ( * ) - atleast length NHI            C
C                  Pointer array to finite difference coefficients.  C
C                  MHPI( ih ) = is, which is the 4th index of        C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     Note that INARRO, MHTO, MHLO and MHMO correspond to            C
C    exactly the above variables but for the output and not input    C
C    vectors. NRI must equal NRO and must both equal NR.             C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NDRVS     : Highest derivative stored in SVFDC.                C
C                (Not needed if we are doing tor --> pol but         C
C                 must be atleast 2 if doing pol --> tor )           C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C     MHII     : Dim ( * ) - atleast atleast length NHI.             C
C                For each temperature harmonic, IH, MHII( IH ) gives C
C                the index of the array CAFIT which stores the       C
C                coefficients for that radial function.              C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C                                                                    C
C                If MHII( IH ) = IITH, then CA, CB and CC are        C
C                respectively stored in CAFIT( 1, IITH ),            C
C                CAFIT( 2, IITH ) and CAFIT( 3, IITH ).              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VI        : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NRI*NHI                     C
C     VO        : Solution vector. Dim ( * )  ( output )             C
C                 Length must be atleast NRO*NHO                     C
C     FAC       : Multiplier of curl to be added.                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     CAFIT     : Dimension ( 3, * ). See MHII.                      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ASVITL( NR, NDCS, VI, INARRI, MHTI, MHLI, MHMI,
     1                  MHPI, VO, INARRO, MHTO, MHLO, MHMO, FAC, NBN,
     2                  NDRVS, NDRVM, ILNR, IRNR, SVFDC, XARR, NFDCM,
     3                  MHII, CAFIT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NDCS, INARRI( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        MHPI( * ), INARRO( * ), MHTO( * ), MHLO( * ), MHII( * ),
     2        MHMO( * ), NBN, NDRVS, NDRVM, ILNR, IRNR, NFDCM
      DOUBLE PRECISION VI( * ), VO( * ), FAC, XARR( NR ),
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     2                 CAFIT( 3, * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IFORMI, IFORMO, NRI, NRO, NHI, NHO, IHI, IHO, IR,
     1        IS, INDO, IHD, INDFUN, L, ICMP, IITH
      DOUBLE PRECISION LOW, RAD, DL, D0F, D1F, D2F, DERV( 3 ),
     1                 RI, RO, CA, CB, CC
      PARAMETER ( LOW = 1.0d-9, ICMP = 3 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Don't bother checking ILNR, IRNR, NDRVS or NBN
C as these errors will all be trapped by ASVDR
C
      IFORMI = INARRI( 1 )
      NRI    = INARRI( 2 )
      NHI    = INARRI( 3 )
C
      IFORMO = INARRO( 1 )
      NRO    = INARRO( 2 )
      NHO    = INARRO( 3 )
C
      IF ( IFORMI.EQ.1 .OR. IFORMI.EQ.2 .OR.
     1     IFORMO.EQ.1 .OR. IFORMO.EQ.2 ) THEN
         PRINT *,' Subroutine ASVITL.'
         PRINT *,' IFORMI = ', IFORMI
         PRINT *,' IFORMO = ', IFORMO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( NRI.NE.NR .OR. NRO.NE.NR ) THEN
         PRINT *,' Subroutine ASVITL.'
         PRINT *,' NR   = ', NR
         PRINT *,' NRI  = ', NRI
         PRINT *,' NRO  = ', NRO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C     . early exit?
C
      IF ( ABS(FAC).LT.LOW ) RETURN
      IHD = 2
C     .
      RI = XARR(  1 )
      RO = XARR( NR )
C     .
C     . Loop around in harmonics
C     .
      DO IHI = 1, NHI
        IF ( MHTI( IHI ).NE.ICMP ) GOTO 50
C
C ok - so our harmonic is of the correct type
C retrieve inhomogeneous temperature coefficients
C
        IITH   = MHII( IHI )
        CA     = CAFIT( 1, IITH )
        CB     = CAFIT( 2, IITH )
        CC     = CAFIT( 3, IITH )
C
C so let's look for its corresponding harmonic
C in VO
C
        DO IHO = 1, NHO
          IF (   MHTO( IHO ).EQ.ICMP          .AND.
     1           MHLO( IHO ).EQ.MHLI( IHI )   .AND.
     2           MHMO( IHO ).EQ.MHMI( IHI )         ) THEN
C
C o.k. we've found the corresponding harmonic
C
            DO IR = ILNR, IRNR
C
C Find locations of source and destination vectors ...
C
             INDO = INDFUN( IR, IHO, INARRO )
C            .
             L   = MHLI( IHI )
             IS  = MHPI( IHI )
             RAD = XARR( IR )
             CALL ASVDR( VI, IR, IS, IHI, NBN, IHD, NFDCM, NR,
     1                   NDRVS, NDRVM, DERV, INARRI, SVFDC, NDCS )
C
C Add on inhomogeneous part of temperature
C
             CALL ITFA( RAD, RI, RO, CA, CB, CC, DERV, IHD )
C
             D0F = DERV( 1 )
             D1F = DERV( 2 )
             D2F = DERV( 3 )
             VO( INDO ) = VO( INDO )
     1                       + FAC*DL( L, RAD, D0F, D1F, D2F )
C            .
            ENDDO
C
          ENDIF 
        ENDDO
C
 50   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
