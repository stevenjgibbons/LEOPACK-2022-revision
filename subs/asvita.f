C*********************************************************************
C subroutine Adapted Solution Vector Inhomogeneous Temperature  Add **
C            -       -        -      -             -            -   **
C Steve Gibbons Wed Jan 26 10:16:17 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Simply adds a factor of FAC * VI( i ) to VO( j ) where i is every  C
C harmonic, IH, such that MHTI( ih ) = 3 and grid node such that     C
C ILNR .le. ir .le. IRNR - and j is every harmonic, JH, such that    C
C MHTO( jh ) = 2 and grid node, jr, such that                        C
C  ILNR .le. ir .le. IRNR                                            C
C                                                                    C
C In addition to Theta, the inhomogeneous temperature function is    C
C added by the use of the routine ITFA.                              C
C                                                                    C
C The inhomog. temp. coef.s are passed in through the array CAFIT    C
C (see ITHCAR) and the indices are passed in by the array MHII.      C
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
C                  See below for key. (corresponds to input vec.)    C
C                                                                    C
C     MHTI( i ) = 1 for a poloidal velocity harmonic, i.             C
C     MHTI( i ) = 2 for a toroidal velocity harmonic, i.             C
C     MHTI( i ) = 3 for a temperature harmonic, i.                   C
C     MHTI( i ) = 4 for a poloidal magnetic field harmonic, i.       C
C     MHTI( i ) = 5 for a toroidal magnetic field harmonic, i.       C
C                                                                    C
C     MHLI      : Array length ( * ) - atleast length NHI            C
C                  Sph. harm. degree, l.                             C
C     MHMI      : Array length ( * ) - atleast length NHI            C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     Note that INARRO, MHTO, MHLO and MHMO correspond to            C
C    exactly the above variables but for the output and not input    C
C    vectors. NRI must equal NRO and must both equal NR.             C
C                                                                    C
C     MHPI      : Array length ( * ) - atleast length NHI            C
C                  Pointer array to finite difference coefficients.  C
C                  MHPI( ih ) = is, which is the 4th index of        C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C     ICMPI     : The type of harmonic we wish to act upon.          C
C     ICMPO     : The type of harmonic we wish add term to.          C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NDRVS     : Highest derivative stored in SVFDC.                C
C                 Must be atleast zero.                              C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
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
C     FAC       : Multiplication factor.                             C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     CAFIT     : Dimension ( 3, * ). See MHII.                      C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ASVITA( NR, NDCS, VI, INARRI, MHTI, MHLI, MHMI,
     1                  MHPI, VO, INARRO, MHTO, MHLO, MHMO,
     2                  FAC, ILNR, IRNR, NBN, NDRVS, NDRVM, NFDCM,
     3                  SVFDC, MHII, CAFIT, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NDCS, INARRI( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        MHPI( * ), INARRO( * ), MHTO( * ), MHLO( * ), MHII( * ),
     2        MHMO( * ), ILNR, IRNR, NBN, NDRVS, NDRVM, NFDCM
      DOUBLE PRECISION VI( * ), VO( * ), FAC, XARR( NR ),
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     2                 CAFIT( 3, * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER INDO, IHI, IHO, NRI, NHI, NRO, NHO, IR,
     1        INDFUN, IS, IHD, ICMPI, ICMPO, IITH
      DOUBLE PRECISION DERV( 1 ), D0F, RI, RO, CA, CB, CC, RAD
      PARAMETER ( ICMPI = 3, ICMPO = 2 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IHD = 0
C
C Check grid nodes ...
C
      NRI    = INARRI( 2 )
      NHI    = INARRI( 3 )
C
      NRO    = INARRO( 2 )
      NHO    = INARRO( 3 )
C     .
      IF ( NRI.NE.NR .OR. NRO.NE.NR ) THEN
         PRINT *,' Subroutine ASVITA.'
         PRINT *,' NR   = ', NR
         PRINT *,' NRI  = ', NRI
         PRINT *,' NRO  = ', NRO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
      RI = XARR(  1 )
      RO = XARR( NR )
C     .
      DO IHI = 1, NHI
        IF ( MHTI( IHI ).NE.ICMPI ) GOTO 50
C
C ok - so our harmonic is of the correct type
C retrieve inhomogeneous temperature coefficients
C
        IITH   = MHII( IHI )
        CA     = CAFIT( 1, IITH )
        CB     = CAFIT( 2, IITH )
        CC     = CAFIT( 3, IITH )
C
C so let's look for its corresponding harmonic in VO
C
        DO IHO = 1, NHO
          IF (   MHTO( IHO ).EQ.ICMPO          .AND.
     1           MHLO( IHO ).EQ.MHLI( IHI )    .AND.
     2           MHMO( IHO ).EQ.MHMI( IHI )         ) THEN
C
C o.k. we've found the corresponding harmonic
C
            IS = MHPI( IHI )
            DO IR = ILNR, IRNR
C
C Find locations of source and destination vectors ...
C
             CALL ASVDR( VI, IR, IS, IHI, NBN, IHD, NFDCM, NR,
     1                   NDRVS, NDRVM, DERV, INARRI, SVFDC, NDCS )
C
C Add on inhomogeneous part of temperature
C
             RAD = XARR( IR )
             CALL ITFA( RAD, RI, RO, CA, CB, CC, DERV, IHD )
C
             INDO = INDFUN( IR, IHO, INARRO )
             D0F = DERV( 1 )
             VO( INDO ) = VO( INDO ) + FAC*D0F
C
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

