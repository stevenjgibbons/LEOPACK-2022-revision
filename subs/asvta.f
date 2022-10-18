C*********************************************************************
C subroutine Adapted Solution Vector Term Addition *******************
C            -       -        -      -    -        *******************
C Steve Gibbons Tue Oct 26 13:10:21 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Simply adds a factor of FAC * VI( i ) to VO( j ) where i is every  C
C harmonic, IH, such that MHTI( ih ) = ICMPI and grid node such that C
C ILNR .le. ir .le. IRNR - and j is every harmonic, JH, such that    C
C MHTO( jh ) = ICMPO and grid node, jr, such that                    C
C  ILNR .le. ir .le. IRNR                                            C
C                                                                    C
C Useful maybe if adding time derivatives to the forcing term.       C
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
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ASVTA( NR, NDCS, VI, INARRI, MHTI, MHLI, MHMI, ICMPI,
     1                  MHPI, VO, INARRO, MHTO, MHLO, MHMO, ICMPO,
     2                  FAC, ILNR, IRNR, NBN, NDRVS, NDRVM, NFDCM,
     3                  SVFDC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NDCS, INARRI( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        ICMPI, MHPI( * ), INARRO( * ), MHTO( * ), MHLO( * ),
     2        MHMO( * ), ICMPO, ILNR, IRNR, NBN, NDRVS, NDRVM, NFDCM
      DOUBLE PRECISION VI( * ), VO( * ), FAC,
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER INDO, IHI, IHO, NRI, NHI, NRO, NHO, IR,
     1        INDFUN, IS, IHD
      DOUBLE PRECISION DERV( 1 ), D0F
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
         PRINT *,' Subroutine ASVTA.'
         PRINT *,' NR   = ', NR
         PRINT *,' NRI  = ', NRI
         PRINT *,' NRO  = ', NRO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
      DO IHI = 1, NHI
        IF ( MHTI( IHI ).NE.ICMPI ) GOTO 50
C
C ok - so our harmonic is of the correct type
C so let's look for its corresponding harmonic
C in VO
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

