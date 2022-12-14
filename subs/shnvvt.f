C*********************************************************************
C double precision function Single Harmonic and Node curl V x curl V *
C                           -      -            -         -        - *
C                                                    Terms function. *
C                                                    -               *
C Steve Gibbons Wed Jan 26 15:02:35 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C For a velocity harmonic number IHA (velocity stored in VECA)       C
C and velocity harmonic number IHB (stored in VECB) - VECA           C
C and VECB, and INARRA and INARRB may be identical ofcourse -        C
C SHNVVT will evaluate the necessary radial functions and their      C
C derivatives to evaluate the contributions as defined in equations  C
C (B.63) to (B.72) of my thesis.                                     C
C                                                                    C
C This is evaluated at the grid node IR.                             C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Number of radial grid node.                        C
C     NBN       : Number of bounding nodes. See above.               C
C     NR        : Number of radial grid nodes in each function.      C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C     NFDCM     : Leading dimension of SVFDC. At least (2*NBN+1)     C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     LA        : Spherical harmonic degree of velocity harm. alpha  C
C     IHA       : Number of harmonic in velocity expansion.          C
C     ISA       : Number of finite difference scheme used for        C
C                  velocity harmonic.                                C
C     INARRA    : Format of VECA containing velocity radial function C
C                 Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                 INARRA( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                 INARRA( 2 ) = NRR. Must be consistent with NR.     C
C                 INARRA( 3 ) = NHA = total number of radial func.s  C
C                                                                    C
C     LB        : Spherical harmonic degree of velocity harm. beta   C
C     IHB       : Number of harmonic in velocity expansion.          C
C     ISB       : Number of finite difference scheme used for        C
C                  velocity harmonic.                                C
C     INARRB    : Format of VECB containing beta  radial function    C
C                 Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                 INARRB( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                 INARRB( 2 ) = NRR. Must be consistent with NR.     C
C                 INARRB( 3 ) = NHB = total number of radial func.s  C
C                                                                    C
C     LG        : Spherical harmonic degree of velocity harm. gamma  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Value of radius.                                   C
C     VECA      : Solution vector for velocity. Dim ( * )            C
C     VECB      : Solution vector for velocity. Dim ( * )            C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TYPE      : *(4) One of 'CQTT', 'CQTS' or 'CSTQ'               C
C                 Indicates the interaction type required.           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION SHNVVT( IR, RAD, NBN, TYPE, NR, NDRVS, NDRVM, NFDCM,
     1                 NDCS, SVFDC, LA, VECA, IHA, ISA, INARRA, LB,
     2                 VECB, IHB, ISB, INARRB, LG )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, NBN, NR, NDRVS, NDRVM, NFDCM, NDCS, LA, IHA, ISA,
     1        INARRA( 3 ), LB, IHB, ISB, INARRB( 3 ), LG
      DOUBLE PRECISION RAD, SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                 VECA( * ), VECB( * ), SHNVVT
      CHARACTER *(4) TYPE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IHDA, IHDB
      DOUBLE PRECISION FLA, FLB, SQRLL1, WORKA( 4 ), WORKB( 4 ),
     1                 TEMP, LOW, FLG, DLINE1, DLINE2, DLINE3,
     2                 DLINE4, RAD2, RAD3, RAD4, DL
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( RAD.LT.LOW ) THEN
        PRINT *,' Function SHNVVT.'
        PRINT *,' RAD = ', RAD
        PRINT *,' A division by zero will occur.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RAD2 = RAD*RAD
      RAD3 = RAD2*RAD
      RAD4 = RAD3*RAD
C
      IF ( TYPE.EQ.'CQTT' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Take derivatives of poloidal alpha and 
C       . poloidal beta radial functions. 
C       .
        IHDA = 0
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains P_{ALPHA}( r )
C       .
        IHDB = 2
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains P_{BETA}( r )
C       WORKB( 2 ) now contains d P_{BETA}( r )/ dr
C       WORKB( 3 ) now contains d^2 P_{BETA}( r )/ dr^2
C       .
        TEMP = DL( LB, RAD, WORKB( 1 ), WORKB( 2 ), WORKB( 3 ) )
        SHNVVT = FLA*FLA*FLB*WORKA( 1 )*TEMP/(RAD*FLG)
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CQTS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Take derivatives of poloidal alpha and 
C       . poloidal beta radial functions.
C       .
        IHDA = 1
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains P_{ALPHA}( r )
C       WORKA( 2 ) now contains d P_{ALPHA}( r )/ dr
C       .
        IHDB = 3
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains P_{BETA}( r )
C       WORKB( 2 ) now contains d P_{BETA}( r )/ dr
C       WORKB( 3 ) now contains d^2 P_{BETA}( r )/ dr^2
C       WORKB( 4 ) now contains d^3 P_{BETA}( r )/ dr^3
C       .
        DLINE1 = 2.0d0*FLB*FLB*WORKA( 1 )*WORKB( 1 )/RAD4 -
     1           FLB*FLB*WORKA( 2 )*WORKB( 1 )/RAD3
        DLINE2 = ( -2.0d0 - FLB*FLB )*WORKA( 1 )*WORKB( 2 )/RAD3
        DLINE3 = (2.0d0/RAD2)*(   WORKA( 2 )*WORKB( 2 )
     1                          + WORKA( 1 )*WORKB( 3 )   )
        DLINE4 = WORKA( 2 )*WORKB( 3 )/RAD +
     1           WORKA( 1 )*WORKB( 4 )/RAD
        TEMP   = DLINE1 + DLINE2 + DLINE3 + DLINE4
        SHNVVT = FLA*FLA*FLB*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CSTQ' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Take derivatives of poloidal alpha and
C       . poloidal beta radial functions.
C       .
        IHDA = 1
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains P_{ALPHA}( r )
C       WORKA( 2 ) now contains d P_{ALPHA}( r )/ dr
C       .
        IHDB = 2
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains P_{BETA}( r )
C       WORKB( 2 ) now contains d P_{BETA}( r )/ dr
C       WORKB( 3 ) now contains d^2 P_{BETA}( r )/ dr^2
C       .
        DLINE1 = (-1.0d0*FLB*FLB)*( WORKA( 1 )*WORKB( 1 )/RAD4
     1                            + WORKA( 2 )*WORKB( 1 )/RAD3   )
        DLINE2 = 2.0d0*WORKB( 2 )*( WORKA( 1 )/RAD3 + WORKA( 2 )/RAD2)
        DLINE3 = WORKB( 3 )*( WORKA( 1 )/RAD2 + WORKA( 2 )/RAD )
        TEMP   = DLINE1 + DLINE2 + DLINE3
        SHNVVT = (-1.0d0)*FLA*FLB*TEMP
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CQST' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Take derivatives of poloidal alpha and
C       . toroidal beta radial functions.
C       .
        IHDA = 0
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains P_{ALPHA}( r )
C       .
        IHDB = 1
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains T_{BETA}( r )
C       WORKB( 2 ) now contains d T_{BETA}( r )/ dr
C       .
        TEMP   = ( WORKB( 1 )/RAD + WORKB( 2 ) )/RAD
        SHNVVT = FLA*FLA*FLB*WORKA( 1 )*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CSQT' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Take derivatives of poloidal alpha and
C       . toroidal beta radial functions.
C       .
        IHDA = 1
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains P_{ALPHA}( r )
C       WORKA( 2 ) now contains d P_{ALPHA}( r )/ dr
C       .
        IHDB = 0
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains T_{BETA}( r )
C       .
        TEMP   = ( WORKA( 1 )/RAD + WORKA( 2 ) )/RAD
        SHNVVT = FLB*FLB*FLA*WORKB( 1 )*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CSSQ' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Take derivatives of poloidal alpha and
C       . toroidal beta radial functions.
C       .
        IHDA = 1
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains P_{ALPHA}( r )
C       WORKA( 2 ) now contains d P_{ALPHA}( r )/ dr
C       .
        IHDB = 1
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains T_{BETA}( r )
C       WORKB( 2 ) now contains d T_{BETA}( r )/ dr
C       .
        DLINE1 = WORKA( 1 )*WORKB( 1 )/RAD3 +
     1           WORKA( 1 )*WORKB( 2 )/RAD2
        DLINE2 = WORKA( 2 )*WORKB( 1 )/RAD2 +
     1           WORKA( 2 )*WORKB( 2 )/RAD
        TEMP   = DLINE1 + DLINE2
        SHNVVT = (-1.0d0)*FLA*FLB*TEMP
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CQSS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Take derivatives of poloidal alpha and
C       . toroidal beta radial functions.
C       .
        IHDA = 1
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains P_{ALPHA}( r )
C       WORKA( 2 ) now contains d P_{ALPHA}( r )/ dr
C       .
        IHDB = 2
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains T_{BETA}( r )
C       WORKB( 2 ) now contains d T_{BETA}( r )/ dr
C       WORKB( 3 ) now contains d^2 T_{BETA}( r )/ dr^2
C       .
        DLINE1 = WORKA( 2 )*WORKB( 1 )/RAD2 +
     1           WORKA( 2 )*WORKB( 2 )/RAD  +
     2           WORKA( 1 )*WORKB( 3 )/RAD 
        DLINE2 = WORKA( 1 )*WORKB( 2 )/RAD2 -
     1           WORKA( 1 )*WORKB( 1 )/RAD3
        TEMP   = DLINE1 + DLINE2
        SHNVVT = FLA*FLA*FLB*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CSQS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Take derivatives of poloidal alpha and
C       . toroidal beta radial functions.
C       .
        IHDA = 2
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains P_{ALPHA}( r )
C       WORKA( 2 ) now contains d P_{ALPHA}( r )/ dr
C       WORKA( 3 ) now contains d^2 P_{ALPHA}( r )/ dr^2
C       .
        IHDB = 1
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains T_{BETA}( r )
C       WORKB( 2 ) now contains d T_{BETA}( r )/ dr
C       .
        DLINE1 = WORKB( 2 )*WORKA( 1 )/RAD2 +
     1           WORKB( 2 )*WORKA( 2 )/RAD  +
     2           WORKB( 1 )*WORKA( 3 )/RAD 
        DLINE2 = WORKB( 1 )*WORKA( 2 )/RAD2 -
     1           WORKB( 1 )*WORKA( 1 )/RAD3
        TEMP   = DLINE1 + DLINE2
        SHNVVT = FLB*FLB*FLA*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CTTQ' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Take derivatives of toroidal alpha and
C       . poloidal beta radial functions.
C       .
        IHDA = 0
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains T_{ALPHA}( r )
C       .
        IHDB = 2
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains P_{BETA}( r )
C       WORKB( 2 ) now contains d P_{BETA}( r )/ dr
C       WORKB( 3 ) now contains d^2 P_{BETA}( r )/ dr^2
C       .
        TEMP   = DL( LB, RAD, WORKB( 1 ), WORKB( 2 ), WORKB( 3 ) )
        SHNVVT = FLA*FLB*TEMP*WORKA( 1 )/RAD
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CTQT' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Take derivatives of toroidal alpha and
C       . toroidal beta radial functions.
C       .
        IHDA = 0
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains T_{ALPHA}( r )
C       .
        IHDB = 0
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains T_{BETA}( r )
C       .
        TEMP   = WORKA( 1 )*WORKB( 1 )/RAD
        SHNVVT = (-1.0d0)*FLA*FLB*FLB*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CTQS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Take derivatives of toroidal alpha and
C       . toroidal beta radial functions.
C       .
        IHDA = 1
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains T_{ALPHA}( r )
C       WORKA( 2 ) now contains d T_{ALPHA}( r )/dr
C       .
        IHDB = 1
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains T_{BETA}( r )
C       WORKB( 2 ) now contains d T_{BETA}( r )/dr
C       .
        TEMP   = WORKA( 1 )*WORKB( 2 )/RAD +
     1           WORKA( 2 )*WORKB( 1 )/RAD
        SHNVVT = (-1.0d0)*FLA*FLB*FLB*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CTSQ' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Take derivatives of toroidal alpha and
C       . toroidal beta radial functions.
C       .
        IHDA = 0
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains T_{ALPHA}( r )
C       .
        IHDB = 1
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains T_{BETA}( r )
C       WORKB( 2 ) now contains d T_{BETA}( r )/dr
C       .
        TEMP   = WORKB( 1 )/RAD + WORKB( 2 )
        SHNVVT = FLA*FLB*TEMP*WORKA( 1 )/RAD
        RETURN
C       .
      ENDIF
C     .
      PRINT *,' Function SHNVVT.'
      PRINT *,' TYPE = ', TYPE
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************

