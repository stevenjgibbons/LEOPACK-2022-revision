C*********************************************************************
C double precision function Single Harmonic and Node Velocity dot    *
C                           -      -            -    -               *
C                                        Gradient of Theta function. *
C                                        -           -               *
C Steve Gibbons Wed Jan 19 11:31:14 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C For a velocity harmonic number IHA (velocity stored in VECA)       C
C and temperature harmonic number IHB (temp. stored in VECB) - VECA  C
C and VECB, and INARRA and INARRB may be identical ofcourse -        C
C SHNVGT will evaluate the necessary radial functions and their      C
C derivatives to evaluate the contributions as defined in equations  C
C (B.51) to (B.53) of my thesis.                                     C
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
C     LB        : Spherical harmonic degree of temperature harm. betaC
C     IHB       : Number of harmonic in temperature expansion.       C
C     ISB       : Number of finite difference scheme used for        C
C                  velocity harmonic.                                C
C     INARRB    : Format of VECB containing theta radial function    C
C                 Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                 INARRB( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                 INARRB( 2 ) = NRR. Must be consistent with NR.     C
C                 INARRB( 3 ) = NHB = total number of radial func.s  C
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
C     TYPE      : *(4) One of 'SPQQ', 'SPSS' or 'SPTS'               C
C                 Indicates the interaction type required.           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION SHNVGT( IR, RAD, NBN, TYPE, NR, NDRVS, NDRVM, NFDCM,
     1                 NDCS, SVFDC, LA, VECA, IHA, ISA, INARRA, LB,
     2                 VECB, IHB, ISB, INARRB )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, NBN, NR, NDRVS, NDRVM, NFDCM, NDCS, LA, IHA, ISA,
     1        INARRA( 3 ), LB, IHB, ISB, INARRB( 3 )
      DOUBLE PRECISION RAD, SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                 VECA( * ), VECB( * ), SHNVGT
      CHARACTER *(4) TYPE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IHDA, IHDB
      DOUBLE PRECISION FLA, FLB, SQRLL1, WORKA( 2 ), WORKB( 2 ),
     1                 TEMP, LOW
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( RAD.LT.LOW ) THEN
        PRINT *,' Function SHNVGT.'
        PRINT *,' RAD = ', RAD
        PRINT *,' A division by zero will occur.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( TYPE.EQ.'SPQQ' ) THEN
        FLA = SQRLL1( LA )
C       .
C       . Take derivatives of poloidal and temperature 
C       . radial functions. 
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
C       WORKB( 2 ) now contains d psi_{BETA}( r )/ dr
C       .
        SHNVGT = FLA*FLA*WORKA( 1 )*WORKB( 2 )/RAD
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'SPSS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Take derivatives of poloidal and temperature
C       . radial functions.
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
C       WORKB( 1 ) now contains psi_{BETA}( r )
C       .
        TEMP   = WORKA( 2 ) + WORKA( 1 )/RAD
        SHNVGT = FLA*FLB*TEMP*WORKB( 1 )/RAD
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'SPTS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Take derivatives of toroidal and temperature
C       . radial functions.
C       .
        IHDA = 0
        CALL ASVDR( VECA, IR, ISA, IHA, NBN, IHDA, NFDCM, NR,
     1            NDRVS, NDRVM, WORKA, INARRA, SVFDC, NDCS )
C       .
C       WORKA( 1 ) now contains tor_{ALPHA}( r )
C       .
        IHDB = 0
        CALL ASVDR( VECB, IR, ISB, IHB, NBN, IHDB, NFDCM, NR,
     1            NDRVS, NDRVM, WORKB, INARRB, SVFDC, NDCS )
C       .
C       WORKB( 1 ) now contains psi_{BETA}( r )
C       .
        SHNVGT = (-1.0d0)*FLA*FLB*WORKA( 1 )*WORKB( 1 )/RAD
        RETURN
C       .
      ENDIF
C     .
      PRINT *,' Function SHNVGT.'
      PRINT *,' TYPE = ', TYPE
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************

