C*********************************************************************
C double precision function Stored Derivative Single harmonic and  ***
C                           -      -          -                    ***
C                          node Magnetic Advection Terms function. ***
C                               -        -         -               ***
C Steve Gibbons Sat Feb  5 10:53:01 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C For a velocity harmonic number IHA (velocity stored in V0A, V1A    C
C and V2A) and mag. field harmonic number IHB (stored in V0B, V1B    C
C and V2B), vectors and INARRA and INARRB may be identical ofcourse  C
C SDSMAT will evaluate the contributions as defined in equations     C
C (B.74) to (B.83) of my thesis.                                     C
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
C     LA        : Spherical harmonic degree of velocity harm. alpha  C
C     IHA       : Number of harmonic in velocity expansion.          C
C     INARRA    : Format of VECA containing velocity radial function C
C                 Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                 INARRA( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                 INARRA( 2 ) = NRR. Must be consistent with NR.     C
C                 INARRA( 3 ) = NHA = total number of radial func.s  C
C                                                                    C
C     LB        : Spherical harmonic degree of mag. fld. harm. beta  C
C     IHB       : Number of harmonic in mag. fld. expansion.         C
C     INARRB    : Format of VECB containing mag. f. radial function  C
C                 Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                 INARRB( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                 INARRB( 2 ) = NRR. Must be consistent with NR.     C
C                 INARRB( 3 ) = NHB = total number of radial func.s  C
C                                                                    C
C     LG        : Spherical harmonic degree of mag. fld. harm. gamma C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Value of radius.                                   C
C     V0A       : 0 derivatives of alpha velocity. Dim ( * )         C
C     V1A       : 1 derivatives of alpha velocity. Dim ( * )         C
C     V2A       : 2 derivatives of alpha velocity. Dim ( * )         C
C     V0B       : 0 derivatives of beta velocity. Dim ( * )          C
C     V1B       : 1 derivatives of beta velocity. Dim ( * )          C
C     V2B       : 2 derivatives of beta velocity. Dim ( * )          C
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
      FUNCTION SDSMAT( IR, RAD, TYPE, LA, IHA, INARRA, LB, IHB,
     1                 INARRB, LG, V0A, V1A, V2A, V0B, V1B, V2B )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, LA, IHA, INARRA( 3 ), LB, IHB, INARRB( 3 ), LG
      DOUBLE PRECISION RAD, SDSMAT, V0A( * ), V1A( * ), V2A( * ),
     1                 V0B( * ), V1B( * ), V2B( * )
      CHARACTER *(4) TYPE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER INDA, INDB, INDFUN
      DOUBLE PRECISION FLA, FLB, SQRLL1, TEMP, LOW, FLG, DLINE1,
     1                 DLINE2, RAD2, RAD3
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( RAD.LT.LOW ) THEN
        PRINT *,' Function SDSMAT.'
        PRINT *,' RAD = ', RAD
        PRINT *,' A division by zero will occur.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RAD2 = RAD*RAD
      RAD3 = RAD2*RAD
C
      INDA = INDFUN( IR, IHA, INARRA )
      INDB = INDFUN( IR, IHB, INARRB )
C
      IF ( TYPE.EQ.'CQST' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Derivatives of poloidal alpha and 
C       . poloidal beta radial functions. 
C       .
C       V0A( INDA )  contains P_{ALPHA}( r )
C       .
C       V0B( INDB )  contains P_{BETA}( r )
C       V1B( INDB )  contains d P_{BETA}( r )/ dr
C       .
        TEMP = V1B( INDB ) + V0B( INDB )/RAD
        SDSMAT = (-1.0d0)*FLA*FLA*FLB*V0A( INDA )*TEMP/(RAD*FLG)
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CSQT' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Derivatives of poloidal alpha and
C       . poloidal beta radial functions.
C       .
C       V0A( INDA )  contains P_{ALPHA}( r )
C       V1A( INDA )  contains d P_{ALPHA}( r )/ dr
C       .
C       V0B( INDB )  contains P_{BETA}( r )
C       .
        TEMP = V1A( INDA ) + V0A( INDA )/RAD
        SDSMAT = (-1.0d0)*FLB*FLB*FLA*V0B( INDB )*TEMP/(RAD*FLG)
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CSSQ' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Derivatives of poloidal alpha and
C       . poloidal beta radial functions.
C       .
C       V0A( INDA )  contains P_{ALPHA}( r )
C       V1A( INDA )  contains d P_{ALPHA}( r )/ dr
C       .
C       V0B( INDB )  contains P_{BETA}( r )
C       V1B( INDB )  contains d P_{BETA}( r )/ dr
C       .
        DLINE1 = V0A( INDA )*V0B( INDB )/RAD3 +
     1           V0A( INDA )*V1B( INDB )/RAD2
        DLINE2 = V1A( INDA )*V0B( INDB )/RAD2 +
     1           V1A( INDA )*V1B( INDB )/RAD
        TEMP   = DLINE1 + DLINE2
        SDSMAT = FLA*FLB*TEMP
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CQSS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Derivatives of poloidal alpha and
C       . poloidal beta radial functions.
C       .
C       V0A( INDA )  contains P_{ALPHA}( r )
C       V1A( INDA )  contains d P_{ALPHA}( r )/ dr
C       .
C       V0B( INDB )  contains P_{BETA}( r )
C       V1B( INDB )  contains d P_{BETA}( r )/ dr
C       V2B( INDB )  contains d^2 P_{BETA}( r )/ dr^2
C       .
        DLINE1 = V1A( INDA )*V0B( INDB )/RAD2 +
     1           V1A( INDA )*V1B( INDB )/RAD  +
     2           V0A( INDA )*V2B( INDB )/RAD 
        DLINE2 = V0A( INDA )*V1B( INDB )/RAD2 -
     1           V0A( INDA )*V0B( INDB )/RAD3
        TEMP   = DLINE1 + DLINE2
        SDSMAT = (-1.0d0)*FLA*FLA*FLB*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CSQS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Derivatives of poloidal alpha and
C       . poloidal beta radial functions.
C       .
C       V0A( INDA )  contains P_{ALPHA}( r )
C       V1A( INDA )  contains d P_{ALPHA}( r )/ dr
C       V2A( INDA )  contains d^2 P_{ALPHA}( r )/ dr^2
C       .
C       V0B( INDB )  contains P_{BETA}( r )
C       V1B( INDB )  contains d P_{BETA}( r )/ dr
C       .
        DLINE1 = V1B( INDB )*V0A( INDA )/RAD2 +
     1           V1B( INDB )*V1A( INDA )/RAD  +
     2           V0B( INDB )*V2A( INDA )/RAD 
        DLINE2 = V0B( INDB )*V1A( INDA )/RAD2 -
     1           V0B( INDB )*V0A( INDA )/RAD3
        TEMP   = DLINE1 + DLINE2
        SDSMAT = (-1.0d0)*FLA*FLB*FLB*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CQTT' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Derivatives of poloidal alpha and
C       . toroidal beta radial functions.
C       .
C       V0A( INDA )  contains P_{ALPHA}( r )
C       .
C       V0B( INDB )  contains T_{BETA}( r )
C       .
        TEMP   = V0A( INDA )*V0B( INDB )/RAD
        SDSMAT = FLA*FLA*FLB*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CQTS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Derivatives of poloidal alpha and
C       . toroidal beta radial functions.
C       .
C       V0A( INDA )  contains P_{ALPHA}( r )
C       V1A( INDA )  contains d P_{ALPHA}( r )/ dr
C       .
C       V0B( INDB )  contains T_{BETA}( r )
C       V1B( INDB )  contains d T_{BETA}( r )/ dr
C       .
        TEMP   = V0A( INDA )*V1B( INDB )/RAD +
     1           V1A( INDA )*V0B( INDB )/RAD
        SDSMAT = FLA*FLA*FLB*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CSTQ' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Derivatives of poloidal alpha and
C       . toroidal beta radial functions.
C       .
C       V0A( INDA )  contains P_{ALPHA}( r )
C       V1A( INDA )  contains d P_{ALPHA}( r )/ dr
C       .
C       V0B( INDB )  contains T_{BETA}( r )
C       .
        TEMP   = V0A( INDA )/RAD + V1A( INDA )
        SDSMAT = (-1.0d0)*FLA*FLB*TEMP*V0B( INDB )/RAD
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CTQT' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Derivatives of toroidal alpha and
C       . poloidal beta radial functions.
C       .
C       V0A( INDA )  contains T_{ALPHA}( r )
C       .
C       V0B( INDB )  contains P_{BETA}( r )
C       .
        TEMP   = V0A( INDA )*V0B( INDB )/RAD
        SDSMAT = FLA*FLB*FLB*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CTQS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
        FLG = SQRLL1( LG )
C       .
C       . Derivatives of toroidal alpha and
C       . poloidal beta radial functions.
C       .
C       V0A( INDA )  contains T_{ALPHA}( r )
C       V1A( INDA )  contains d T_{ALPHA}( r )/ dr
C       .
C       V0B( INDB )  contains P_{BETA}( r )
C       V1B( INDB )  contains d P_{BETA}( r )/ dr
C       .
        TEMP   = V0A( INDA )*V1B( INDB )/RAD +
     1           V1A( INDA )*V0B( INDB )/RAD
        SDSMAT = FLA*FLB*FLB*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CTSQ' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Derivatives of toroidal alpha and
C       . poloidal beta radial functions.
C       .
C       V0A( INDA )  contains T_{ALPHA}( r )
C       .
C       V0B( INDB )  contains P_{BETA}( r )
C       V1B( INDB )  contains d P_{BETA}( r )/ dr
C       .
        TEMP   = V1B( INDB ) + V0B( INDB )/RAD
        SDSMAT = (-1.0d0)*FLA*FLB*TEMP*V0A( INDA )/RAD
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'CTTQ' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Derivatives of toroidal alpha and
C       . toroidal beta radial functions.
C       .
C       V0A( INDA )  contains T_{ALPHA}( r )
C       .
C       V0B( INDB )  contains T_{BETA}( r )
C       .
        TEMP   = V0A( INDA )*V0B( INDB )/RAD
        SDSMAT = FLA*FLB*TEMP
        RETURN
C       .
      ENDIF
C     .
      PRINT *,' Function SDSMAT.'
      PRINT *,' TYPE = ', TYPE
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************

