C*********************************************************************
C double precision function Stored Derivatives Single harmonic and ***
C                           -      -           -                   ***
C                    node Velocity dot Gradient of Theta function. ***
C                         -            -           -               ***
C Steve Gibbons Sat Feb  5 11:20:48 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C For a velocity harmonic number IHA (velocity stored in V0A, V1A)   C
C and temperature harmonic number IHB (temp. stored in V0B and V1B)  C
C -  vectors and INARRA and INARRB may be identical ofcourse -       C
C SDSVGT will evaluate the  contributions as defined in equations    C
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
C     LB        : Spherical harmonic degree of temperature harm. betaC
C     IHB       : Number of harmonic in temperature expansion.       C
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
C     V0A       : 0 derivatives of alpha velocity. Dim ( * )         C
C     V1A       : 1 derivatives of alpha velocity. Dim ( * )         C
C     V0B       : 0 derivatives of beta velocity. Dim ( * )          C
C     V1B       : 1 derivatives of beta velocity. Dim ( * )          C
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
      FUNCTION SDSVGT( IR, RAD, TYPE, LA, IHA, INARRA, LB, IHB,
     1                 INARRB, V0A, V1A, V0B, V1B )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, LA, IHA, INARRA( 3 ), LB, IHB, INARRB( 3 )
      DOUBLE PRECISION RAD, SDSVGT, V0A( * ), V1A( * ),
     1                 V0B( * ), V1B( * )
      CHARACTER *(4) TYPE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER INDA, INDB, INDFUN
      DOUBLE PRECISION FLA, FLB, SQRLL1, TEMP, LOW
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( RAD.LT.LOW ) THEN
        PRINT *,' Function SDSVGT.'
        PRINT *,' RAD = ', RAD
        PRINT *,' A division by zero will occur.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INDA = INDFUN( IR, IHA, INARRA )
      INDB = INDFUN( IR, IHB, INARRB )
C
      IF ( TYPE.EQ.'SPQQ' ) THEN
        FLA = SQRLL1( LA )
C       .
C       . Derivatives of poloidal and temperature 
C       . radial functions. 
C       .
C       V0A( INDA ) now contains P_{ALPHA}( r )
C       .
C       V1B( INDB ) now contains d psi_{BETA}( r )/ dr
C       .
        SDSVGT = FLA*FLA*V0A( INDA )*V1B( INDB )/RAD
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'SPSS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Derivatives of poloidal and temperature
C       . radial functions.
C       .
C       V0A( INDA ) now contains P_{ALPHA}( r )
C       V1A( INDA ) now contains d P_{ALPHA}( r )/ dr
C       .
C       V0B( INDB ) now contains psi_{BETA}( r )
C       .
        TEMP   = V1A( INDA ) + V0A( INDA )/RAD
        SDSVGT = FLA*FLB*TEMP*V0B( INDB )/RAD
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'SPTS' ) THEN
        FLA = SQRLL1( LA )
        FLB = SQRLL1( LB )
C       .
C       . Derivatives of toroidal and temperature
C       . radial functions.
C       .
C       V0A( INDA ) now contains tor_{ALPHA}( r )
C       .
C       V0B( INDB ) now contains psi_{BETA}( r )
C       .
        SDSVGT = (-1.0d0)*FLA*FLB*V0A( INDA )*V0B( INDB )/RAD
        RETURN
C       .
      ENDIF
C     .
      PRINT *,' Function SDSVGT.'
      PRINT *,' TYPE = ', TYPE
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************

