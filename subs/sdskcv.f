C*********************************************************************
C double precision function Stored Derivatives Single harmonic and ***
C                           -      -           -                   ***
C                                             node K Cross V func. ***
C                                                  - -     -       ***
C Steve Gibbons Wed Jan 26 11:40:35 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C For a velocity harmonic number IHA (velocity stored in V0A, V1A    C
C and V2A) SDSKCV will evaluate the contributions as defined in      C
C equations (B.55) to (B.59) of my thesis.                           C
C                                                                    C
C This is evaluated at the grid node IR.                             C
C                                                                    C
C Only the degree, l, of the output (gamma) harmonic is required     C
C and is supplied in the integer variable LG.                        C
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
C     LG        : Spherical harmonic degree of velocity harm. gamma  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Value of radius.                                   C
C     V0A       : 0 derivatives of alpha velocity. Dim ( * )         C
C     V1A       : 1 derivatives of alpha velocity. Dim ( * )         C
C     V2A       : 2 derivatives of alpha velocity. Dim ( * )         C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TYPE      : *(4) One of 'KCQT', 'KCST', 'KCSQ', etc.           C
C                 Indicates the interaction type required.           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION SDSKCV( IR, RAD, TYPE, LA, IHA, INARRA, LG,
     1                 V0A, V1A, V2A )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, LA, IHA, INARRA( 3 ), LG
      DOUBLE PRECISION RAD, SDSKCV, V0A( * ), V1A( * ), V2A( * )
      CHARACTER *(4) TYPE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER INDA, INDFUN
      DOUBLE PRECISION FLA, FLG, SQRLL1, TEMP, LOW
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( RAD.LT.LOW ) THEN
        PRINT *,' Function SDSKCV.'
        PRINT *,' RAD = ', RAD
        PRINT *,' A division by zero will occur.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INDA = INDFUN( IR, IHA, INARRA )
C
      IF ( TYPE.EQ.'KCQT' ) THEN
        FLA = SQRLL1( LA )
        FLG = SQRLL1( LG )
C       .
C       . Derivative of poloidal radial function. 
C       .
C       V0A( INDA ) now contains P_{ALPHA}( r )
C       .
        SDSKCV = (-1.0d0)*FLA*FLA*V0A( INDA )/(RAD*FLG)
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'KCST' ) THEN
        FLA = SQRLL1( LA )
        FLG = SQRLL1( LG )
C       .
C       . Derivative of poloidal radial function.
C       .
C       V0A( INDA ) now contains P_{ALPHA}( r )
C       V1A( INDA ) now contains d P_{ALPHA}( r )/ dr
C       .
        TEMP   = V1A( INDA ) + V0A( INDA )/RAD
        SDSKCV = (-1.0d0)*FLA*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'KCSQ' ) THEN
        FLA = SQRLL1( LA )
C       .
C       . Derivative of poloidal radial function.
C       .
C       V0A( INDA ) now contains P_{ALPHA}( r )
C       V1A( INDA ) now contains d P_{ALPHA}( r )/ dr
C       .
        TEMP   = V1A( INDA ) + V0A( INDA )/RAD
        SDSKCV = FLA*TEMP/RAD
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'KCSS' ) THEN
        FLA = SQRLL1( LA )
        FLG = SQRLL1( LG )
C       .
C       . Derivative of poloidal radial function.
C       .
C       V1A( INDA ) now contains d P_{ALPHA}( r )/ dr
C       V2A( INDA ) now contains d^2 P_{ALPHA}( r )/ dr^2
C       .
        TEMP   = V2A( INDA ) + 2.0d0*V1A( INDA )/RAD
        SDSKCV = (-1.0d0)*FLA*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'KCQS' ) THEN
        FLA = SQRLL1( LA )
        FLG = SQRLL1( LG )
C       .
C       . Derivative of poloidal radial function.
C       .
C       V1A( INDA ) now contains d P_{ALPHA}( r )/ dr
C       .
        TEMP   = V1A( INDA )/RAD
        SDSKCV = (-1.0d0)*FLA*FLA*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'KCTT' ) THEN
        FLA = SQRLL1( LA )
        FLG = SQRLL1( LG )
C       .
C       . Derivative of toroidal radial function.
C       .
C       V0A( INDA ) now contains tor_{ALPHA}( r )
C       .
        SDSKCV = FLA*V0A( INDA )/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'KCTS' ) THEN
        FLA = SQRLL1( LA )
        FLG = SQRLL1( LG )
C       .
C       . Derivative of toroidal radial function.
C       .
C       V0A( INDA ) now contains tor_{ALPHA}( r )
C       V1A( INDA ) now contains d tor_{ALPHA}( r )/ dr
C       .
        TEMP   = V0A( INDA )/RAD + V1A( INDA )
        SDSKCV = FLA*TEMP/FLG
        RETURN
C       .
      ENDIF
C     .
      IF ( TYPE.EQ.'KCTQ' ) THEN
        FLA = SQRLL1( LA )
C       .
C       . Derivative of toroidal radial function.
C       .
C       V0A( INDA ) now contains tor_{ALPHA}( r )
C       .
        SDSKCV = (-1.0d0)*FLA*V0A( INDA )/RAD
        RETURN
C       .
      ENDIF
C     .
      PRINT *,' Function SDSKCV.'
      PRINT *,' TYPE = ', TYPE
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************

