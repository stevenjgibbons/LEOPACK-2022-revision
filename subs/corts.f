C*********************************************************************
C subroutine CORiolis Term Subroutine ********************************
C            ---      -    -          ********************************
C Steve Gibbons  20.11.98                                            C
C____________________________________________________________________C
C All the input variables are as described in the above routine and  C
C so will not be described in detail here. The inputs which differ   C
C are INTPAR and DPPARS which depend on the term of curl (K cross V) C
C which is being applied.                                            C
C                                                                    C
C   In ALL cases, INTPAR( 1 ) = IFLAG.                               C
C  IFLAG determines which term is being applied.                     C
C                                                                    C
C  INTPAR( 2 ) is the degree ( L_{\alpha} ) of the input harmonic    C
C  INTPAR( 3 ) is the degree ( L_{\gamma} ) of the output harmonic   C
C                                                                    C
C  Cases                                                             C
C  -----                                                             C
C   IFLAG = 1. This is a poloidal interaction                        C
C              resulting from T_ag^q ( or K_{qt}^{ag} in thesis )    C
C              In my thesis this is equation B.55 - first term.      C
C              DPPARS( 1 ) contains the value K_{qt}^{ag}.           C
C                                                                    C
C   IFLAG = 2. This is a toroidal from poloidal term resulting       C
C              from the S_{ag}^q term ( K_{qs}^{ag} in thesis )      C
C              In my thesis this is equation B.56 - last term.       C
C              DPPARS( 1 ) contains the value K_{qs}^{ag}.           C
C                                                                    C
C   IFLAG = 3. This is a toroidal from poloidal term resulting       C
C              from the Q_{ag}^s term ( K_{sq}^{ag} in thesis )      C
C              In my thesis this is equation B.56 - first term.      C
C              DPPARS( 1 ) contains the value K_{sq}^{ag}.           C
C                                                                    C
C   IFLAG = 4. This is a toroidal from poloidal term resulting       C
C              from the S_{ag}^s term ( K_{ss}^{ag} in thesis )      C
C              In my thesis this is equation B.56 - second term.     C
C              DPPARS( 1 ) contains the value K_{ss}^{ag}.           C
C                                                                    C
C   IFLAG = 5. This is a toroidal from poloidal term resulting       C
C              from the T_{ag}^s term ( K_{st}^{ag} in thesis )      C
C              In my thesis this is equation B.55 - second term.     C
C              DPPARS( 1 ) contains the value K_{st}^{ag}.           C
C                                                                    C
C   IFLAG = 6. This is a toroidal from toroidal term resulting       C
C              from the Q_{ag}^s term ( K_{tq}^{ag} in thesis )      C
C              In my thesis this is equation B.59 - second term.     C
C              DPPARS( 1 ) contains the value K_{tq}^{ag}.           C
C                                                                    C
C   IFLAG = 7. This is a toroidal from toroidal term resulting       C
C              from the S_{ag}^t term ( K_{ts}^{ag} in thesis )      C
C              In my thesis this is equation B.59 -  first term.     C
C              DPPARS( 1 ) contains the value K_{ts}^{ag}.           C
C                                                                    C
C   IFLAG = 8. This is a poloidal from toroidal term resulting       C
C              from the T_{ag}^t term ( K_{tt}^{ag} in thesis )      C
C              In my thesis this is equation B.58 -  first term.     C
C              DPPARS( 1 ) contains the value K_{tt}^{ag}.           C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     INTPAR    : Array of dimension (*)                             C
C                 (see above)                                        C
C     IRN       : Number of the radial grid node.                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Value of the radius.                               C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     VEC       : Array (*) - not referred to.                       C
C     DPPARS    : Array (*) -                                        C
C  Character                                                         C
C  ---------                                                         C
C     ORD       : Order of deriv. accuracy                           C
C____________________________________________________________________C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C     D0FAC     : Multiplies zeroth derivative in matrix             C
C     D1FAC     : Multiplies first derivative in matrix              C
C     D2FAC     : Multiplies second derivative in matrix             C
C     D3FAC     : Multiplies third derivative in matrix              C
C     D4FAC     : Multiplies fourth derivative in matrix             C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Note IRN, RI, RO, VEC and ORD or not referred to by CORTS          C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CORTS ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD,
     1                   RI, RO, INTPAR, DPPARS, VEC, ORD, IRN )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IRN, INTPAR( * )
      DOUBLE PRECISION D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD,
     1                   RI, RO, DPPARS( * ), VEC( * )
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFLAG, LA, LG
      DOUBLE PRECISION FLA, FLG, FAC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFLAG = INTPAR( 1 )
      LA    = INTPAR( 2 )
      LG    = INTPAR( 3 )
C
      FLA = DSQRT( 1.0d0*( LA*LA + LA )  )
      FLG = DSQRT( 1.0d0*( LG*LG + LG )  )
      FAC = DPPARS( 1 )
C
      IF ( IFLAG.EQ.1 ) THEN
         D0FAC = (-1.0d0)*FLA*FLA*FAC/(RAD*FLG)
         D1FAC =  0.0d0
         D2FAC =  0.0d0
         D3FAC =  0.0d0
         D4FAC =  0.0d0
      ENDIF
C
      IF ( IFLAG.EQ.2 ) THEN
         D0FAC =  0.0d0
         D1FAC = (-1.0d0)*FLA*FLA*FAC/(RAD*FLG)
         D2FAC =  0.0d0
         D3FAC =  0.0d0
         D4FAC =  0.0d0
      ENDIF
C
      IF ( IFLAG.EQ.3 ) THEN
         D0FAC =  FLA*FAC/(RAD*RAD)
         D1FAC =  FLA*FAC/RAD
         D2FAC =  0.0d0
         D3FAC =  0.0d0
         D4FAC =  0.0d0
      ENDIF
C
      IF ( IFLAG.EQ.4 ) THEN
         D0FAC =  0.0d0
         D1FAC =  (-2.0d0)*FLA*FAC/(RAD*FLG)
         D2FAC =  (-1.0d0)*FLA*FAC/FLG
         D3FAC =  0.0d0
         D4FAC =  0.0d0
      ENDIF
C
      IF ( IFLAG.EQ.5 ) THEN
         D0FAC =  (-1.0d0)*FLA*FAC/(FLG*RAD)
         D1FAC =  (-1.0d0)*FLA*FAC/FLG
         D2FAC =  0.0d0
         D3FAC =  0.0d0
         D4FAC =  0.0d0
      ENDIF
C
      IF ( IFLAG.EQ.6 ) THEN
         D0FAC =  (-1.0d0)*FLA*FAC/RAD
         D1FAC =  0.0d0
         D2FAC =  0.0d0
         D3FAC =  0.0d0
         D4FAC =  0.0d0
      ENDIF
C
      IF ( IFLAG.EQ.7 ) THEN
         D0FAC =  FLA*FAC/(FLG*RAD)
         D1FAC =  FLA*FAC/FLG
         D2FAC =  0.0d0
         D3FAC =  0.0d0
         D4FAC =  0.0d0
      ENDIF
C
      IF ( IFLAG.EQ.8 ) THEN
         D0FAC =  FLA*FAC/FLG
         D1FAC =  0.0d0
         D2FAC =  0.0d0
         D3FAC =  0.0d0
         D4FAC =  0.0d0
      ENDIF
C
C
      RETURN
      END
C*********************************************************************
