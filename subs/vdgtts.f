C*********************************************************************
C subroutine Velocity Dot Gradiert of Theta Term Subroutine **********
C            -        -   -           -     -    -          **********
C Steve Gibbons  30.06.99                                            C
C____________________________________________________________________C
C                                                                    C
C  DPPARS( 1 ) = Scalar spherical harmonic interaction ( c.f.        C
C                 equations B.39, B.40 and B.41 in my thesis )       C
C  INTPAR( 1 ) = IMODE                                               C
C                                                                    C
C   IMODE = 1. Simply does l=0 temperature harmonic old theta term   C
C  (no inhomogeneous part to monopole term) to new poloidal vel.     C
C  This is eq B.52 in thesis first term.                             C
C                                                                    C
C   IMODE = 2. Does v_old . grad ( theta_new ) term for the          C
C  monopole temperature harmonic.                                    C
C  This is eq B.52 in thesis first term.                             C
C                                                                    C
C   IMODE = 3. Does v_new . grad ( theta_old + fg ) term for the     C
C  Q_{ab}^g term ( S_{qq}^{abg} in thesis - this is B.52 first term) C
C                                                                    C
C   IMODE = 4. Does v_old . grad ( theta_new ) term for the          C
C  Q_{ab}^g term ( S_{qq}^{abg} in thesis - this is B.52 first term) C
C                                                                    C
C   IMODE = 5. Does v_new . grad ( theta_old + fg ) term for the     C
C  S_{ab}^g term ( S_{ss}^{abg} in thesis - this is B.52 2nd term)   C
C                                                                    C
C   IMODE = 6. Does v_old . grad ( theta_new ) term for the          C
C  S_{ab}^g term ( S_{ss}^{abg} in thesis - this is B.52 2nd term)   C
C                                                                    C
C   IMODE = 7. Does v_new . grad ( theta_old + fg ) term for the     C
C  T_{ab}^g term ( S_{ts}^{abg} in thesis - this is equation B.53    C
C                                                                    C
C   IMODE = 8. Does v_old . grad ( theta_new ) term for the          C
C  T_{ab}^g term ( S_{ts}^{abg} in thesis - this is equation B.53    C
C                                                                    C
C  INTPAR( 2 ) = LA                                                  C
C  INTPAR( 3 ) = LB                                                  C
C  INTPAR( 4 ) = NR - number of radial grid nodes.                   C
C  INTPAR( 5 ) = NH - number of harmonics in vector, VEC.            C
C  INTPAR( 6 ) = IH - Harmonic number in vector, VEC.                C
C  INTPAR( 7 ) = NDIM - Dimension of VEC; (not neceassarily          C
C                  the same dimension as the matrix in a stability   C
C                     calculation.                                   C
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
C     VEC       : Array (*) - vector containing v_0 and theta_0      C
C                                                                    C
C                                                                    C
C                                                                    C
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
C
C*********************************************************************
      SUBROUTINE VDGTTS ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD,
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
      INTEGER IMODE, LA, LB, NR, NH, IH, NDIM, IND
      DOUBLE PRECISION FAC, SQRLL1, D0F, D1F, D2F, D3F, D4F
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IMODE  = INTPAR( 1 )
      LA     = INTPAR( 2 )
      LB     = INTPAR( 3 )
      NR     = INTPAR( 4 )
      NH     = INTPAR( 5 )
      IH     = INTPAR( 6 )
      NDIM   = INTPAR( 7 )
C
      FAC    = DPPARS( 1 )
C
C Check on the dimensions of VEC
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine VDGTTS. Bad dimensions '
         PRINT *,' Array length = ',NDIM
         PRINT *,' NR = ', NR,' NH = ',NH
         STOP
      ENDIF
C
      IF ( IMODE.EQ.1 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         D0FAC = D1F*LA*(LA+1.0d0)/RAD
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF

      IF ( IMODE.EQ.2 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         D0FAC = 0.0d0
         D1FAC = D0F*LA*(LA+1.0d0)/RAD
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF

      IF ( IMODE.EQ.3 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         D0FAC = FAC*D1F*LA*(LA+1.0d0)/RAD
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF

      IF ( IMODE.EQ.4 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         D0FAC = 0.0d0
         D1FAC = FAC*D0F*LA*(LA+1.0d0)/RAD
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF

      IF ( IMODE.EQ.5 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         D0FAC = FAC*D0F*SQRLL1(LA)*SQRLL1(LB)/(RAD*RAD)
         D1FAC = FAC*D0F*SQRLL1(LA)*SQRLL1(LB)/RAD
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF

      IF ( IMODE.EQ.6 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         D0FAC = FAC*(D1F+D0F/RAD)*SQRLL1(LA)*SQRLL1(LB)/RAD
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF

      IF ( IMODE.EQ.7 ) THEN
         IND = ( IRN - 1 )*NH + IH
         D0F = VEC( IND )
         D0FAC = (-1.0d0)*FAC*D0F*SQRLL1(LA)*SQRLL1(LB)/RAD
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF

      IF ( IMODE.EQ.8 ) THEN
         IND = ( IRN - 1 )*NH + IH
         D0F = VEC( IND )
         D0FAC = (-1.0d0)*FAC*D0F*SQRLL1(LA)*SQRLL1(LB)/RAD
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      RETURN
      END
C*********************************************************************
