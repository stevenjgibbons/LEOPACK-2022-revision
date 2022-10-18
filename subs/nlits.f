C*********************************************************************
C subroutine Non Linear Inertial Term Subroutine *********************
C            -   -      -        -    -          *********************
C Steve Gibbons  30.06.99                                            C
C____________________________________________________________________C
C                                                                    C
C  DPPARS( 1 ) = Scalar spherical harmonic interaction ( c.f.        C
C                 equations B.39, B.40 and B.41 in my thesis )       C
C  INTPAR( 1 ) = IMODE                                               C
C                                                                    C
C    IMODE = 1. Inputs - Old toroidal alpha                          C
C                        New poloidal beta                           C
C               Output -     toroidal gamma                          C
C     This follows equation ( B.69 ) in my thesis.                   C
C                                                                    C
C    IMODE = 2. Inputs - Old toroidal alpha                          C
C                        New toroidal beta                           C
C               Output -     toroidal gamma                          C
C     This is the term from C_{tqs}^{abg} interaction and is         C
C     the first term of ( B.72 ) in my thesis.                       C
C                                                                    C
C    IMODE = 3. Inputs - Old toroidal alpha                          C
C                        New toroidal beta                           C
C               Output -     poloidal gamma                          C
C     This is the term from C_{tqt}^{abg} interaction and is         C
C     the first term of ( B.71 ) in my thesis.                       C
C                                                                    C
C    IMODE = 4. Inputs - Old toroidal alpha                          C
C                        New toroidal beta                           C
C               Output -     toroidal gamma                          C
C     This is the term from C_{tsq}^{abg} interaction and is         C
C     the second term of ( B.72 ) in my thesis.                      C
C                                                                    C
C    IMODE = 5. Inputs - Old poloidal alpha                          C
C                        New poloidal beta                           C
C               Output -     toroidal gamma                          C
C     This is the term from C_{tsq}^{abg} interaction and is         C
C     the first term of ( B.64 ) in my thesis.                       C
C                                                                    C
C    IMODE = 6. Inputs - Old poloidal alpha                          C
C                        New poloidal beta                           C
C               Output -     poloidal gamma                          C
C     This is the term from C_{qtt}^{abg} interaction and is         C
C     ( B.63 ) in my thesis.                                         C
C                                                                    C
C    IMODE = 7. Inputs - Old poloidal alpha                          C
C                        New toroidal beta                           C
C               Output -     toroidal gamma                          C
C     This is the term from C_{qss}^{abg} interaction and is         C
C     the second term of ( B.67 ) in my thesis.                      C
C                                                                    C
C    IMODE = 8. Inputs - Old poloidal alpha                          C
C                        New toroidal beta                           C
C               Output -     poloidal gamma                          C
C     This is the term from C_{qst}^{abg} interaction and is         C
C     the first term of ( B.66 ) in my thesis.                       C
C                                                                    C
C    IMODE = 9. Inputs - Old poloidal alpha                          C
C                        New poloidal beta                           C
C               Output -     toroidal gamma                          C
C     This is the term from C_{stq}^{abg} interaction and is         C
C     the second term of ( B.64 ) in my thesis.                      C
C                                                                    C
C    IMODE = 10. Inputs - Old poloidal alpha                         C
C                        New toroidal beta                           C
C               Output -     toroidal gamma                          C
C     This is the term from C_{sqs}^{abg} interaction and is         C
C     the third term of ( B.67 ) in my thesis.                       C
C                                                                    C
C    IMODE = 11. Inputs - Old poloidal alpha                         C
C                         New toroidal beta                          C
C               Output -     poloidal gamma                          C
C     This is the term from C_{sqt}^{abg} interaction and is         C
C     the second term of ( B.66 ) in my thesis.                      C
C                                                                    C
C    IMODE = 12. Inputs - Old poloidal alpha                         C
C                         New toroidal beta                          C
C               Output -     toroidal gamma                          C
C     This is the term from C_{ssq}^{abg} interaction and is         C
C     the first term of ( B.67 ) in my thesis.                       C
C                                                                    C
C    IMODE = 13. Inputs - New toroidal alpha                         C
C                         Old poloidal beta                          C
C               Output -     toroidal gamma                          C
C     This is the term from C_{ttq}^{abg} interaction and is         C
C     the first term of ( B.69 ) in my thesis.                       C
C                                                                    C
C    IMODE = 14. Inputs - New toroidal alpha                         C
C                         Old toroidal beta                          C
C               Output -     toroidal gamma                          C
C     This is the term from C_{tqs}^{abg} interaction and is         C
C     the first term of ( B.72 ) in my thesis.                       C
C                                                                    C
C    IMODE = 15. Inputs - New toroidal alpha                         C
C                         Old toroidal beta                          C
C               Output -     poloidal gamma                          C
C     This is the term from C_{tqt}^{abg} interaction and is         C
C     the first term of ( B.71 ) in my thesis.                       C
C                                                                    C
C    IMODE = 16. Inputs - New toroidal alpha                         C
C                         Old toroidal beta                          C
C               Output -     toroidal gamma                          C
C     This is the term from C_{tsq}^{abg} interaction and is         C
C     the second term of ( B.72 ) in my thesis.                      C
C                                                                    C
C    IMODE = 17. Inputs - New poloidal alpha                         C
C                         Old poloidal beta                          C
C               Output -     toroidal gamma                          C
C     This is the term from C_{qts}^{abg} interaction and is         C
C     the first term of ( B.64 ) in my thesis.                       C
C                                                                    C
C    IMODE = 18. Inputs - New poloidal alpha                         C
C                         Old poloidal beta                          C
C               Output -     poloidal gamma                          C
C     This is the term from C_{qtt}^{abg} interaction and is         C
C     the first term of ( B.63 ) in my thesis.                       C
C                                                                    C
C    IMODE = 19. Inputs - New poloidal alpha                         C
C                         Old toroidal beta                          C
C               Output -     toroidal gamma                          C
C     This is the term from C_{qss}^{abg} interaction and is         C
C     the second term of ( B.67 ) in my thesis.                      C
C                                                                    C
C    IMODE = 20. Inputs - New poloidal alpha                         C
C                         Old toroidal beta                          C
C               Output -     poloidal gamma                          C
C     This is the term from C_{qst}^{abg} interaction and is         C
C     the first term of ( B.66 ) in my thesis.                       C
C                                                                    C
C    IMODE = 21. Inputs - New poloidal alpha                         C
C                         Old poloidal beta                          C
C               Output -     toroidal gamma                          C
C     This is the term from C_{stq}^{abg} interaction and is         C
C     the second term of ( B.64 ) in my thesis.                      C
C                                                                    C
C    IMODE = 22. Inputs - New poloidal alpha                         C
C                         Old toroidal beta                          C
C               Output -     toroidal gamma                          C
C     This is the term from C_{sqs}^{abg} interaction and is         C
C     the third term of ( B.67 ) in my thesis.                       C
C                                                                    C
C    IMODE = 23. Inputs - New poloidal alpha                         C
C                         Old toroidal beta                          C
C               Output -     poloidal gamma                          C
C     This is the term from C_{sqt}^{abg} interaction and is         C
C     the second term of ( B.66 ) in my thesis.                      C
C                                                                    C
C    IMODE = 24. Inputs - New poloidal alpha                         C
C                         Old toroidal beta                          C
C               Output -     toroidal gamma                          C
C     This is the term from C_{ssq}^{abg} interaction and is         C
C     the first term of ( B.67 ) in my thesis.                       C
C                                                                    C
C  INTPAR( 2 ) = LA                                                  C
C  INTPAR( 3 ) = LB                                                  C
C  INTPAR( 4 ) = LG                                                  C
C  INTPAR( 5 ) = NR - number of radial grid nodes.                   C
C  INTPAR( 6 ) = NH - number of harmonics in vector, VEC.            C
C  INTPAR( 7 ) = IH - Harmonic number in vector, VEC.                C
C  INTPAR( 8 ) = NDIM - Dimension of VEC; (not neceassarily          C
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
      SUBROUTINE NLITS ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD,
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
      INTEGER IMODE, LA, LB, LG, NR, NH, IH, NDIM, IND
      DOUBLE PRECISION FLA, FLB, FLG, FAC, SQRLL1,
     1                 D0F, D1F, D2F, D3F, D4F,
     2                 RAD2, RAD3, RAD4, TEMP
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IMODE  = INTPAR( 1 )
      LA     = INTPAR( 2 )
      LB     = INTPAR( 3 )
      LG     = INTPAR( 4 )
      NR     = INTPAR( 5 )
      NH     = INTPAR( 6 )
      IH     = INTPAR( 7 )
      NDIM   = INTPAR( 8 )
C
      FAC    = DPPARS( 1 )
C
      FLA = SQRLL1( LA )
      FLB = SQRLL1( LB )
      FLG = SQRLL1( LG )
C
      RAD2 = RAD*RAD
      RAD3 = RAD*RAD*RAD
      RAD4 = RAD*RAD*RAD*RAD
C
C Check on the dimensions of VEC
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine NLITS. Bad dimensions '
         PRINT *,' Array length = ',NDIM
         PRINT *,' NR = ', NR,' NH = ',NH
         STOP
      ENDIF
C
      IF ( IMODE.EQ.1 ) THEN
         IND = ( IRN - 1 )*NH + IH
         D0F = VEC( IND )
         TEMP = D0F * FLA * FLB * FAC / RAD
         D0FAC = (-1.0d0)*TEMP*FLB*FLB/(RAD*RAD)
         D1FAC = (2.0d0)*TEMP/RAD
         D2FAC = TEMP
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.2 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = (-1.0d0)*FAC*FLA*FLB*FLB/FLG
         D0FAC = TEMP*D1F/RAD
         D1FAC = TEMP*D0F/RAD
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.3 ) THEN
         IND = ( IRN - 1 )*NH + IH
         D0F = VEC( IND )
         TEMP = (-1.0d0)*FAC*FLA*FLB*FLB/FLG
         D0FAC = TEMP*D0F/RAD
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.4 ) THEN
         IND = ( IRN - 1 )*NH + IH
         D0F = VEC( IND )
         TEMP = FAC*FLA*FLB/RAD
         D0FAC = TEMP*D0F/RAD
         D1FAC = TEMP*D0F
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.5 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLA*FLB/FLG
         D0FAC = TEMP*( 2.0d0*(LB*LB + LB)*D0F/RAD4 -
     1              D1F*(LB*LB + LB)/RAD3 )
         D1FAC = TEMP*( 2.0d0*D1F/RAD2 - 2.0d0*D0F/RAD3 -
     1              D0F*(LB*LB + LB)/RAD3 )
         D2FAC = TEMP*( 2.0d0*D0F/RAD2 + D1F/RAD )
         D3FAC = TEMP*D0F/RAD
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.6 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLA*FLB*D0F/(FLG*RAD)
         D0FAC = TEMP*(-1.0d0)*FLB*FLB/RAD2
         D1FAC = TEMP*(2.0d0)/RAD
         D2FAC = TEMP
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.7 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLA*FLB/FLG
         D0FAC = TEMP*( D1F/RAD2 - D0F/RAD3 )
         D1FAC = TEMP*( D0F/RAD2 + D1F/RAD )
         D2FAC = TEMP*D0F/RAD
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.8 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLA*FLB*D0F/(FLG*RAD)
         D0FAC = TEMP/RAD
         D1FAC = TEMP
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.9 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = (-1.0d0)*FAC*FLA*FLB
         D0FAC = TEMP*(-1.0d0)*FLB*FLB*( D0F/RAD4 + D1F/RAD3 )
         D1FAC = TEMP*(2.0d0)*(D0F/RAD3 + D1F/RAD2)
         D2FAC = TEMP*(D0F/RAD2 + D1F/RAD)
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.10 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLB*FLB/FLG
         D0FAC = TEMP*( D2F/RAD + D1F/RAD2 - D0F/RAD3 )
         D1FAC = TEMP*( D0F/RAD2 + D1F/RAD )
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.11 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLB*FLB/FLG
         D0FAC = TEMP*( D1F + D0F/RAD )/RAD
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.12 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = (-1.0d0)*FAC*FLA*FLB
         D0FAC = TEMP*( D0F/RAD3 + D1F/RAD2 )
         D1FAC = TEMP*( D0F/RAD2 + D1F/RAD )
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.13 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLB/RAD
         D0FAC = TEMP*( D2F + 2.0d0*D1F/RAD - FLB*FLB*D0F/RAD2 )
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.14 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = (-1.0d0)*FAC*FLA*FLB*FLB/(RAD*FLG)
         D0FAC = TEMP*D1F
         D1FAC = TEMP*D0F
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.15 ) THEN
         IND = ( IRN - 1 )*NH + IH
         D0F = VEC( IND )
         TEMP = (-1.0d0)*FAC*FLA*FLB*FLB/(RAD*FLG)
         D0FAC = TEMP*D0F
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.16 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLB/RAD
         D0FAC = TEMP*( D0F/RAD + D1F )
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.17 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLA*FLB/FLG
         D0FAC = TEMP*( 2.0d0*(LB*LB + LB)*D0F/RAD4
     1                  - (LB*LB + LB)*D1F/RAD3
     2                  - 2.0d0*D1F/RAD3
     3                  + 2.0d0*D2F/RAD2
     4                  + D3F/RAD )
         D1FAC = TEMP*( D2F/RAD + 2.0d0*D1F/RAD2
     1                  - (LB*LB + LB)*D0F/RAD3 )
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.18 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLA*FLB/FLG
         D0FAC = TEMP*( D2F/RAD + 2.0d0*D1F/RAD2
     1                  - D0F*( LB*LB + LB )/RAD3 )
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.19 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLA*FLB/FLG
         D0FAC = TEMP*( D2F/RAD - D0F/RAD3 + D1F/RAD2 )
         D1FAC = TEMP*( D0F/RAD2 + D1F/RAD )
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.20 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLA*FLB/FLG
         D0FAC = TEMP*( D0F/(RAD*RAD) + D1F/RAD )
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.21 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = (-1.0d0)*FAC*FLA*FLB/RAD
         D0FAC = TEMP*( D2F/RAD + 2.0d0*D1F/RAD2 -
     1               DBLE( LB*LB + LB )*D0F/RAD3 )
         D1FAC = TEMP*( D2F + 2.0d0*D1F/RAD -
     1               DBLE( LB*LB + LB )*D0F/RAD2 )
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.22 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = FAC*FLA*FLB*FLB/FLG
         D0FAC = TEMP*( D1F/RAD2 - D0F/RAD3 )
         D1FAC = TEMP*( D1F/RAD + D0F/RAD2 )
         D2FAC = TEMP*( D0F/RAD )
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.23 ) THEN
         IND = ( IRN - 1 )*NH + IH
         D0F = VEC( IND )
         TEMP = FAC*FLA*FLB*FLB*D0F/FLG
         D0FAC = TEMP/RAD2
         D1FAC = TEMP/RAD
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
      IF ( IMODE.EQ.24 ) THEN
         CALL SVDERF ( NDIM, NH, IH, NR, IRN, VEC, ORD,
     1                 D0F, D1F, D2F, D3F, D4F, RI, RO )
         TEMP = (-1.0d0)*FAC*FLA*FLB
         D0FAC = TEMP*( D0F/RAD3 + D1F/RAD2 )
         D1FAC = TEMP*( D0F/RAD2 + D1F/RAD )
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
      ENDIF
C
C
      RETURN
      END
C*********************************************************************
