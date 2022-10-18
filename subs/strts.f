C*********************************************************************
C subroutine STRatification Term Subroutine **************************
C            ---            -    -          **************************
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C Returns a value to add the stratification term to the matrix when  C
C called by DPMES for example.                                       C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     INTPAR    : Array of dimension (*)                             C
C         INTPAR( 1 ) must contain the spherical harmonic degree, L  C
C     IRN       : Number of the radial grid node.                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Value of the radius.                               C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     VEC       : Array (*) - not referred to.                       C
C     DPPARS    : Array (*) - DPPARS( 1 ) = CB1                      C
C                             DPPARS( 2 ) = CB2                      C
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
C Note IRN, RI, RO, VEC and ORD or not referred to by STRTS          C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE STRTS ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD,
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
      INTEGER L
      DOUBLE PRECISION CB1, CB2
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      L = INTPAR( 1 )
      CB1 = DPPARS( 1 )
      CB2 = DPPARS( 2 )
C
      D0FAC =  DBLE(L*L + L)*( CB1 + CB2/RAD**3 )
      D1FAC =  0.0d0
      D2FAC =  0.0d0
      D3FAC =  0.0d0
      D4FAC =  0.0d0
C
      RETURN
      END
C*********************************************************************
