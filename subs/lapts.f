C*********************************************************************
C subroutine LAPlacian Term  Subroutine ******************************
C            ---       -     -          ******************************
C Steve Gibbons  20.11.98                                            C
C____________________________________________________________________C
C Fills in the terms for the Laplacian when called by for example    C
C dpmes. Can currently return 1.0d0 only ( INTPAR( 2 )=0 ),          C
C  D_l ( INTPAR( 2 )=1 ) or -D_l^2 ( INTPAR( 2 )=2 ).                C
C Most of the arguments are not required but the form is             C
C maintained for generality and compatibility.                       C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     INTPAR    : Array of dimension (*)                             C
C         INTPAR( 1 ) must contain the spherical harmonic degree, L  C
C         INTPAR( 2 ) must contain the operation flag, IOP.          C
C          IF IOP = 0,  d0fac is returned as 1.0d0, all deriv.s zero C
C          IF IOP = 1,  D_l is returned.                             C
C          IF IOP = 2,  -D_l^2 is returned.                          C
C     IRN       : Number of the radial grid node.                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Value of the radius.                               C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     VEC       : Array (*) - not referred to.                       C
C     DPPARS    : Array (*) - not referred to.                       C
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
C Note IRN, RI, RO, VEC, DPPARS and ORD or not referred to by LAPTS  C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LAPTS ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD,
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
      INTEGER L, IOP
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      L = INTPAR( 1 )
      IOP = INTPAR( 2 )
C
      IF ( IOP.EQ.0 ) THEN
         D0FAC = 1.0d0
         D1FAC = 0.0d0
         D2FAC = 0.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
         RETURN
      ENDIF
C
      IF ( IOP.EQ.1 ) THEN
         D0FAC = (-1.0d0)*DBLE(L*L+L)/(RAD*RAD)
         D1FAC = 2.0d0/RAD
         D2FAC = 1.0d0
         D3FAC = 0.0d0
         D4FAC = 0.0d0
         RETURN
      ENDIF
C
      IF ( IOP.EQ.2 ) THEN
         D0FAC = (L+2.0d0)*(1.0d0-L)*(L*L+L)/(RAD*RAD*RAD*RAD)
         D1FAC = 0.0d0
         D2FAC = (2.0d0)*DBLE(L*L+L)/(RAD*RAD)
         D3FAC = (-4.0d0)/RAD
         D4FAC = (-1.0d0)
         RETURN
      ENDIF
C
      PRINT *,' Subroutine LAPTS. INTPAR( 2 ) = ', INTPAR( 2 )
      PRINT *,' Current options are 0, 1 or 2.'
      STOP
      END
C*********************************************************************
