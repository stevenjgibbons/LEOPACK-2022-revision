C*********************************************************************
C subroutine Double precision DERivative Calulate - 3 coefficients ***
C            -                ---        -          -              ***
C Steve Gibbons 15.3.98                                              C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     COEF1     : First array of finite difference coeff.s dim(4)    C
C     F1        :  ... corresponding function value                  C
C     COEF2     : Second array of finite difference coeff.s dim(4)   C
C     F2        :  ... corresponding function value                  C
C     COEF3     : Third array of finite difference coeff.s dim(4)    C
C     F3        :  ... corresponding function value                  C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     D1F       : First derivative.                                  C
C     D2F       : Second derivative.                                 C
C     D3F       : Third derivative.                                  C
C     D4F       : Fourth derivative.                                 C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DDERC3 ( COEF1, F1, COEF2, F2, COEF3, F3,
     1                    D1F, D2F, D3F, D4F )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION COEF1( 4 ), F1, COEF2( 4 ), F2,
     1                 COEF3( 4 ), F3,
     2                 D1F, D2F, D3F, D4F
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      D1F = COEF1( 1 )*F1 + COEF2( 1 )*F2 + COEF3( 1 )*F3
C
      D2F = COEF1( 2 )*F1 + COEF2( 2 )*F2 + COEF3( 2 )*F3
C
      D3F = COEF1( 3 )*F1 + COEF2( 3 )*F2 + COEF3( 3 )*F3
C
      D4F = COEF1( 4 )*F1 + COEF2( 4 )*F2 + COEF3( 4 )*F3
C
      RETURN
      END
C*********************************************************************

