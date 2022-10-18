C*********************************************************************
C function DPMM1 *****************************************************
C Steve Gibbons 18.4.97                                              C
C____________________________________________________________________C
C Calculates d(P_(m+1)^m/d(theta) by equation 185 in my notes        C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     S         : Sin(theta)                                         C
C     PMM0      : P_m^m(X) as evaluated by Function PMM.             C
C     DPMM0     : d(P_m^m(X))/d(theta) as evaluated by Function DPMM C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPMM1 ( M, X, S, PMM0, DPMM0)
      IMPLICIT NONE
      DOUBLE PRECISION DPMM1
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION X, S, PMM0, DPMM0
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function DPMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = DBLE( M )
      DPMM1 = DSQRT( 2.0d0*RM+1.0d0 )*( X*DPMM0 - S*PMM0 )
      RETURN
      END
C*********************************************************************
