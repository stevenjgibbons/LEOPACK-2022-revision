C*********************************************************************
C function PMM1 ******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Evaluates the Schmidt Normalised Legendre Function P_(m+1)^m (X)   C
C according to equation 179 in my notes ; i.e.                       C
C                                                                    C
C    P_(m+1)^m (X) = SQRT( 2m+1 ).X.P^m_m(X)                         C
C                                                                    C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     PMM0      : P_m^m(X) as evaluated by Function PMM.             C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMM1 ( M, X, PMM0 )
      IMPLICIT NONE
      DOUBLE PRECISION PMM1
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION X, PMM0
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = DBLE( M )
      PMM1 = X*PMM0*DSQRT( 2.0d0*RM+1.0d0 )
      RETURN
      END
C*********************************************************************
