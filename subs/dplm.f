C*********************************************************************
C function DPLM *******************************************************
C Steve Gibbons 18.4.97                                              C
C____________________________________________________________________C
C                                                                    C
C Calculates general P_l^m derivative from eqn 187 in my notes       C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C     L         : Well it's L isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     S         : Sin(theta)                                         C
C     PLMIN1    : P_(l-1)^m ( X )                                    C
C     DPLMN1   : P_(l-1)^m ( X ) derivative                         C
C     DPLMN2   : P_(l-2)^m ( X ) derivative                         C
C____________________________________________________________________C
C Local Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     RL        : DBLE  ( L )                                        C
C     RM        : DBLE  ( M )                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPLM ( L, M, X, S, PLMIN1, DPLMN1, DPLMN2 )
      IMPLICIT NONE
      DOUBLE PRECISION DPLM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M,L
      DOUBLE PRECISION X, S, PLMIN1, DPLMN1, DPLMN2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM,RL
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function DPLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function DPLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' DPLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = DBLE( M )
      RL = DBLE( L )
      DPLM = ( 2.0d0*RL - 1.0d0 )*(X*DPLMN1-S*PLMIN1)
      DPLM = DPLM - DPLMN2*DSQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      DPLM = DPLM/DSQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
