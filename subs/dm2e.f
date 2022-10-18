C*********************************************************************
C subroutine DM2E                                                    C
C            ----                                                    C
C Steve Gibbons Wed Feb 28 14:29:43 MET 2001                         C
C                                                                    C
C This is a Dave Gubbins routine which I am converting               C
C to double precision and adding safety checks.                      C
C____________________________________________________________________C
C                                                                    C
C  Converts a D and M value to EPS0, EPS1, EPS2 and EPS3 for the     C
C Kumar Roberts dynamo problem. EPS2 and EPS3 both set to 3.0d0      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
        SUBROUTINE DM2E( D, M, EPS0, EPS1, EPS2, EPS3 )
        IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
        DOUBLE PRECISION D, M, EPS0, EPS1, EPS2, EPS3
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
        DOUBLE PRECISION ALPHA, BETA, GAMMA, DELTA, DZERO,
     1                   FAC1, FAC2, FAC3, FAC4, FAC5, FAC6, DTOL
        INTEGER          SIGND, SIGNM
C
        PARAMETER ( ALPHA = 16.0D0/315.0D0,
     1              BETA  = 248832.0D0/26558675.0D0,
     2              GAMMA = 0.258985312D0,
     3              DELTA = 0.2686211516D0,
     4              DTOL  = 1.0d-9, DZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
        EPS2 = 3.0d0
        EPS3 = 3.0d0
C
        SIGND = 1
        IF ( D.LT.0 ) SIGND = -1
        SIGNM = 1
        IF ( M.LT.0 ) SIGNM = -1
C
        FAC1 = GAMMA*D*SIGNM*EPS2*EPS2 + DELTA*D*SIGNM*EPS3*EPS3
        FAC2 = ALPHA*(-SIGNM*D-SIGND*M+SIGNM*SIGND)
        FAC3 = GAMMA*M*SIGND*EPS2*EPS2 + DELTA*M*SIGND*EPS3*EPS3
        FAC4 = BETA*(-SIGND*M-SIGNM*D+SIGNM*SIGND)
C
        IF ( DABS( FAC2 ).LT.DTOL .OR. DABS( FAC4 ).LT.DTOL ) THEN
          PRINT *,' Subroutine DM2E'
          PRINT *,' Division by zero imminent.'
          PRINT *,' D = ', D,' M = ', M
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
        FAC5 = FAC1/FAC2
        FAC6 = FAC3/FAC4
C
        IF ( FAC5.LT.DZERO .OR. FAC6.LT.DZERO ) THEN
          PRINT *,' Subroutine DM2E'
          PRINT *,' Root of negative number imminent.'
          PRINT *,' D = ', D,' M = ', M
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
        EPS0 = SIGND*DSQRT(FAC5)
        EPS1 = SIGNM*DSQRT(FAC6)
C
        RETURN
        END
C*********************************************************************
