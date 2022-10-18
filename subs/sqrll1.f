C*********************************************************************
C function SQuare Root of L*(L+1) ************************************
C          --     -       -  - -  ************************************
C Steve Gibbons 25.4.97                                              C
C____________________________________________________________________C
C
      FUNCTION SQRLL1 ( L )
      IMPLICIT NONE
      DOUBLE PRECISION SQRLL1,Q
      INTEGER L
   
      IF ( L.LT.0 ) THEN
         PRINT *,' Function SQRLL1. L less than 0. Program aborted.'
         STOP
      ENDIF
      Q = DBLE( L )
      SQRLL1 = DSQRT( Q*Q + Q )
      RETURN
      END
C*********************************************************************
