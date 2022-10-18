C*********************************************************************
C function PMM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Gives the Schmidt Normalised Legendre Function P_m^m ( X )         C
C from eqn. 175 in my notes ie                                       C
C                                                                    C
C                   ( 2m - 1)!! (1- XX)^(m/2) * SQRT (2.0d0 )        C
C   P_m^m( X ) =  ---------------------------------------------      C
C                        SQRT ( (2m)! )                              C
C                                                                    C
C       for m non-zero and                                           C
C                                                                    C
C   P_0^0( X ) = 1.0d0                                               C
C                                                                    C
C N.B. The double factorial sign means the product of all ODD        C
C integers between 1 and ( 2m - 1 ).                                 C
C Best to use the form in eq 184 to calculate PMM.                   C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SINE      : Sin (theta)                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMM ( M, S )
      IMPLICIT NONE
      DOUBLE PRECISION PMM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION S
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
      DOUBLE PRECISION TWOI
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM. M < 0 error. Stopped.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         PMM = 1.0d0
         RETURN
      ENDIF
C ................................. so M is greater than 0
      PMM = DSQRT ( 2.0d0 )
      DO I = 1, M
         TWOI = 2.0d0*DBLE(I)
         PMM = PMM * S * DSQRT ( (TWOI - 1.0d0 )/TWOI )
      ENDDO
      RETURN
      END
C*********************************************************************
