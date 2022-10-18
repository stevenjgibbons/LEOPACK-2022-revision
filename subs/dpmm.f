C*********************************************************************
C function DPMM ******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Calculates the derivative of P_m^m(theta) according to equation    C
C 185 in my notes i.e.                                               C
C                                                                    C
C d P_m^m(theta)/ d(theta) = sqrt(2.0)*M*cos(theta)/sin(theta) * A   C
C  with A = product_(i=1)^m ( SQRT( (2i-1)/2i ) * sin(theta) )       C
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
      FUNCTION DPMM ( M , C , S )
      IMPLICIT NONE
      DOUBLE PRECISION DPMM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION C, S
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
      DOUBLE PRECISION TWOI
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function DPMM. M < 0 error. Stopped.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         DPMM = 0.0d0
         RETURN
      ENDIF
C ................................. so M is greater than 0
      DPMM = DSQRT ( 2.0d0 )*DBLE(M)*C/S
      DO I = 1, M
         TWOI = 2.0d0*DBLE(I)
         DPMM = DPMM * S * DSQRT ( (TWOI - 1.0d0 )/TWOI )
      ENDDO
      RETURN
      END
C*********************************************************************
