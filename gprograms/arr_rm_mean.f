C*********************************************************************
C Steve Gibbons 15 January 2004
C If RMFLAG is not equal to 1, return with no action.
C If RMFLAG is equal to 1, then calculate the mean of the real
C array and subtract it from each element.
C We also return with no action if N1 or N2 is less than 1.
C--------------------------------------------------
      SUBROUTINE ARR_RM_MEAN( N1, N2, ARR, RMFLAG )
      IMPLICIT NONE
      INTEGER        RMFLAG, N1, N2
      REAL           ARR( N1, N2 )
C
      INTEGER        I1, I2
      REAL           RTOTEL, RMEANV
C
      IF ( RMFLAG.NE.1 ) RETURN
      IF ( N1.LT.1 )     RETURN
      IF ( N2.LT.1 )     RETURN
C
      RTOTEL = 1.0/FLOAT( N1*N2 )
      RMEANV = 0.0
      DO I2 = 1, N2
        DO I1 = 1, N1
          RMEANV = RMEANV + ARR( I1, I2 )*RTOTEL
        ENDDO
      ENDDO
C
C Have now calculated the mean - now subtract it ...
C
      DO I2 = 1, N2
        DO I1 = 1, N1
          ARR( I1, I2 ) = ARR( I1, I2 ) - RMEANV
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
