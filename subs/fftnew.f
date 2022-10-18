C*********************************************************************
C subroutine Fast Fourier Transform NEW version **********************
C            -    -       -         ---         **********************
C Steve Gibbons 22.4.97 (Based on FOUR1 from Numerical Recipes       C
C      I've made slight alterations for double precision and error   C
C      checking at the start.)					     C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NN	: Number of data to be transformed.		     C
C                  NN MUST be a power of two. This is checked for by C
C                  the subroutine POWTWO.                            C
C     ISIGN	: ISIGN = 1. DATA is replaced by it's discrete       C
C                   Fourier Transform.                               C
C                 ISIGN = -1. DATA is replaced by NN * it's inverse  C
C                   discrete Fourier Transform.                      C
C  Double Precision                                                  C
C  ----------------                                                  C
C     DATA	: Array of dimension (2*NN).                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FFTNEW ( DATA, NN, ISIGN )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NN, ISIGN
      DOUBLE PRECISION DATA( 2 * NN )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TEMPR,TEMPI,THETA,WPR,WPI,WR,WI,WTEMP
      INTEGER I,J,N,M,MMAX,ISTEP
      LOGICAL POWT
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C 
C First check that the number NN supplied is infact a power of 2.
      CALL POWTWO ( NN, POWT)
      IF ( .NOT. POWT) THEN
         PRINT *,' Subroutine FFTNEW. NN is not a power of 2.'
         PRINT *,'Program aborted.'
         STOP
      ENDIF
C Now check that ISIGN = 1 or -1.
      IF ( ISIGN*ISIGN .NE. 1) THEN
         PRINT *,' Subroutine FFTNEW. ISIGN must be 1 or -1.'
         PRINT *,'Program aborted.'
         STOP
      ENDIF
C____________________________________________________________________C
C
      N = 2 * NN
      J = 1
      DO I = 1, N, 2
         IF ( J.GT.I ) THEN
            TEMPR = DATA( J )
            TEMPI = DATA( J+1 )
            DATA( J ) = DATA( I )
            DATA( J+1 ) = DATA( I+1 )
            DATA( I ) = TEMPR
            DATA( I+1 ) = TEMPI
         ENDIF
C ... now calculate inverse bit map for next value of I.
         M = N/2
 500     IF ( M.GE.2 .AND. J.GT.M ) THEN
            J = J-M
            M = M/2
            GOTO 500
         ENDIF
         J = J + M
      ENDDO
C ............. now do the real part of transform.
      MMAX = 2
 502  IF ( N.GT.MMAX ) THEN
         ISTEP = 2*MMAX
         THETA = 2.0d0*PI/DBLE(ISIGN*MMAX)
         WPR = -2.0d0*DSIN( 0.5d0*THETA)**2
         WPI = DSIN(THETA)
         WR = 1.0d0
         WI = 0.0d0
         DO M = 1, MMAX, 2
            DO I = M, N, ISTEP
               J = I + MMAX
               TEMPR = WR*DATA( J ) - WI*DATA( J + 1 )
               TEMPI = WR*DATA( J + 1 ) + WI*DATA( J )
               DATA( J ) = DATA( I ) - TEMPR
               DATA( J + 1 ) = DATA( I + 1 ) - TEMPI
               DATA( I ) = DATA( I ) + TEMPR
               DATA( I + 1 ) = DATA( I + 1 ) + TEMPI
            ENDDO
            WTEMP = WR
            WR = WR*WPR - WI*WPI + WR
            WI = WI*WPR + WTEMP*WPI + WI
         ENDDO
         MMAX = ISTEP
         GOTO 502
      ENDIF

      RETURN
      END
C*********************************************************************

