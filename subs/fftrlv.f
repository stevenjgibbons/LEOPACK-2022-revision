C*********************************************************************
C subroutine Fast Fourier Transform ReaL Version *********************
C            -    -       -         -  - -       *********************
C Steve Gibbons 22.4.97 (Based on FOUR1 from Numerical Recipes       C
C      I've made slight alterations for double precision, error      C
C      checking at the start and rescaling the coefficients after    C
C      the forward transform.)                                       C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NN	: Number of data to be transformed.		     C
C                  NN MUST be a power of two. This is checked for by C
C                  the subroutine POWTWO.                            C
C                                                                    C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       C
C                                                                    C
C    The real function f( PHI ) must be periodic, satisfying         C
C    f( PHI ) = f( PHI + 2*PI )                                      C
C                                                                    C
C    f( PHI ) has the formulation                                    C
C                                                                    C
C    f( PHI ) = sum_{m=0}^M [ c_m cos(m PHI) + s_m sin( m PHI) ]     C
C                                                                    C
C    FFTRLV will transform between the coefficients {c_m, s_m}       C
C    and discretised values of f at NN equally spaced points PHI_i   C
C    with                                                            C
C                                                                    C
C    PHI_i = (i-1)*DELTA_phi          and                            C
C                                                                    C
C    DELTA_phi = 2*PI/dble( NN )                                     C
C                                                                    C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       C
C                                                                    C
C     ISIGN	:                                                    C
C                                                                    C
C    If ISIGN = 1, on entry DATA( 2*I - 1 ) must contain the value   C
C    of f at the point PHI = PHI_i as defined above.                 C
C                                                                    C
C    On exit, DATA( 2*m + 1 ) will contain the coeff. c_m and        C
C             DATA( 2*m + 2 ) will contain the coeff. s_m and        C
C                                                                    C
C    for m = 0, M.                                                   C
C      For accuracy, M **MUST** be less than NN/2                    C
C    (This is the Nyquist criterion - see Numerical Recipes ).       C
C                                                                    C
C    If ISIGN = -1, DATA must be zero except for the values          C
C    c_m contained in DATA( 2*m + 1 )  and                           C
C    s_m contained in DATA( 2*m + 2 ) on input.                      C
C                                                                    C
C    On output, f( PHI_i ) will be contained in DATA( 2*I - 1 ).     C
C                                                                    C
C    Note that this routine SCALES the coefficients after            C
C    transforming, which the Numerical Recipes version doesn't.      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     DATA	: Array of dimension (2*NN). See above ...           C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FFTRLV ( DATA, NN, ISIGN )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NN, ISIGN
      DOUBLE PRECISION DATA( 2 * NN )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TEMPR,TEMPI,THETA,WPR,WPI,WR,WI,WTEMP,FAC
      INTEGER I,J,N,M,MMAX,ISTEP
c     LOGICAL POWT
      DOUBLE PRECISION PI, HTHETA
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C 
C  SJG Mon Feb  5 18:18:33 WET 2001
C  I have commented out the power of two checking
C  for speed: reinstate if necessary:
C 
C First check that the number NN supplied is infact a power of 2.
c     CALL POWTWO ( NN, POWT)
c     IF ( .NOT. POWT) THEN
c        PRINT *,' Subroutine FFTRLV. NN is not a power of 2.'
c        PRINT *,'Program aborted.'
c        STOP
c     ENDIF
C Now check that ISIGN = 1 or -1.
      IF ( ISIGN*ISIGN .NE. 1) THEN
         PRINT *,' Subroutine FFTRLV. ISIGN must be 1 or -1.'
         PRINT *,'Program aborted.'
         STOP
      ENDIF
C____________________________________________________________________C
C Here we will zero the imaginary part of the function
C just incase it contains non-zero data which may harm the
C outcome.
C
      IF ( ISIGN.EQ.1 ) THEN
        DO I = 1, NN
          DATA( 2*I ) = 0.0d0
        ENDDO
      ENDIF
C
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
         HTHETA = PI/DBLE(ISIGN*MMAX)
c        THETA = 2.0d0*PI/DBLE(ISIGN*MMAX)
C set THETA temporarily
         THETA = DSIN(HTHETA)
         WPR = -2.0d0*THETA*THETA
c        WPR = -2.0d0*DSIN(HTHETA)**2
         THETA = 2.0d0*HTHETA
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
C
      IF ( ISIGN.EQ.1 ) THEN
        FAC = 2.0d0/DBLE( NN )
        DO I = 1, NN
          DATA( I ) = FAC*DATA( I )
        ENDDO
        DO I = NN+1, 2*NN
          DATA( I ) = 0.0d0
        ENDDO
      ENDIF
C
      RETURN
      END
C*********************************************************************

