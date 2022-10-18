C*********************************************************************
C subroutine Oscillatory Function Evaluate (with OPTION) *************
C            ---         ----     ----           ------  *************
C Steve Gibbons Tue Feb 27 11:02:07 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C F, FREAL and FIMAG are real arrays of dimensions ( N1, N2 ).       C
C FREAL has been evaluated from the Real part of a pair of           C
C eigenvectors and FIMAG is the imaginary part.                      C
C                                                                    C
C If we wish to look at a total of NPOC points in time on a half or  C
C full cycle, then OSC_FUNC_EVAL evaluates, in the array F, the      C
C function at time IPOC.                                             C
C                                                                    C
C To look at a half-cycle, set ICYCLE = 1                            C
C To look at a full-cycle, set ICYCLE = 2                            C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N1        : First dimension of arrays F, FREAL and FIMAG       C
C     N2        : Second dimension of arrays F, FREAL and FIMAG      C
C                                                                    C
C     IPOC      : Number of "time-step" on the half cycle.           C
C     NPOC      : Total number of "time-steps" on the half cycle.    C
C                                                                    C
C     ICYCLE    : 1 --> half-cycle, 2 --> full cycle                 C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     FREAL     : Real part of function. Dim ( N1, N2 )              C
C     FIMAG     : Imag part of function. Dim ( N1, N2 )              C
C                                                                    C
C     F         : Function eval.ted at time IPOC. Dim ( N1, N2 )     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE OSC_FUNC_EVAL_OPTION( N1, N2, IPOC, NPOC, FREAL,
     1                                 FIMAG, F, ICYCLE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, IPOC, NPOC, ICYCLE
      REAL    FREAL( N1, N2 ), FIMAG( N1, N2 ), F( N1, N2 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I1, I2
      REAL    PI, TAU, CT, ST
      PARAMETER ( PI=3.1415926 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( ICYCLE.NE.1 .AND. ICYCLE.NE.2 ) THEN
        PRINT *,' Subroutine OSC_FUNC_EVAL_OPTION.'
        PRINT *,' ICYCLE = ', ICYCLE
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      TAU = REAL( ICYCLE*IPOC )*PI/REAL( NPOC )
      CT  = COS( TAU )
      ST  = SIN( TAU )
C
      DO I2 = 1, N2
        DO I1 = 1, N1
          F( I1, I2 ) = FREAL( I1, I2 )*CT - FIMAG( I1, I2 )*ST
        ENDDO
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
