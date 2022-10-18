C*********************************************************************
C function ScHMidt normalised Deriv of (Polynomial) LeGendre function*
C          - --               -         -           - -      *********
C Steve Gibbons 8.5.97						     C
C____________________________________________________________________C
C Adapted from SCHNLF to give a single derivative of a               C
C Schmidt Normalised Assoc.                                          C
C Legendre Function as a func. of integers L and M.                  C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     L		: Level of Assoc. Legendre Polynomial.               C
C     M		: Order of Assoc. Legendre Polynomial.               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     COSTH	: Cosine of theta.				     C
C____________________________________________________________________C
C Functions called:-						     C
C ----------------                                                   C
C PMM ( M, S)                                                        C
C PMM1 ( M, X, PMM0 )                                                C
C PLM ( L, M, X, PLMIN1, PLMIN2 )                                    C
C DPMM ( M, C, S)                                                    C
C DPMM1 ( M, C, S, PMM0, DPMM )                                      C
C DPLM ( L, M, C, S, PMMIN1 , DPLMIN1, DPLMIN2 )                     C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION SHDPLG ( L, M, COSTH )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER L,M
      DOUBLE PRECISION SHDPLG, COSTH
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION SINE,DPOLD,DPOLD1,POLD1,POLD,PNEW
      INTEGER I
C____________________________________________________________________C
C Variable declarations - Functions called ..........................C
      DOUBLE PRECISION PMM, PMM1, PLM
      DOUBLE PRECISION DPMM, DPMM1, DPLM
C
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check validity of arguments ....
      IF (L.LT.0) THEN
         PRINT *,' Function SHDPLG. L.LT.0 Program Stopped. '
         STOP
      ENDIF
      IF (M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' Function SHDPLG.'
         PRINT *,' M must be between 0 and L.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( COSTH.LT.-1.0d0 .OR. COSTH.GT.1.0d0 ) THEN
         PRINT *,' Function SHDPLG.'
         PRINT *,' Cos(theta) out of range.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C ....................................................................
      SINE = (1.0d0 - COSTH)*(1.0d0 + COSTH )
      SINE = DSQRT( SINE )
C ....................................................................
C .... first caluclate DPMM
      SHDPLG = DPMM ( M, COSTH, SINE )
      IF ( L.EQ.M ) RETURN
C ...... so L is atleast M + 1 ... so calculate PMM and therefore
C DPMM1 ................
      DPOLD = SHDPLG
      POLD = PMM( M , SINE )
      SHDPLG = DPMM1 ( M, COSTH, SINE, POLD, DPOLD )
      IF ( L.EQ.M+1 ) RETURN
C ..... so L is atleast M + 2 ...
      POLD1 = POLD
      POLD = PMM1( M , COSTH, POLD1 )
      DO I = M + 2, L
         DPOLD1 = DPOLD
         DPOLD = SHDPLG
         SHDPLG = DPLM ( I, M, COSTH, SINE, POLD, DPOLD, DPOLD1 )
         PNEW = PLM ( I , M, COSTH, POLD, POLD1 )
         POLD1 = POLD
         POLD = PNEW
      ENDDO
      
      RETURN
      END
C*********************************************************************
