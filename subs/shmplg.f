C*********************************************************************
C function ScHMidt normalised (Polynomial) LeGendre function *********
C          - --                -           - -               *********
C Steve Gibbons 8.5.97						     C
C____________________________________________________________________C
C Adapted from SCHNLF to give a single Schmidt Normalised Assoc.     C
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
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION SHMPLG ( L, M, COSTH )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER L,M
      DOUBLE PRECISION SHMPLG, COSTH
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION SINE,POLD,POLD1
      INTEGER I
C____________________________________________________________________C
C Variable declarations - Functions called ..........................C
      DOUBLE PRECISION PMM, PMM1, PLM
C
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check validity of arguments ....
      IF (L.LT.0) THEN
         PRINT *,' Function SHMPLG. L.LT.0 Program Stopped. '
         STOP
      ENDIF
      IF (M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' Function SHMPLG.'
         PRINT *,' M must be between 0 and L.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( COSTH.LT.-1.0d0 .OR. COSTH.GT.1.0d0 ) THEN
         PRINT *,' Function SHMPLG.'
         PRINT *,' Cos(theta) out of range.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C ....................................................................
      SINE = (1.0d0 - COSTH)*(1.0d0 + COSTH )
      SINE = DSQRT( SINE )
      SHMPLG = PMM ( M, SINE )
      IF ( L.EQ.M ) RETURN
      POLD = SHMPLG
      SHMPLG = PMM1 ( M, COSTH, POLD )
      IF ( L.EQ.M+1 ) RETURN
C ............................ so L is greater than M + 1.
      DO I = M + 2, L
         POLD1 = POLD
         POLD = SHMPLG
         SHMPLG = PLM( I, M, COSTH, POLD, POLD1 )
      ENDDO

      RETURN
      END
C*********************************************************************
