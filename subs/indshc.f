C*********************************************************************
C function INDex for Spherical Harmonic Coefficient ******************
C          ---       -         -        -           ******************
C Steve Gibbons 25.4.97                                              C
C____________________________________________________________________C
C Inputs are integers, L and M are obvious ICS = 1 for a cosine harm C
C and ICS = 2 for a sine harmonic.                                   C
C____________________________________________________________________C
      FUNCTION INDSHC ( L, M, ICS )
      IMPLICIT NONE
      INTEGER INDSHC, L, M, ICS
C 
      IF ( L.EQ.0 .AND. M.EQ.0 .AND. ICS.EQ.1 ) THEN
         INDSHC = 0
         RETURN
      ENDIF
C 
      IF ( L.LT.1 ) THEN
        PRINT *,' Function INDSHC. L = ', L
        STOP
      ENDIF
C 
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' Function INDSHC. M invalid. Program Aborted.'
         PRINT *,' L = ', L
         PRINT *,' M = ', M
         STOP
      ENDIF
      IF ( ICS.NE.1 .AND. ICS.NE.2 ) THEN
         PRINT *,' Function INDSHC. ICS = ', ICS
         PRINT *,' ICS must be 1 or 2.'
         PRINT *,' Program Aborted.'
         STOP
      ENDIF
      IF ( M.EQ.0 .AND. ICS.EQ.1 ) THEN
         INDSHC = L*L
         RETURN
      ENDIF
      IF ( M.EQ.0 .AND. ICS.NE.1 ) THEN
         PRINT *,' Function INDSHC. M = ', M
         PRINT *,' ICS = ', ICS,' Program aborted.'
         STOP
      ENDIF
      INDSHC = L*L + 2*M - 2 + ICS
C
      RETURN
      END
C*********************************************************************

