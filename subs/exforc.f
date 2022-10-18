C*********************************************************************
C function EXtended FORmula Coefficients *****************************
C          --       ---     -            *****************************
C Steve Gibbons  3. 8.99                                             C
C____________________________________________________________________C
C                                                                    C
C  This function returns the double precision coefficient for        C
C performing integrals on equally spaced abscissae as described in   C
C Numerical Recipes chapter 4.                                       C
C  Current alternatives are IOPT = 1; extended trapezoidal rule.     C
C                           IOPT = 2; ext. formula to O(1/N^3).      C
C                           IOPT = 3; ext. Simpson's rule.           C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NABS      : Total number of abscissae. Must be atleast 10.     C
C     IABS      : Number of current abscissa.                        C
C     IOPT      : See above. Determines scheme for integration.      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     H         : Distance between two adjacent abscissae.           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION EXFORC ( NABS, IABS, IOPT, H )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NABS, IABS, IOPT
      DOUBLE PRECISION H, EXFORC
C____________________________________________________________________C
C Variable declarations - Working Variables .........................C
      DOUBLE PRECISION FAC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NABS.LT.10 ) THEN
        PRINT *,' Function EXFORC. NABS = ', NABS
        PRINT *,' It should be 10 or greater.'
        STOP
      ENDIF
C
      IF ( IABS.LT.1 .OR. IABS.GT.NABS ) THEN
        PRINT *,' Function EXFORC. IABS = ', IABS
        PRINT *,' NABS = ', NABS,' Sort this out.'
        STOP
      ENDIF
C
C First do the case of extended trapezoidal rule
C
      IF ( IOPT.EQ.1 ) THEN
        IF ( IABS.EQ.1 .OR. IABS.EQ.NABS ) THEN
          EXFORC = 0.5d0*H
        ELSE
          EXFORC = H
        ENDIF
        RETURN
      ENDIF
C
C Now do the case of extended trapezoidal rule O(1/N^3)
C
      IF ( IOPT.EQ.2 ) THEN
        FAC = 1.0d0/12.0d0
        IF ( IABS.EQ.1 .OR. IABS.EQ.NABS ) THEN
          EXFORC = 5.0d0*FAC*H
          RETURN
        ENDIF
        IF ( IABS.EQ.2 .OR. IABS.EQ.(NABS-1) ) THEN
          EXFORC = 13.0d0*FAC*H
        ELSE
          EXFORC = H
        ENDIF
        RETURN
      ENDIF
C
C Now do the case of alternative extended Simpson's rule
C
      IF ( IOPT.EQ.3 ) THEN
        FAC = 1.0d0/48.0d0
        IF ( IABS.EQ.1 .OR. IABS.EQ.NABS ) THEN
          EXFORC = 17.0d0*FAC*H
          RETURN
        ENDIF
        IF ( IABS.EQ.2 .OR. IABS.EQ.(NABS-1) ) THEN
          EXFORC = 59.0d0*FAC*H
          RETURN
        ENDIF
        IF ( IABS.EQ.3 .OR. IABS.EQ.(NABS-2) ) THEN
          EXFORC = 43.0d0*FAC*H
          RETURN
        ENDIF
        IF ( IABS.EQ.4 .OR. IABS.EQ.(NABS-3) ) THEN
          EXFORC = 49.0d0*FAC*H
        ELSE
          EXFORC = H
        ENDIF
        RETURN
      ENDIF
C
      PRINT *,' Function EXFORC. IOPT = ', IOPT
      PRINT *,' IOPT must be set to 1, 2 or 3 '
      STOP
      END
C*********************************************************************

