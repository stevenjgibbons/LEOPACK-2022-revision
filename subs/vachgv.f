C*********************************************************************
C subroutine VAlue of CHebyshev polynomial Gibbons Version ***********
C            --       --                   -       -       ***********
C Steve Gibbons Fri Apr 28 12:38:41 BST 2000                         C
C                                                                    C
C This routine is almost entirely derived from VACHGA written        C
C by D. Funaro of Consiglio Nazionale Delle Ricerche,                C
C Via Abbiategrasso, 209 - 27100 Pavia, Italy.                       C
C                                                                    C
C I have merely changed the IMPLICIT DOUBLE PRECISION statement to   C
C an implicit none statement and have declared all subsequent        C
C variables.                                                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N         : Degree of Chebshev polynomial, T_N.                C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     X         : value of x in interval [-1,1]                      C
C     Y         : value of Cheb. polynomial T_N                      C
C     DY        : 1st derivative of Cheb. polynomial T_N             C
C     D2Y       : 2nd derivative of Cheb. polynomial T_N             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VACHGV( N, X, Y, DY, D2Y )
C**************************************************************
C Following seven lines are Funaro notes.
C   COMPUTES THE VALUE OF THE CHEBYSHEV POLYNOMIAL OF DEGREE N
C   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
C   N  = DEGREE OF THE POLYNOMIAL
C   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
C   Y  = VALUE OF THE POLYNOMIAL IN X
C   DY = VALUE OF THE FIRST DERIVATIVE IN X
C   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
C**************************************************************
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X, Y, DY, D2Y
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER K
      DOUBLE PRECISION YP, DYP, D2YP, YM, DYM, D2YM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      Y   = 1.0D0
      DY  = 0.0D0
      D2Y = 0.0D0
      IF (N .EQ. 0) RETURN
      Y   = X
      DY  = 1.0D0
      D2Y = 0.0D0
      IF (N .EQ. 1) RETURN
      YP   = 1.0D0
      DYP  = 0.0D0
      D2YP = 0.0D0
      DO 1 K = 2, N
         YM  = Y
         Y   = 2.0D0*X*Y-YP
         YP  = YM
         DYM = DY
         DY  = 2.0D0*X*DY+2.0D0*YP-DYP
         DYP = DYM
         D2YM= D2Y
         D2Y = 2.0D0*X*D2Y+4.0D0*DYP-D2YP
         D2YP= D2YM
 1    CONTINUE
C
      RETURN
      END
C*********************************************************************
