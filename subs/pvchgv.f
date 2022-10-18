C*********************************************************************
C subroutine Polynomial eValuate CHebyshev Gibbons Version ***********
C            -           -       --        -       -       ***********
C Steve Gibbons Fri Apr 28 12:38:41 BST 2000                         C
C                                                                    C
C This routine is almost entirely derived from PVCHEX written        C
C by D. Funaro of Consiglio Nazionale Delle Ricerche,                C
C Via Abbiategrasso, 209 - 27100 Pavia, Italy.                       C
C                                                                    C
C I have merely changed the IMPLICIT DOUBLE PRECISION statement to   C
C an implicit none statement and have declared all subsequent        C
C variables.                                                         C
C                                                                    C
C If a function q( x ) = \sum_{k=0}^{n} c_k T_k( x ),                C
C and c_k is contained in CO( k ), then PVCHGV returns q, dq/dx and  C
C q''(x) in Y, DY and D2Y respectively.                              C
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
C     CO        : CO( k ) = c_k. Dim (0:N)                           C
C     Y         : value of Cheb. polynomial T_N                      C
C     DY        : 1st derivative of Cheb. polynomial T_N             C
C     D2Y       : 2nd derivative of Cheb. polynomial T_N             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PVCHGV( N, X, CO, Y, DY, D2Y )
C*********************************************************************
C   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND
C   SECOND DERIVATIVES BY KNOWING THE CHEBYSHEV FOURIER COEFFICIENTS
C   N  = THE DEGREE OF THE POLYNOMIAL
C   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
C   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
C   Y  = VALUE OF THE POLYNOMIAL IN X
C   DY = VALUE OF THE FIRST DERIVATIVE IN X
C   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
C*********************************************************************
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X, CO(0:N), Y, DY, D2Y
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER K
      DOUBLE PRECISION P, DP, D2P, PP, DPP, D2PP, PM, DPM, D2PM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      Y   = CO(0)
      DY  = 0.0D0
      D2Y = 0.0D0
      IF (N .EQ. 0) RETURN
      P   = X
      DP  = 1.0D0
      D2P = 0.0D0
      Y   = CO(0)+P*CO(1)
      DY  = DP*CO(1)
      D2Y = 0.0D0
      IF (N .EQ. 1) RETURN
      PP  = 1.0D0
      DPP = 0.0D0
      D2PP = 0.0D0
      DO 1 K = 2, N
          PM = P
          P  = 2.0D0*X*P-PP
          Y  = Y+P*CO(K)
          PP = PM
          DPM = DP
          DP  = 2.0D0*X*DP+2.0D0*PP-DPP
          DY  = DY+DP*CO(K)
          DPP  = DPM
          D2PM = D2P
          D2P  = 2.0D0*X*D2P+4.0D0*DPP-D2PP
          D2Y  = D2Y+D2P*CO(K)
          D2PP = D2PM
 1    CONTINUE
C
      RETURN
      END
C*********************************************************************
