C*********************************************************************
C subroutine SiMple ITerative SoLve **********************************
C            - -    --        - -   **********************************
C Steve Gibbons Thu Oct  7 16:26:59 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C If the double precision function FUNC, which has been declared     C
C EXTERNAL in the (sub)program which calls SMITSL, gives f( X )      C
C for a given value of X then SMITSL will iterate linearly towards   C
C a solution of f( X ) = 0.                                          C
C                                                                    C
C FUNC must have the argument list ( X, INTARR, DPRARR )             C
C                                                                    C
C with INTARR and DPRARR being respectively integer and double       C
C precision arrays of undefined length. INTARR and DPRARR are not    C
C referred to by SMITSL other than in the call to FUNC.              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FUNC      : Declared EXTERNALly. Has the form                  C
C                 FUNC( X, INTPAR, DPPAR )                           C
C                 (See above).                                       C
C                                                                    C
C     GUESS     : First value of X to be used.                       C
C     FACTOR    : Second value of X to be used to begin the          C
C                 iteration is FACTOR*GUESS                          C
C                                                                    C
C                 So if FACTOR = 1.1, then X2 = 1.1*X1               C
C                                                                    C
C     ERR       : User specified tolerance of error of the solution. C
C                 If f( X ) has a modulus less than ERR then X is    C
C                 returned as SOL.                                   C
C                                                                    C
C     SOL       : Approximation of solution.                         C
C                                                                    C
C     DPRARR    : Double precision array of parameters for FUNC.     C
C                                                  Dim (*)           C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ITMX      : Maximum number of iterations permitted.            C
C                                                                    C
C     INFO      : Output variable. If solution converges             C
C                 succesfully then INFO returns the number of        C
C                 iterations required to converge - this is          C
C                 an integer between 1 and ITMX.                     C
C                                                                    C
C                 A failiure will result in a negative value         C
C                 of INFO.                                           C
C                                                                    C
C                 INFO = -1 means two consecutive X values were      C
C                                                      identical.    C
C                                                                    C
C                 INFO = -2 means two consecutive f(X) values were   C
C                                                      identical.    C
C                                                                    C
C                 INFO = -3 means the maximum number of iterations   C
C                                          were exceeded.            C
C                                                                    C
C     INTPAR    : Integer array of parameters for FUNC. Dim (*)      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SMITSL( FUNC, GUESS, FACTOR, ITMX, INFO, ERR, SOL,
     1                   DPRARR, INTARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION FUNC, GUESS, FACTOR, ERR, SOL, DPRARR( * )
      INTEGER ITMX, INFO, INTARR( * )
      EXTERNAL FUNC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION X, FX, XOLD, FXOLD, DPM, DPC, FXSUB, XSUB,
     1                 TOL
      PARAMETER ( TOL = 1.0d-8 )
      INTEGER NOIT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Perform initial stage
C
      XOLD = GUESS
      FXOLD = FUNC( XOLD, INTARR, DPRARR )
C
C Now set X as a factor of FACTOR * GUESS
C
      X = FACTOR * GUESS
      FX = FUNC( X, INTARR, DPRARR )
C
      NOIT = 0
C
C Now begin iterative loop
C
 50   CONTINUE
      NOIT = NOIT + 1
      XSUB = X - XOLD
C
C Error trapping. If we have two identical consecutive
C values for x, this will result in a division by zero
C Set INFO = -1
C
      IF ( ABS( XSUB ).LT.TOL ) THEN
        INFO = -1
        RETURN
      ENDIF
C
C Set the gradient of our straight line to FXSUB/XSUB
C
      FXSUB = FX - FXOLD
      DPM = FXSUB/XSUB
C
C Error trapping. If gradient is zero then the function
C has had the same value for the previous iterations.
C Now this might be a (very unfortunate!) coincidence,
C or it could mean that the function is constant and
C no root can be found. In any case, this algorithm
C is unable to get out of this fix ...
C Set INFO = -2
C
      IF ( ABS( DPM ).LT.TOL ) THEN
        INFO = -2
        RETURN
      ENDIF
C
C Ok so gradient is non-zero and we can take the y
C intercept, DPC.
C
      DPC = FX - DPM*X
C
C So calculate our linearly extrapolated solution
C
      SOL = (-1.0d0)*DPC/DPM
C
C Set XOLD and FXOLD to X and FX
C
      XOLD = X
      FXOLD = FX
C
C Set X to SOL and calculate FX
C
      X = SOL
      FX = FUNC( X, INTARR, DPRARR )
C
C Check to see if we have converged
C If so then set INFO to the number of iterations taken
C and return
C
      IF ( ABS( FX ).LT.ERR ) THEN
        INFO = NOIT
        RETURN
      ENDIF
C
C See if we have reached the maximum number of iterations ...
C Set INFO = -3 and return
C
      IF ( NOIT.EQ.ITMX ) THEN
        INFO = -3
        RETURN
      ENDIF
C
C No problems, but no solution either, so go back
C to beginning of loop
C
      GOTO 50
C
      END
C*********************************************************************
