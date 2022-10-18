C*********************************************************************
C subroutine NEWTon Raphson Solve ************************************
C            ----   -       -     ************************************
C Steve Gibbons Sat Oct  9 15:14:31 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C If the double precision function FUNC, which has been declared     C
C EXTERNAL in the (sub)program which calls NEWTRS, gives either      C
C f( x ) of df/dx( x ) depending on whether the integer NDER = 0 or  C
C NDER = 1 for a given value of X, then NEWTRS will iterate linearly C
C towards a solution of f( X ) = 0.                                  C
C                                                                    C
C FUNC must have the argument list ( X, NDER, INTARR, DPRARR )       C
C                                                                    C
C with INTARR and DPRARR being respectively integer and double       C
C precision arrays of undefined length. INTARR and DPRARR are not    C
C referred to by NEWTRS other than in the call to FUNC.              C
C                                                                    C
C The root must be contained in the interval [ LOW, HIGH ]           C
C                                                                    C
C This routine was guided by the C routine rtsafe in Numerical       C
C Recipes in C.                                                      C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FUNC      : Declared EXTERNALly. Has the form                  C
C                 FUNC( X, NDER, INTPAR, DPPAR )                     C
C                 (See above).                                       C
C                                                                    C
C     DLOW      : Lower bound to interval containing root.           C
C     HIGH      : Upper bound to interval containing root.           C
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
C                 INFO = -1 means no root is located between DLOW    C
C                           and HIGH. There may be many roots in     C
C                           this interval, but Newton-Raphson will   C
C                           not find them easily.                    C
C                                                                    C
C                 INFO = -2 means a zero value of derivative was     C
C                                                   encountered.     C
C                                                                    C
C                 INFO = -3 means the maximum number of iterations   C
C                                          were exceeded.            C
C                                                                    C
C     INTPAR    : Integer array of parameters for FUNC. Dim (*)      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NEWTRS( FUNC, DLOW, HIGH, ITMX, INFO, ERR,
     1                   SOL, DPRARR, INTARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION FUNC, DLOW, HIGH, ERR, SOL, DPRARR( * )
      INTEGER ITMX, INFO, INTARR( * )
      EXTERNAL FUNC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION XL, FL, XH, FH, DX, DXOLD,
     1                 FX, DFX, VAL1, VAL2, TOL, X
      PARAMETER ( TOL = 1.0d-8 )
      INTEGER NOIT, NDER
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First check that there is a root between DLOW and HIGH
C (the following test does not prove there isn't a root
C but it proves that Newton-Raphson is not the way to find it!)
C
C Find function at DLOW and HIGH
C
      NDER = 0
      FL   = FUNC( DLOW, NDER, INTARR, DPRARR )
      FH   = FUNC( HIGH, NDER, INTARR, DPRARR )
C
C Return if FL and FH are of the same sign
C
      VAL1 = FL*FH
      IF ( VAL1.GT.0.0d0 ) THEN
        INFO = -1
        RETURN
      ENDIF
C
C Early exit if we have already found a root ...
C
      IF ( ABS( FL ).LT.ERR ) THEN
        INFO = 0
        SOL = DLOW
        RETURN
      ENDIF
C
      IF ( ABS( FH ).LT.ERR ) THEN
        INFO = 0
        SOL = HIGH
        RETURN
      ENDIF
C
C Orient the search so that f( xl ) < 0
C
      IF ( FL.LT.0.0d0 ) THEN
        XL = DLOW
        XH = HIGH
      ELSE
        XL = HIGH
        XH = DLOW
      ENDIF
C
C Now set X to a guess between DLOW and HIGH
C
      X     = 0.5d0*( DLOW + HIGH )
      DXOLD = ABS( HIGH - DLOW )
      DX    = DXOLD
C
      NOIT = 0
C
C Evaluate the function at this new X value
C
      NDER = 0
      FX = FUNC( X, NDER, INTARR, DPRARR )
      NDER = 1
      DFX = FUNC( X, NDER, INTARR, DPRARR )
C
C Now begin iterative loop
C
 50   CONTINUE
      NOIT = NOIT + 1
C
C Check that maximum number of iterations
C have not been exceeded
C
      IF ( NOIT.GT.ITMX ) THEN
        INFO = -3
        RETURN
      ENDIF
C
      VAL1 = (X-XH)*DFX - FX
      VAL2 = (X-XL)*DFX - FX
C
C Now we bisect if Newton is out of range ...
C
      IF (    (VAL1*VAL2).GT.0.0d0     .OR.
     1        ABS(2.0d0*FX).GT.ABS(DXOLD*DFX) ) THEN
        DXOLD = DX
        DX    = 0.5d0*(XH - XL)
        X     = XL + DX
      ELSE
C
C or Newton step is acceptable ...
C
        DXOLD = DX
        IF ( ABS( DFX ).LT.TOL ) THEN
          INFO = -2
          RETURN
        ENDIF
C       .
        DX = FX/DFX
        X = X - DX
C       .
      ENDIF
C
C Evaluate the function at this new X value
C
      NDER = 0
      FX = FUNC( X, NDER, INTARR, DPRARR )
      NDER = 1
      DFX = FUNC( X, NDER, INTARR, DPRARR )
C
      IF ( ABS( FX ).LT.ERR ) THEN
        INFO = NOIT
        SOL = X
        RETURN
      ENDIF
C
      IF ( FX.LT.0.0d0 ) THEN
        XL = X
      ELSE
        XH = X
      ENDIF
C
      GOTO 50
C
      END
C*********************************************************************
