C*********************************************************************
C subroutine EXTended Formula INtegrate ******************************
C            ---      -       --        ******************************
C Steve Gibbons  3. 8.99                                             C
C____________________________________________________________________C
C                                                                    C
C Integrates a function over a discrete set of values in array VALS. C
C VALS( i ) is the value of a fucntion at X = x_i.                   C
C x_{i+1} - x_i = H.                                                 C
C Formulae according to Numerical Recipes chapter 4.                 C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NABS      : Total number of abscissae. Must be atleast 10.     C
C     IOPT      : Determines scheme for integration.                 C
C                                                                    C
C  Current alternatives are IOPT = 1; extended trapezoidal rule.     C
C                           IOPT = 2; ext. formula to O(1/N^3).      C
C                           IOPT = 3; ext. Simpson's rule.           C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     H         : Distance between two adjacent abscissae.           C
C     VALS      : Array of function values. Dimension ( NABS ).      C
C     APPINT    : Approximation of integral.                         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE EXTFIN ( NABS, IOPT, H, VALS, APPINT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NABS, IOPT
      DOUBLE PRECISION H, VALS( NABS ), APPINT
C____________________________________________________________________C
C Variable declarations - Working Variables .........................C
      DOUBLE PRECISION EXFORC
      INTEGER IABS 
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      APPINT = 0.0d0
      DO IABS = 1, NABS
        APPINT = APPINT +
     1            VALS( IABS ) * EXFORC ( NABS, IABS, IOPT, H )
      ENDDO
C
      RETURN
      END
C*********************************************************************
