C*********************************************************************
C function CHandrasekhar's ORthogonal Cosine Function ****************
C          --              --         -      -        ****************
C Steve Gibbons Thu Oct  7 15:00:53 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C For a value of lambda and x, CHORCF returns the value of           C
C                                                                    C
C          COSH( lambda*x )       COS( lambda*x )                    C
C C( x ) = ----------------   -   ----------------                   C
C          COSH( 0.5 * x )        COS( 0.5 * x )                     C
C                                                                    C
C as defined in Appendix V (page 635) of S. Chandrasekhar,           C
C 'Hydodynamic and Hydromagnetic Stability' 1981, Dover, New York.   C
C                                                                    C
C C( x ) = 0 at x = -0.5 and x = 0.5 and                             C
C dC/dx = 0 at x = -0.5 and x = 0.5 if lambda is the root of         C
C the characteristic equation                                        C
C                                                                    C
C    tanh ( 0.5*lambda ) + tan( 0.5*lambda ) = 0                     C
C                                                                    C
C Function and both input variables are double precision scalars.    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION CHORCF( X, LAMBDA)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION CHORCF, X, LAMBDA
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TOL, Q1, Q2, Q3, Q4, LAMX, HL
      PARAMETER (TOL=1.0d-10)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      LAMX = LAMBDA*X
      HL   = 0.5d0*LAMBDA
      Q1   = COSH( LAMX )
      Q2   = COSH( HL )
      Q3   = DCOS( LAMX )
      Q4   = DCOS( HL )
C
      IF ( ABS( Q2 ).LT.TOL ) THEN
        PRINT *,' Function CHORCF.'
        PRINT *,' COSH(',HL,') = ', Q2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ABS( Q4 ).LT.TOL ) THEN
        PRINT *,' Function CHORCF.'
        PRINT *,' COS(',HL,') = ', Q4
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      CHORCF = Q1/Q2 - Q3/Q4
C
      RETURN
      END
C*********************************************************************
