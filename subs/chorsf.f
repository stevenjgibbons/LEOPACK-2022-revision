C*********************************************************************
C function CHandrasekhar's ORthogonal Sine Function ******************
C          --              --         -    -        ******************
C Steve Gibbons Thu Oct  7 15:00:53 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C For a value of myu and x, CHORSF returns the value of              C
C                                                                    C
C          SINH( myu*x )           SIN( myu*x )                      C
C C( x ) = ----------------   -   ----------------                   C
C          SINH( 0.5 * x )        SIN( 0.5 * x )                     C
C                                                                    C
C as defined in Appendix V (page 635) of S. Chandrasekhar,           C
C 'Hydodynamic and Hydromagnetic Stability' 1981, Dover, New York.   C
C                                                                    C
C C( x ) = 0 at x = -0.5 and x = 0.5 and                             C
C dC/dx = 0 at x = -0.5 and x = 0.5 if myu is the root of            C
C the characteristic equation                                        C
C                                                                    C
C    coth ( 0.5*myu ) - cot( 0.5*myu ) = 0                           C
C                                                                    C
C Function and both input variables are double precision scalars.    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION CHORSF( X, MYU)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION CHORSF, X, MYU
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TOL, Q1, Q2, Q3, Q4, MYUX, HL
      PARAMETER (TOL=1.0d-10)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      MYUX = MYU*X
      HL   = 0.5d0*MYU
      Q1   = SINH( MYUX )
      Q2   = SINH( HL )
      Q3   = DSIN( MYUX )
      Q4   = DSIN( HL )
C
      IF ( ABS( Q2 ).LT.TOL ) THEN
        PRINT *,' Function CHORSF.'
        PRINT *,' SINH(',HL,') = ', Q2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ABS( Q4 ).LT.TOL ) THEN
        PRINT *,' Function CHORSF.'
        PRINT *,' SIN(',HL,') = ', Q4
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      CHORSF = Q1/Q2 - Q3/Q4
C
      RETURN
      END
C*********************************************************************
