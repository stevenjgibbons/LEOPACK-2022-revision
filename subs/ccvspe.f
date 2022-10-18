C*********************************************************************
C subroutine Chebyshev Coefficient Vector Single Point Evaluate ******
C            -         -           -      -      -     -        ******
C Steve Gibbons Sat Apr 29 08:55:27 BST 2000                         C
C                                                                    C
C This routine is almost entirely derived from PVCHEX written        C
C by D. Funaro of Consiglio Nazionale Delle Ricerche,                C
C Via Abbiategrasso, 209 - 27100 Pavia, Italy.                       C
C                                                                    C
C The vector CCV contains NH radial functions, f(r),                 C
C                                                                    C
C If a function f( r ) = \sum_{k=0}^{n-1} c_k T_k( x ),              C
C (where x is the radius, transformed to be in the interval [-1,1])  C
C                                                                    C
C and c_k is contained in CCV( ind ), where ind is given by the      C
C integer function CHINDF, then PVCHGV returns f, df/dr and f''(r)   C
C in Y, DY and D2Y respectively.                                     C
C                                                                    C
C n in the above expansion is given by INARR( 2 ) - 2.               C
C INARR( 2 ) is the number of radial grid nodes that the finite      C
C difference solution was originally stored in.                      C
C                                                                    C
C There are no coefficients c_n and c_{n+1} in the expansion of f(r) C
C and the locations in CCV where these coefficients would be stored  C
C (according to the function CHINDF) hold RI and RO respectively.    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Dim (3). INARR( 1 ) is the format (see CHINDF).    C
C                          INARR( 2 ) = n+2                          C
C                          INARR( 3 ) is the number of radial func.s C
C                                                                    C
C     IH        : Number of radial function to be considered.        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     R         : value of r in interval [RI,RO]                     C
C     CCV       : Dim (*). Chebyshev coefficient vector. See above.  C
C     Y         : value of Cheb. polynomial T_N                      C
C     DY        : 1st derivative of Cheb. polynomial T_N             C
C     D2Y       : 2nd derivative of Cheb. polynomial T_N             C
C                                                                    C
C     the derivatives DY and D2Y are scaled appropriately for the    C
C     functions of r over (ri,ro). Using the Chebyshev recurrence    C
C     relations directly give derivatives with respect to x over     C
C     the interval (-1,1).                                           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CCVSPE( INARR, IH, R, CCV, Y, DY, D2Y )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), IH
      DOUBLE PRECISION R, CCV( * ), Y, DY, D2Y
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER K, N, INDCCV, CHINDF, IDIR
      DOUBLE PRECISION P, DP, D2P, PP, DPP, D2PP, PM, DPM, D2PM,
     1                 RI, RO, X, FAC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      N             = INARR( 2 ) - 2
C
      K             = N
      INDCCV        = CHINDF( K, IH, INARR )
      RI            = CCV( INDCCV )
C
      K             = N+1
      INDCCV        = CHINDF( K, IH, INARR )
      RO            = CCV( INDCCV )
C
      IF ( R.LT.RI .OR. R.GT.RO ) THEN
        PRINT *,' Subroutine CCVSPE.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' R  = ', R
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C------------------------
C FAC is necessary to scale the derivatives
C as they are of functions between RI and RO
C and not -1 and 1.
C------------------------
      FAC = 2.0d0/(RO-RI)
C------------------------
C Now calculate the value of X in the
C interval [-1,1] which corresponds to 
C R in the interval [RI,RO].
C------------------------
      IDIR = 1
      CALL CPCOVR( IDIR, R, X, RI, RO )
C------------------------
C
      K             = 0
      INDCCV        = CHINDF( K, IH, INARR )
      Y   = CCV( INDCCV )
      DY  = 0.0D0
      D2Y = 0.0D0
      P   = X
      DP  = 1.0D0
      D2P = 0.0D0
      K             = 1
      INDCCV        = CHINDF( K, IH, INARR )
      Y   = Y+P*CCV( INDCCV )
      DY  = DP*CCV( INDCCV )
      D2Y = 0.0D0
      PP  = 1.0D0
      DPP = 0.0D0
      D2PP = 0.0D0
      DO 1 K = 2, N-1
          INDCCV = CHINDF( K, IH, INARR )
          PM = P
          P  = 2.0D0*X*P-PP
          Y  = Y+P*CCV( INDCCV )
          PP = PM
          DPM = DP
          DP  = 2.0D0*X*DP+2.0D0*PP-DPP
          DY  = DY+DP*CCV( INDCCV )
          DPP  = DPM
          D2PM = D2P
          D2P  = 2.0D0*X*D2P+4.0D0*DPP-D2PP
          D2Y  = D2Y+D2P*CCV( INDCCV )
          D2PP = D2PM
 1    CONTINUE
C
C Now simply rescale derivatives so they are
C functions with respect to R and not to X.
C
      DY  = DY*FAC
      D2Y = D2Y*FAC*FAC
C
      RETURN
      END
C*********************************************************************
