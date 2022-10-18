C*********************************************************************
C function Poloidal Velocity No Slip Radial Function *****************
C          -        -        -  -    -      -        *****************
C Steve Gibbons Fri Oct  8 09:46:45 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C This function is a series of functions as detailed in              C
C Chandrasekhar Appendix V, page 635                                 C
C 'Hydodynamic and Hydromagnetic Stability' 1981, Dover, New York.   C
C                                                                    C
C PVNSRF( r ) = \sum_{m=1}^{m=MAXM} \left(                           C
C                 c_m C_m( x ) + s_m S_m( x ) \right)                C
C                                                                    C
C where c_m = COEFS( 2m - 1 ), s_m = COEFS( 2m ),                    C
C                                                                    C
C x = ( r - ri )/(ro - ri ) - 0.5d0,                                 C
C                                                                    C
C            COSH( lambda_m*x )       COS( lambda_m*x )              C
C C_m( x ) = ----------------     -    ----------------              C
C             COSH( 0.5 * x )           COS( 0.5 * x )               C
C                                                                    C
C             SINH( myu_m*x )           SIN( myu_m*x )               C
C S_m( x ) = ----------------     -    ----------------              C
C             SINH( 0.5 * x )           SIN( 0.5 * x )               C
C                                                                    C
C  lambda_m and myu_m are respectively the m^{th} non-zero roots of  C
C                                                                    C
C tanh ( 0.5*lambda ) + tan( 0.5*lambda ) = 0 and                    C
C                                                                    C
C coth ( 0.5*myu ) - cot( 0.5*myu ) = 0                              C
C                                                                    C
C which are calculated by calls to SMITSL with different options     C
C on the EXTERNAL function CHORCH.                                   C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     R         : Radius                                             C
C     RI        : Radius at inner boundary                           C
C     RO        : Radius at outer boundary                           C
C     COEFS     : c_m in COEFS( 2m - 1 )                             C
C                 s_m in COEFS( 2m )                                 C
C     DPRPA     : D.p. array( * ). Not referred to                   C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     MMAX      : Highest value of M                                 C
C     INTPA     : Int array( * ). Not referred to                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PVNSRF( R, RI, RO, MAXM, COEFS, INTPA, DPRPA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION PVNSRF, R, RI, RO, COEFS( * ), DPRPA( * )
      INTEGER MAXM, INTPA( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER M, IT, IND, INTARR( 1 ), ITMX, INFO
      DOUBLE PRECISION SOL, ERR, X, CHORCH, CHORCF, CHORSF, CCF, SCF,
     1                 GUESS, FACTOR, PI
      PARAMETER (PI=3.14159265358979312D0, ERR = 1.0d-9, ITMX = 40 )
      EXTERNAL CHORCH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      FACTOR = 1.02d0
      IF ( ABS(RO-RI).LT.ERR ) THEN
        PRINT *,' Function PVNSRF.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      IF ( R.LT.RI .OR. R.GT.RO ) THEN
        PRINT *,' Function PVNSRF. R = ', R
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Convert R (between RI and RO) to X
C     . (between -0.5 and 0.5 )
C     .
      X = (R-RI)/(RO-RI) - 0.5d0
C     .
      PVNSRF = 0.0d0
C     .
C     . Loop around values of M
C     .
      DO M = 1, MAXM
C       .
C       . it = 1  --> C_m( x )  (see CHORCF)
C       .
        IT = 1
        INTARR( 1 ) = IT
        IND = 2*M - 1
        CCF = COEFS( IND )
C
C Try to locate the M^th zero of 
C tanh( 0.5*x ) + tan( 0.5*x )
C initial guess for this is (2m-0.5)*pi
C
        GUESS = (2.0d0*DBLE( M ) - 0.5d0 )*PI
        CALL SMITSL( CHORCH, GUESS, FACTOR, ITMX, INFO, ERR,
     1               SOL, DPRPA, INTARR )
C
C Check that we have successfully located a root
C
        IF ( INFO.LT.0 ) THEN
          PRINT *,' Function PVNSRF.'
          PRINT *,' SMITSL has returned INFO = ', INFO
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
C So SOL now contains lambda_m; so add contribution
C to PVNSRF
C
        PVNSRF = PVNSRF + CCF*CHORCF( X, SOL )
C       .
C       . it = 2  --> S_m( x )  (see CHORSF)
C       .
        IT = 2
        INTARR( 1 ) = IT
        IND = 2*M
        SCF = COEFS( IND )
C       .
C Try to locate the M^th zero of 
C coth( 0.5*x ) - cot( 0.5*x )
C initial guess for this is (2m+0.5)*pi
C
        GUESS = (2.0d0*DBLE( M ) + 0.5d0 )*PI
        CALL SMITSL( CHORCH, GUESS, FACTOR, ITMX, INFO, ERR,
     1               SOL, DPRPA, INTARR )
C
C Check that we have successfully located a root
C
        IF ( INFO.LT.0 ) THEN
          PRINT *,' Function PVNSRF.'
          PRINT *,' SMITSL has returned INFO = ', INFO
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
C So SOL now contains lambda_m; so add contribution
C to PVNSRF
C
        PVNSRF = PVNSRF + SCF*CHORSF( X, SOL )
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
