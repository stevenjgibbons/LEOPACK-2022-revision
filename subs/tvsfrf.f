C*********************************************************************
C function Toroidal Velocity Stress Free Radial Function *************
C          -        -        -      -    -      -        *************
C Steve Gibbons Fri Oct  8 16:57:14 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Fills a function f( r ) = \sum_{m = 0}^{MAXM}                      C
C                 r cos( m pi (r - ri)/(ro - ri) )                   C
C                                                                    C
C This satisfies the boundary conditions                             C
C                                                                    C
C       d [ f ]                                                      C
C      -- [---] = 0     at    r = ri    and    r = ro                C
C      dr [ r ]                                                      C
C                                                                    C
C as required by the toroidal velocity condition with stress free    C
C boundary conditions.                                               C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     R         : Radius                                             C
C     RI        : Radius at inner boundary                           C
C     RO        : Radius at outer boundary                           C
C     COEFS     : c_m in COEFS( m + 1 )                              C
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
      FUNCTION TVSFRF( R, RI, RO, MAXM, COEFS, INTPA, DPRPA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION TVSFRF, R, RI, RO, COEFS( * ), DPRPA( * )
      INTEGER MAXM, INTPA( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER M
      DOUBLE PRECISION PI, FAC, ERR, F2
      PARAMETER (PI=3.14159265358979312D0, ERR = 1.0d-9)
      EXTERNAL CHORCH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( ABS(RO-RI).LT.ERR ) THEN
        PRINT *,' Function TVSFRF.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      IF ( R.LT.RI .OR. R.GT.RO ) THEN
        PRINT *,' Function TVSFRF. R = ', R
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      FAC = PI*(R - RI)/(RO-RI)
C     .
      TVSFRF = 0.0d0
C     .
C     . Loop around values of M
C     .
      DO M = 0, MAXM
C       .
        F2 = FAC*DBLE( M )
        TVSFRF = TVSFRF + R*DCOS( F2 )*COEFS( M + 1)
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
