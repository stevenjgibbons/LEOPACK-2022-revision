C*********************************************************************
C function Poloidal Velocity Stress Free Radial Function *************
C          -        -        -      -    -      -        *************
C Steve Gibbons Fri Oct  8 15:24:51 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Fills a function f( r ) = \sum_{m = 0}^{MAXM}                      C
C                   sin( m pi (r - ri)/(ro - ri) )                   C
C                                                                    C
C This satisfies the boundary conditions                             C
C                                                                    C
C p( ri ) = p( ro ) = p''( ri ) = p''( ro ) = 0                      C
C                                                                    C
C as required by the poloidal velocity condition with stress free    C
C boundary conditions.                                               C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     R         : Radius                                             C
C     RI        : Radius at inner boundary                           C
C     RO        : Radius at outer boundary                           C
C     COEFS     : c_m in COEFS( m )                                  C
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
      FUNCTION PVSFRF( R, RI, RO, MAXM, COEFS, INTPA, DPRPA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION PVSFRF, R, RI, RO, COEFS( * ), DPRPA( * )
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
        PRINT *,' Function PVSFRF.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      IF ( R.LT.RI .OR. R.GT.RO ) THEN
        PRINT *,' Function PVSFRF. R = ', R
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      FAC = PI*(R - RI)/(RO-RI)
C     .
      PVSFRF = 0.0d0
C     .
C     . Loop around values of M
C     .
      DO M = 1, MAXM
C       .
        F2 = FAC*DBLE( M )
        PVSFRF = PVSFRF + DSIN( F2 )*COEFS( M )
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
