C*********************************************************************
C function Poloidal Magnetic Field (Insulated) Radial Function *******
C          -        -        -      -          -      -        *******
C Steve Gibbons Mon Oct 11 11:59:12 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  If a poloidal radial function is expanded in the form             C
C                                                                    C
C  P( r ) = sum_i c_i cos ( a_i r - b_i )                            C
C                                                                    C
C  when both inner and outer boundaries are insulating,              C
C  then for L, RI and RO (spherical harmonic degree, radius of       C
C  inner and outer boundary respectively) then PMFIRF returns the    C
C  value of P( r ) where c_i is given by COEFS( i ) and the a_i and  C
C  b_i are determined by the subroutine PMFABF.                      C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     R         : Radius                                             C
C     RI        : Radius at inner boundary                           C
C     RO        : Radius at outer boundary                           C
C     COEFS     : c_i in COEFS( i )                                  C
C     DPRPA     : D.p. array( * ). Not referred to                   C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IMAX      : Highest value of I                                 C
C     INTPA     : Int array( * ). INTPA( 1 ) contains sph. harm      C
C                                              degree, L.            C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMFIRF( R, RI, RO, IMAX, COEFS, INTPA, DPRPA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION PMFIRF, R, RI, RO, COEFS( * ), DPRPA( * )
      INTEGER          IMAX, INTPA( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, L, ITMX
      DOUBLE PRECISION COEF, ALPHA, BETA, PI, ERR, BM1, BM2, DIST
      PARAMETER ( PI=3.14159265358979312D0, ERR = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( ABS(RO-RI).LT.ERR ) THEN
        PRINT *,' Function PMFIRF.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( R.LT.RI .OR. R.GT.RO ) THEN
        PRINT *,' Function PMFIRF. R = ', R
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DIST = RO - RI
C
      L      = INTPA( 1 )
      ITMX   = 60
      PMFIRF = 0.0d0
C
C Initialise guesses
C
      ALPHA = PI*0.75d0/DIST
      BETA  = PI*0.5d0
C
C Loop around the i
C
      DO I = 1, IMAX
C
        COEF = COEFS( I )
C
        IF ( I.GT.2 ) BETA = 2.0d0*BM1 - BM2
C
C       Calculate the constants a_i and b_i
C
        CALL PMFABF( ALPHA, BETA, I, RI, RO, L, ITMX )
C
C       Now add onto our total for PMFIRF
C
        PMFIRF = PMFIRF + COEF*COS( ALPHA*R - BETA )
C
C Update guesses
C
        BM2 = BM1
        BM1 = BETA
        ALPHA = ALPHA + PI/DIST
C
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
