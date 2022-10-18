C*********************************************************************
C subroutine Zonal Velocity Defining Function ************************
C            -     -        -        -        ************************
C Steve Gibbons Tue Dec  7 12:25:20 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C The zonal velocity for jupiter is defined by                       C
C                                                                    C
C   v_r         = 0                                                  C
C   v_{\theta}  = 0                                                  C
C   v_{\phi}    = s^{nn1} {cos/sin} ( nn2 * \pi ( s - r_i ) )        C
C                                                                    C
C  ZVDF returns the value of v_{\phi} multiplied by the              C
C  normalising factor DPNF. Program aborts if DPNF is zero as        C
C  this is clearly pointless.                                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INTARR    : Integer parameter array. Dimension ( * )           C
C                                                                    C
C                   INTARR( 1 ) = NN1 in above equation.             C
C                   INTARR( 2 ) = NN2 in above equation.             C
C                   INTARR( 3 ) = 1 for COS version                  C
C                                 2 for SIN version                  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Distance from centre of sphere.                    C
C     THE       : Co-latitude in radians.                            C
C                                                                    C
C     DPRARR    : Double precision parameter array. Dimension ( * )  C
C                                                                    C
C                   DPRARR( 1 ) = RI (inner radius value)            C
C                   DPRARR( 2 ) = DPNF (normalising factor).         C
C                                                                    C
C                                                                    C
C                                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION ZVDF ( RAD, THE, DPRARR, INTARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION ZVDF, RAD, THE, DPRARR( * )
      INTEGER INTARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NN1, NN2, ICS
      DOUBLE PRECISION RI, DPNF, LOW, PI, FAC, S
      PARAMETER ( LOW = 1.0d-8, PI = 3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      RI     = DPRARR( 1 )
      DPNF   = DPRARR( 2 )
C
      NN1    = INTARR( 1 )
      NN2    = INTARR( 2 )
      ICS    = INTARR( 3 )
C
      IF ( ICS.NE.1 .AND. ICS.NE.2 ) THEN
        PRINT *,' Function ZVDF.'
        PRINT *,' ICS = ', ICS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DABS( DPNF ).LT.LOW ) THEN
        PRINT *,' Function ZVDF.'
        PRINT *,' DPNF = ', DPNF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( RAD.LT.RI ) THEN
        PRINT *,' Function ZVDF.'
        PRINT *,' RAD = ', RAD,' RI = ', RI
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      S   = RAD * DSIN( THE )
      FAC = DBLE( NN2 ) * PI * ( S - RI )
      IF ( ICS.EQ.1 ) ZVDF = DPNF * S**NN1 * DCOS( FAC )
      IF ( ICS.EQ.2 ) ZVDF = DPNF * S**NN1 * DSIN( FAC )
C
      RETURN
      END
C*********************************************************************
