C*********************************************************************
C double pr. func. Jupiter Velocity Kinetic Energy Integrand Func. ***
C                  -       -        -       -      -         -     ***
C Steve Gibbons Tue Jun 27 10:07:43 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C The zonal velocity for jupiter is defined by                       C
C                                                                    C
C   v_r         = 0                                                  C
C   v_{\theta}  = 0                                                  C
C   v_{\phi}    = dpnf s^{nn1} {cos/sin} ( nn2 * \pi ( s - r_i ) )   C
C                                                                    C
C If we want to calculate the integral in the spherical shell of     C
C the kinetic energy of this flow then we need to integrate, over s, C
C from s = 0 to s = RO, the integrand:                               C
C                                                                    C
C  2 \pi v_{\phi}^2 h( s ) s                                         C
C                                                                    C
C where h( s ) = sqrt( r_o^2 - s^2 ) for s.gt.r_i                    C
C and                                                                C
C       h( s ) = sqrt( r_o^2 - s^2 ) - sqrt( r_i^2 - s^2 )           C
C                                       for s.le.r_i                 C
C                                                                    C
C The factor 2 \pi is for the azimuthal average.                     C
C                                                                    C
C dpnf is a normalisation factor calculated by ZVNFF.                C
C                                                                    C
C We integrate only one hemisphere and then do not divide the k.e.   C
C by 2. The function is designed to be used with the subroutine      C
C NUMINT with input parameters X1 = 0 and X2 = RO                    C
C                                                                    C
C The integer arrays INTPAR and DPPAR must contain the following:-   C
C                                                                    C
C INTPAR( 1 ) = NN1 (in above equation)                              C
C INTPAR( 2 ) = NN2 (in above equation)                              C
C INTPAR( 3 ) = 1 for COS version, 2 for SIN version.                C
C                                                                    C
C DPPAR( 1 ) = RI                                                    C
C DPPAR( 2 ) = DPNF (normalisation factor)                           C
C DPPAR( 3 ) = RO                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INTPAR    : Integer parameter array. Dimension ( * )           C
C                  See above.                                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     S         : Distance from rotation axis.                       C
C     DPPAR     : D.p. parameter array. Dimension ( * )              C
C                  See above.                                        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION JVKEIF( S, INTPAR, DPPAR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INTPAR( * )
      DOUBLE PRECISION S, DPPAR( * ), JVKEIF
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NN1, NN2, ICS
      DOUBLE PRECISION RI, RO, DPNF, VPHI, PI, ZERO, FAC, H
      PARAMETER ( ZERO = 0.0d0, PI = 3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NN1    = INTPAR( 1 )
      NN2    = INTPAR( 2 )
      ICS    = INTPAR( 3 )
C
      IF ( ICS.NE.1 .AND. ICS.NE.2 ) THEN
        PRINT *,' Function JVKEIF.'
        PRINT *,' ICS = ', ICS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RI   = DPPAR( 1 )
      RO   = DPPAR( 3 )
      DPNF = DPPAR( 2 )
C
      IF ( S.LT.ZERO .OR. S.GT.RO ) THEN
        PRINT *,' Function JVKEIF.'
        PRINT *,' S = ', S,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      FAC = DBLE( NN2 ) * PI * ( S - RI )
      IF ( ICS.EQ.1 ) VPHI = DPNF * S**NN1 * DCOS( FAC )
      IF ( ICS.EQ.2 ) VPHI = DPNF * S**NN1 * DSIN( FAC )
C
      H = DSQRT( RO*RO - S*S )
      IF ( S.LT.RI ) THEN
        H = H - DSQRT( RI*RI - S*S )
      ENDIF
C
      JVKEIF = 2.0d0*PI*VPHI*VPHI*H*S
C
      RETURN
      END
C*********************************************************************
