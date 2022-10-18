C*********************************************************************
C subroutine Zonal Velocity Normalising Factor Find ******************
C            -     -        -           -      -    ******************
C Steve Gibbons Thu Dec  9 10:51:44 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C The zonal velocity for jupiter is defined by                       C
C                                                                    C
C   v_r         = 0                                                  C
C   v_{\theta}  = 0                                                  C
C   v_{\phi}    = s^{nn1} {cos/sin} ( nn2 * \pi ( s - r_i ) ) dpnf   C
C                                                                    C
C ZVNFF returns in DPNF the scalar required so that v_{\phi}         C
C attains a maximum value of 1 in the spherical shell.               C
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
C     NP        : Number of points to investigate between            C
C                 S = 0.0d0 and S = RO.                              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RI        : Inner radius of spherical shell                    C
C     RO        : Outer radius of spherical shell                    C
C                                                                    C
C     DPNF      : Output. Normalisation factor.                      C 
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ZVNFF( INTARR, NP, RI, RO, DPNF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INTARR( * ), NP
      DOUBLE PRECISION RI, RO, DPNF
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NN1, NN2, ICS, IP
      DOUBLE PRECISION LOW, PI, FAC, S, ZVDF, DELTAS
      PARAMETER ( LOW = 1.0d-8, PI = 3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NN1    = INTARR( 1 )
      NN2    = INTARR( 2 )
      ICS    = INTARR( 3 )
C
      IF ( ICS.NE.1 .AND. ICS.NE.2 ) THEN
        PRINT *,' Subroutine ZVNFF.'
        PRINT *,' ICS = ', ICS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NP.LT.10 ) THEN
        PRINT *,' Subroutine ZVNFF.'
        PRINT *,' NP  = ', NP ,' This is too low for accuracy.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( RO.LT.LOW ) THEN
        PRINT *,' Subroutine ZVNFF.'
        PRINT *,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DELTAS = RO/DBLE( NP - 1 )
      DPNF   = 0.0d0
      DO IP = 1, NP
        S   = DBLE( IP - 1 )*DELTAS
        FAC = DBLE( NN2 ) * PI * ( S - RI )
        IF ( ICS.EQ.1 ) ZVDF = S**NN1 * DCOS( FAC )
        IF ( ICS.EQ.2 ) ZVDF = S**NN1 * DSIN( FAC )
        IF ( DABS( ZVDF ).GT.DABS( DPNF ) ) DPNF = ZVDF
      ENDDO
C
      IF ( DABS( DPNF ).LT.LOW ) THEN
        PRINT *,' Subroutine ZVNFF.'
        PRINT *,' Maximum value = ', DPNF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DPNF = 1.0d0/DPNF
C
      RETURN
      END
C*********************************************************************
