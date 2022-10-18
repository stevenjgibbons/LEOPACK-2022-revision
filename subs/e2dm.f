C*********************************************************************
C subroutine E2DM                                                    C
C            ----                                                    C
C Steve Gibbons Mon Apr  9 15:44:51 MET DST 2001                     C
C                                                                    C
C This is a Dave Gubbins routine which I am converting               C
C to double precision and adding safety checks.                      C
C____________________________________________________________________C
C                                                                    C
C  Reads in double precision numbers E0, E1, E2 and E3 along with    C
C a magnetic Reynolds number, DRM. (The latter is optional, only the C
C rescaled magnetic reynolds number DRMS depends upon it).           C
C It returns the values D, M and DRMS (scaled Rm value).             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE E2DM( E0, E1, E2, E3, DRM, D, M, DRMS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION E0, E1, E2, E3, DRM, D, M, DRMS
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION ESQ, DALPHA, DBETA, DGAMMA, DDELTA, DLOW,
     1                 DZERO
      INTEGER          ISIGN0, ISIGN1
C
      PARAMETER ( DALPHA = 16.0D0/315.0D0,
     1            DBETA  = 248832.0D0/26558675.0D0,
     2            DGAMMA = 0.258985312D0,
     3            DDELTA = 0.2686211516D0,
     4            DLOW   = 1.0d-9, DZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      ISIGN0=1
      IF ( E0.LT.DZERO ) ISIGN0 = -1
      ISIGN1=1
      IF ( E1.LT.DZERO ) ISIGN1 = -1
C
      ESQ = (DALPHA*E0*E0) + (DBETA*E1*E1)
     1     + (DGAMMA*E2*E2) + (DDELTA*E3*E3)
C
      IF ( ESQ.LT.DLOW ) THEN
        PRINT *,' Subroutine E2DM'
        PRINT *,' ESQ = ', ESQ
        PRINT *,' Trouble ahead! Program aborted.'
        STOP
      ENDIF
C
      D = (ISIGN0*DALPHA*E0*E0)/ESQ
C
      M = (ISIGN1*DBETA*E1*E1)/ESQ
C
      DRMS = DRM*DSQRT( ESQ )
C
      RETURN
      END
C*********************************************************************
