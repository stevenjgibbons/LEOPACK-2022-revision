C*********************************************************************
C subroutine Real Contour Level Determination Routine ****************
C            -    -       -     -             -       ****************
C Steve Gibbons Thu Jan 11 09:10:02 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Given REAL numbers FMIN and FMAX and an integer NLEV, RCLDR fills  C
C the REAL array CONTL with NLEV values such that the contours are   C
C evenly spaced. FMIN must be strictly less than FMAX. If FMIN is    C
C negative and FMAX is positive, then it is arranged that one of the C
C contour lines is identical to zero.                                C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NLEV      : Desired number of contour levels.                  C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     FMIN      : Lowest value in the data set.                      C
C     FMAX      : Highest value in the data set.                     C
C     CONTL     : Dim (NLEV). Resulting contour levels.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RCLDR( NLEV, FMIN, FMAX, CONTL )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NLEV
      REAL    FMIN, FMAX, CONTL( NLEV )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER     ILEV, CINTP, CINTN
      REAL        RINT, RINTP, RINTN, FAC, ZERO
      PARAMETER ( ZERO = 0.0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( FMAX.LE.FMIN ) THEN
        PRINT *,' Subroutine RCLDR.'
        PRINT *,' FMIN = ', FMIN
        PRINT *,' FMAX = ', FMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C So we know now that FMAX is strictly greater than FMIN.
C
      IF (     (FMIN.LT.ZERO .AND. FMAX.LE.ZERO) .OR.
     1         (FMIN.GE.ZERO .AND. FMAX.GT.ZERO)       ) THEN
        IF ( NLEV.LT.2 ) THEN
          PRINT *,' Subroutine RCLDR.'
          PRINT *,' FMAX  = ', FMAX
          PRINT *,' FMIN  = ', FMIN
          PRINT *,' NLEV  = ', NLEV
          PRINT *,' You need bigger NLEV!'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        RINT = (FMAX - FMIN)/REAL( NLEV - 1 )
        DO ILEV = 1, NLEV
          CONTL( ILEV ) = FMIN + RINT*( ILEV - 1 )
        ENDDO
      ENDIF
C
      IF ( FMIN.LT.ZERO .AND. FMAX.GT.ZERO ) THEN
        FAC   = REAL( NLEV - 1 )*FMAX/(FMAX - FMIN )
        CINTP = INT( FAC )
        IF ( CINTP.EQ.0 ) CINTP = 1
        CINTN = NLEV - 1 - CINTP
        IF ( CINTN.LE.0 ) THEN
          PRINT *,' Subroutine RCLDR.'
          PRINT *,' FMAX  = ', FMAX
          PRINT *,' FMIN  = ', FMIN
          PRINT *,' CINTN = ', CINTN
          PRINT *,' NLEV  = ', NLEV
          PRINT *,' You need bigger NLEV!'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        RINTP = FMAX/REAL( CINTP )
        RINTN = -FMIN/REAL( CINTN )
        RINT  = MAX( RINTP, RINTN )
        DO ILEV = 1, CINTN
          CONTL( ILEV ) = RINT*(-1.0)*REAL(CINTN+1-ILEV)
        ENDDO
        CONTL( CINTN + 1 ) = ZERO
        DO ILEV = 1, CINTP
          CONTL( CINTN + 1 + ILEV ) = RINT*REAL(ILEV)
        ENDDO
      ENDIF
C
      RETURN
      END
C*********************************************************************
