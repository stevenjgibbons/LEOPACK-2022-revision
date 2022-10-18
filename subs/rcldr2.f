C*********************************************************************
C subroutine Real Contour Level Determination Routine 2 **************
C            -    -       -     -             -       - **************
C Steve Gibbons Thu Jan 11 09:10:02 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Given REAL numbers FMIN and FMAX and an odd integer NLEV, RCLDR2   C
C fills the REAL array CONTL with NLEV values such that, if          C
C NLEVB2 = (NLEV-1)/2 and EXTRV = MAX( ABS( FMIN ), ABS( FMAX ) )    C
C then DCONT = EXTRV/REAL( NLEVB2 )  and                             C
C                                                                    C
C        CONTL(  NLEVB2 + 1 - NLEVB2 ) = -EXTRV                      C
C        CONTL(  NLEVB2 + 1 - I      ) = -DCONT*I                    C
C                                                                    C
C        CONTL(  NLEVB2 + 1 - 1      ) = -DCONT                      C
C        CONTL(  NLEVB2 + 1          ) =   0.0                       C
C        CONTL(  NLEVB2 + 1 + 1      ) =  DCONT                      C
C                                                                    C
C        CONTL(  NLEVB2 + 1 + I      ) =  DCONT*I                    C
C        CONTL(  NLEVB2 + 1 + NLEVB2 ) =  EXTRV                      C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NLEV      : Desired number of contour levels.                  C
C                  (must be odd).                                    C
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
      SUBROUTINE RCLDR2( NLEV, FMIN, FMAX, CONTL )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NLEV
      REAL    FMIN, FMAX, CONTL( NLEV )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER     ILEV, NLEVM1, NLEVB2
      REAL        ZERO, EXTRV, DCONT, RVAL
      PARAMETER ( ZERO = 0.0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NLEV.LT.3 .OR. NLEV/2*2.EQ.NLEV ) THEN
        PRINT *,' Subroutine RCLDR2.'
        PRINT *,' NLEV = ', NLEV
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NLEVM1 = NLEV - 1
      NLEVB2 = NLEVM1/2
C
      EXTRV  = ABS( FMAX )
      IF ( ABS( FMIN ).GT.EXTRV ) EXTRV = ABS( FMIN )
C
      IF ( EXTRV.EQ.ZERO ) THEN
        PRINT *,' Subroutine RCLDR2.'
        PRINT *,' EXTRV = ', EXTRV
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DCONT = EXTRV/REAL( NLEVB2 )
C
      CONTL( NLEVB2 + 1 ) = ZERO
      DO ILEV = 1, NLEVB2
        RVAL = REAL( ILEV )*DCONT
        CONTL( NLEVB2 + 1 + ILEV ) = RVAL
        CONTL( NLEVB2 + 1 - ILEV ) = RVAL*(-1.0)
      ENDDO
C
      RETURN
      END
C*********************************************************************
