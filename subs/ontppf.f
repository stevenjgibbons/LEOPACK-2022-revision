C*********************************************************************
C subroutine Optimum Number of Theta and Phi Points Find *************
C            -       -         -         -   -      -    *************
C Steve Gibbons 9.9.97                                               C
C____________________________________________________________________C
C For a given level of harmonics LH; this routine will return a      C
C number of theta points ( NTHPTS ) which is greater then LH         C
C and less than or equal to NTHMAX. Also a number of PHI points      C
C ( NPHPTS ) which is greater than 2*MMAX and is also a              C
C power of 2 and is also smaller than NPHMAX. Failiure to do either  C
C will be reported.                                                  C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Highest degree, l, of spherical harmonic.          C
C     MMAX      : Highest order, m, of spherical harmonic.           C
C                                                                    C
C     For a fully 3-D problem, MMAX will equal LH but may be         C
C     less if certain symmetries are applied in the azimuthal        C
C     direction.                                                     C
C                                                                    C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     NTHMAX    : Maximum number of theta points.                    C
C     NPHMAX    : Maximum number of phi points.                      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ONTPPF ( LH, MMAX, NTHPTS, NPHPTS, NTHMAX, NPHMAX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, MMAX, NTHPTS, NPHPTS, NTHMAX, NPHMAX
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( LH.EQ.0 ) THEN
        PRINT *,' Subroutine ONTPPF. LH = 0.'
        STOP
      ENDIF
C
      IF ( MMAX.LT.0 .OR. MMAX.GT.LH ) THEN
        PRINT *,' Subroutine ONTPPF. MMAX = ', MMAX
        PRINT *,' Must be between 0 and LH (= ',LH,')'
        STOP
      ENDIF
C
      NPHPTS = 2
 500  CONTINUE
      IF ( 2*MMAX.GE.NPHPTS ) THEN
         NPHPTS = NPHPTS*2
         GOTO 500
      ENDIF
      IF ( NPHPTS.GT.NPHMAX ) THEN
         PRINT *,' Subroutine ONTPPF.'
         PRINT *,' NPHPTS must be atleast ', NPHPTS
         PRINT *,' NPHMAX = ', NPHMAX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      NTHPTS = LH+1
      IF ( NTHPTS/2*2.NE.NTHPTS ) NTHPTS = NTHPTS + 1
      IF ( NTHPTS.GT.NTHMAX ) THEN
         PRINT *,' Subroutine ONTPPF.'
         PRINT *,' NTHPTS must be atleast ', LH+1
         PRINT *,' NTHMAX = ', NTHMAX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
