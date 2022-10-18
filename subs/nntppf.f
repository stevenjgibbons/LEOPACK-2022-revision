C*********************************************************************
C subroutine Necessary Number of Theta and Phi Points Find ***********
C            -         -         -         -   -      -    ***********
C Steve Gibbons 9.9.97                                               C
C____________________________________________________________________C
C For a given level of harmonics LH; this routine will return a      C
C number of theta points ( NTHP ) which is greater then LH           C
C and less than or equal to NTHMAX. Also a number of PHI points      C
C ( NPHP ) which is greater than 2*MMAX and is also a                C
C power of 2 and is also smaller than NPHMAX. Failiure to do either  C
C will be reported.                                                  C
C                                                                    C
C If IA = 2, then NTHP is the closest integer to 1.5 * LH.           C
C                                                                    C
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
C     NTHP      : Number of theta points.                            C
C     NPHP      : Number of phi points.                              C
C     NTHMAX    : Maximum number of theta points.                    C
C     NPHMAX    : Maximum number of phi points.                      C
C                                                                    C
C     IA        : Aliassing flag.                                    C
C                 If nthp .lt. (1.5*LH) then Aliassing can occur.    C
C                 IA = 1 --> nthp is minimum value.                  C
C                 IA = 2 --> nthp is 1.5*lh, to avoid aliassing.     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NNTPPF ( LH, MMAX, NTHP, NPHP, NTHMAX, NPHMAX, IA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, MMAX, NTHP, NPHP, NTHMAX, NPHMAX, IA
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IA.NE.1 .AND. IA.NE.2 ) THEN
        PRINT *,' Subroutine NNTPPF. IA = ', IA
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LH.EQ.0 ) THEN
        PRINT *,' Subroutine NNTPPF. LH = 0.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( MMAX.LT.0 .OR. MMAX.GT.LH ) THEN
        PRINT *,' Subroutine NNTPPF. MMAX = ', MMAX
        PRINT *,' Must be between 0 and LH (= ',LH,')'
        STOP
      ENDIF
C
      NPHP = 2
 500  CONTINUE
      IF ( 2*MMAX.GE.NPHP ) THEN
         NPHP = NPHP*2
         GOTO 500
      ENDIF
      IF ( NPHP.GT.NPHMAX ) THEN
         PRINT *,' Subroutine NNTPPF.'
         PRINT *,' NPHP must be atleast ', NPHP
         PRINT *,' NPHMAX = ', NPHMAX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IA.EQ.1 ) THEN
        NTHP = LH+1
      ELSE
        NTHP = LH+LH/2
      ENDIF
C
      IF ( NTHP/2*2.NE.NTHP ) NTHP = NTHP + 1
      IF ( NTHP.GT.NTHMAX ) THEN
         PRINT *,' Subroutine NNTPPF.'
         PRINT *,' NTHP must be atleast ', NTHP
         PRINT *,' NTHMAX = ', NTHMAX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
