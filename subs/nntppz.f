C*********************************************************************
C subroutine Necessary Number of Theta and Phi Points Find ***********
C            -         -         -         -   -      -    ***********
C                                             aZimuthal symmetry.    C
C                                              -                     C
C Steve Gibbons Sat Jul 22 17:12:33 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C For a given level of harmonics LH; this routine will return a      C
C number of theta points ( NTHP ) which is greater then LH           C
C and less than or equal to NTHMAX. Also a number of PHI points      C
C ( NPHPR ) which is greater than 2*MMAX and is also a               C
C power of 2 and is also smaller than NPHPRM. Failiure to do either  C
C will be reported.                                                  C
C                                                                    C
C M0 is the lowest non-zero wavenumber in the solution.              C
C If M0 is a power of 2, then NPHPV is set to NPHPR/M0 to allow      C
C a treatment of azimuthal symmetry in the transform routines.       C
C NPHPV is subject to a limit NPHPVM.                                C
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
C     NPHPR     : Total number of phi points needed to perform       C
C                 Fourier transform.                                 C
C     NPHPV     : Number of phi points required to fully recreate    C
C                 the behaviour of the solution in a subspace of the C
C                 total physical space.                              C
C     M0        : Smallest non-zero wavenumber in solution.          C
C                 We check that NPHPV*M0 = NPHPR.                    C
C     NTHMAX    : Maximum number of theta points.                    C
C     NPHPRM    : Maximum value of NPHPR.                            C
C     NPHPVM    : Maximum value of NPHPV.                            C
C                                                                    C
C     IA        : Aliassing flag.                                    C
C                 If nthp .lt. (1.5*LH) then Aliassing can occur.    C
C                 IA = 1 --> nthp is minimum value.                  C
C                 IA = 2 --> nthp is 1.5*lh, to avoid aliassing.     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NNTPPZ( LH, MMAX, NTHP, NPHPR, NPHPV, M0, NTHMAX,
     1                   NPHPRM, NPHPVM, IA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, MMAX, NTHP, NPHPR, NPHPV, M0, NTHMAX, NPHPRM,
     1        NPHPVM, IA
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IA.NE.1 .AND. IA.NE.2 ) THEN
        PRINT *,' Subroutine NNTPPZ. IA = ', IA
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( M0.LE.0 ) THEN
        PRINT *,' Subroutine NNTPPZ. M0 = ', M0
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LH.EQ.0 ) THEN
        PRINT *,' Subroutine NNTPPZ. LH = 0.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( MMAX.LT.0 .OR. MMAX.GT.LH ) THEN
        PRINT *,' Subroutine NNTPPZ. MMAX = ', MMAX
        PRINT *,' Must be between 0 and LH (= ',LH,')'
        STOP
      ENDIF
C
      NPHPR = 2
 500  CONTINUE
      IF ( 2*MMAX.GE.NPHPR ) THEN
         NPHPR = NPHPR*2
         GOTO 500
      ENDIF
      IF ( NPHPR.GT.NPHPRM ) THEN
         PRINT *,' Subroutine NNTPPZ.'
         PRINT *,' NPHPR must be atleast ', NPHPR
         PRINT *,' NPHPRM = ', NPHPRM
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      NPHPV = NPHPR/M0
C
      IF ( M0*NPHPV.NE.NPHPR ) THEN
        NPHPV = NPHPR
        M0    = 1
      ENDIF
C
      IF ( NPHPV.GT.NPHPVM ) THEN
         PRINT *,' Subroutine NNTPPZ.'
         PRINT *,' NPHPV must be atleast ', NPHPV
         PRINT *,' NPHPVM = ', NPHPVM
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
         PRINT *,' Subroutine NNTPPZ.'
         PRINT *,' NTHP must be atleast ', NTHP
         PRINT *,' NTHMAX = ', NTHMAX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
