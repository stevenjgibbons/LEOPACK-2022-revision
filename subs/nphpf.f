C*********************************************************************
C subroutine Number of PHi Points Find *******************************
C            -         --  -      -    *******************************
C Steve Gibbons Mon Nov  6 11:36:13 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Chooses an adequate number of phi points to perform a fast         C
C Fourier transform. Now a function f( phi ) is such that            C
C                                                                    C
C f( phi ) = f( phi + 2.pi/M0 )                                      C
C                                                                    C
C If f requires a maximum wavenumber MMAX then NPHPF finds a number  C
C of points in PHI ( NPHP ) such that  NPHP > 2*MMAX/M0 and NPHP     C
C is a power of 2.                                                   C
C                                                                    C
C The i^{th} point, phi_{i} = (i-1)*deltap                           C
C                                                                    C
C where deltap = 2*pi/(NPHP*M0)                                      C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     MMAX      : Highest wavenumber, m, in spectral expansion.      C
C     M0        : Integer, which divides MMAX such that              C
C                        f( phi ) = f( phi + 2.pi/M0 )               C
C                                                                    C
C     NPHP      : Number of phi points required.                     C
C     NPHMAX    : Maximum number of phi points allowed.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NPHPF( MMAX, M0, NPHP, NPHMAX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER MMAX, M0, NPHP, NPHMAX
C____________________________________________________________________C
C Variable declarations - working variables .........................C
      INTEGER NPHL
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( M0.LE.0 .OR. MMAX.LT.0 ) THEN
        PRINT *,' Subroutine NPHPF.'
        PRINT *,' MMAX = ', MMAX,', M0 = ', M0
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( MMAX/M0*M0.NE.MMAX ) THEN
        PRINT *,' Subroutine NPHPF.'
        PRINT *,' MMAX = ', MMAX,', M0 = ', M0
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NPHL = 2*MMAX/M0
      NPHP = 2
 500  CONTINUE
      IF ( NPHL.GE.NPHP ) THEN
        NPHP = NPHP*2
        GOTO 500
      ENDIF
C
      IF ( NPHP.GT.NPHMAX ) THEN
        PRINT *,' Subroutine NPHPF.'
        PRINT *,' MPHP   = ', NPHP
        PRINT *,' MPHMAX = ', NPHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
