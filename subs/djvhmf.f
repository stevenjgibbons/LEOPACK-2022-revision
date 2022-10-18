C*********************************************************************
C subroutine Dudley James Velocity HarMonic Find *********************
C            -      -     -        -        -    *********************
C Steve Gibbons Tue Mar 13 12:26:55 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Enters the indices to be able to store the Dudley James  velocity  C
C harmonics into the arrays MHT, MHL, MHM. For a supplied number     C
C IHARM; DJVHMF will fill the following array elements ...           C
C                                                                    C
C  MHT( iharm    ) = 2, MHL( iharm    ) = 1 and MHM( iharm    ) = 0  C
C  MHT( iharm + 1) = 1, MHL( iharm + 1) = 2 and MHM( iharm + 1) = 0  C
C                                                                    C
C  It is then checked that NH = INARR( 3 ) is atleast IHARM + 1.     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Int. parameter array corresponding to V.           C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR     See INDFUN for details        C
C                 INARR( 3 ) = NH      nrr must equal nr             C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C     MHT( ih ) = 1 if 'ih' is poloidal velocity harmonic.           C
C     MHT( ih ) = 2 if 'ih' is toroidal velocity harmonic.           C
C     MHT( ih ) = 3 if 'ih' is temperature harmonic.                 C
C     MHT( ih ) = 4 if 'ih' is poloidal magnetic field harmonic.     C
C     MHT( ih ) = 5 if 'ih' is toroidal magnetic field harmonic.     C
C                                                                    C
C     MHL       : Spherical harmonic degree, l, of harmonic 'ih'.    C
C                                                                    C
C     MHM       : Spherical harmonic degree, m, of harmonic 'ih'     C
C                if harmonic has (cos m phi dependence) - otherwise  C
C                MHM( ih ) = -m                                      C
C                                                                    C
C     IHARM     : Harmonic at which DJ velocity must begin.          C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DJVHMF( INARR, MHT, MHL, MHM, IHARM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * ), IHARM
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NH  = INARR( 3 )
C     .
      IF ( NH.LT.(IHARM+1) ) THEN
        PRINT *,' Subroutine DJVHMF.'
        PRINT *,' NH = ', NH,' IHARM = ', IHARM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      MHT( IHARM    ) = 2
      MHL( IHARM    ) = 1
      MHM( IHARM    ) = 0
C     .
      MHT( IHARM + 1) = 1
      MHL( IHARM + 1) = 2
      MHM( IHARM + 1) = 0
C     .
      RETURN
      END
C*********************************************************************
