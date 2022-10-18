C*********************************************************************
C subroutine X value ARRay WriTe *************************************
C            -       ---   -  -  *************************************
C Steve Gibbons Fri Nov 12 08:53:38 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Writes out the XARR array of abscissae to a file.                  C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                 Note that NR is not checked for correspondence to  C
C                 any other value - merely for being not greater     C
C                 than NRMAX.                                        C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C     IFORM     : Specifies how the x values are stored on the file. C
C                 Current values are:-                               C
C                                                                    C
C                   IFORM = 1 --> (5(1PD16.7))                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim ( * ) but length atleast NR. Location of       C
C                  radial grid nodes.                                C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XARRWT( NR, XARR, LU, FNAME, IFORM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, LU, IFORM
      CHARACTER *(*) FNAME
      DOUBLE PRECISION XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IWR
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine XARRWT.'
        PRINT *,' IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Open file for writing
C
      IWR = 3
      CALL FOPEN ( LU, FNAME, IWR )
C
C Write number of radial grid nodes
C
      WRITE ( LU, 40 ) NR, IFORM
C
C OK, so write X values ...
C
      IF ( IFORM.EQ.1 ) WRITE ( LU, 41 ) ( XARR( I ), I = 1, NR )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 40   FORMAT(I5,I5)
 41   FORMAT(5(1PD16.7))
      RETURN
      END
C*********************************************************************
