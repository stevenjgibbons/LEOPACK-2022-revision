C*********************************************************************
C subroutine COEfficient File ReaD ***********************************
C            ---         -    -  - ***********************************
C Steve Gibbons 23.6.99                                              C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NDIMMX    : Maximum dimension of the solution vector.          C
C     NRMAX     : Maximum number of radial grid nodes.               C
C     NH        : Number of harmonics (all types)                    C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     FNAME     : File name (unspecified length)                     C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C  Integer                                                           C
C  -------                                                           C
C     NDIM      : Number of elements in SV vector.                   C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SV        : Length ( * ). Length is NH*NR                      C
C                                                                    C
C                                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE COEFRD ( NDIMMX, NRMAX, NH, FNAME, NDIM, NR, SV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIMMX, NRMAX, NH, NDIM, NR
      DOUBLE PRECISION SV( * )
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working Variables .........................C
      INTEGER LU, IWR, IE, NH1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      LU = 999
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
      READ ( LU, * ) NR, NH1
      IF ( NH.NE.NH1 ) THEN
         PRINT *,' Subroutine COEFRD. NH = ', NH
         PRINT *,' There are ',NH1,' harmonics '
         PRINT *,' in the current file. '
         STOP
      ENDIF
      IF ( NR.GT.NRMAX ) THEN
         PRINT *,' Subroutine COEFRD. NR = ', NR
         PRINT *,' NRMAX = ', NRMAX
         STOP
      ENDIF
      NDIM = NH*NR
      IF ( NDIM.GT.NDIMMX ) THEN
         PRINT *,' Subroutine COEFRD. NDIM = ', NDIM
         PRINT *,' NDIMMX = ', NDIMMX
         STOP
      ENDIF
C
      READ ( LU, * ) ( SV( IE ), IE = 1, NDIM )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'error' )
C
      RETURN
      END
C*********************************************************************

