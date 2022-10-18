C*********************************************************************
C subroutine Old Format Integers File Read ***************************
C            -   -      -        -    -    ***************************
C Steve Gibbons Wed Sep 20 08:31:34 BST 200                          C
C____________________________________________________________________C
C                                                                    C
C Reads in an old Leeds-style .inters file and converts to the       C
C new style MHT, MHL, MHM format filling the appropriate arrays.     C
C                                                                    C
C The old format file contains just NH lines of data with            C
C each line containing 4 integer numbers in (I2,I2,I2,I2) format.    C
C                                                                    C
C ttllmmcc                                                           C
C                                                                    C
C tt = 1 --> MHT( ih ) = 4  (Poloidal magnetic field)                C
C tt = 2 --> MHT( ih ) = 5  (Toroidal magnetic field)                C
C                                                                    C
C ll = spherical harmonic degree, l. = MHL( ih )                     C
C                                                                    C
C cc = 1 --> MHM( ih ) = -mm (sin m phi dependence)                  C
C cc = 2 --> MHM( ih ) =  mm (cos m phi dependence)                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Number of vector spherical harmonics. (Output)     C
C     NHMAX     : Maximum number of vector spherical harmonics.      C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NHMAX          C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NHMAX          C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NHMAX          C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE OFIFR( NH, NHMAX, MHT, MHL, MHM, LU, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, NHMAX, MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ), LU
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITT, ILL, IMM, ICC, IWR
      LOGICAL OK
      PARAMETER ( IWR = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      NH = 0
      OK = .TRUE.
      CALL FOPEN ( LU, FNAME, IWR )
 50   CONTINUE
      READ ( LU, 44, END = 60 ) ITT, ILL, IMM, ICC
      CALL CNTRIC( NH, NHMAX, OK )
      IF ( OK ) THEN
C       .
        IF ( ITT.NE.1 .AND. ITT.NE.2 ) THEN
          PRINT *,' Subroutine OFIFR: ITT = ',ITT
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
        IF ( ICC.NE.1 .AND. ICC.NE.2 ) THEN
          PRINT *,' Subroutine OFIFR: ICC = ',ICC
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
        IF ( ITT.EQ.1 ) MHT( NH ) = 4
        IF ( ITT.EQ.2 ) MHT( NH ) = 5
        IF ( ICC.EQ.2 ) MHM( NH ) = IMM
        IF ( ICC.EQ.1 ) MHM( NH ) = -IMM
        MHL( NH ) = ILL
C       .
      ENDIF
C
C Close file
C
      GOTO 50
 60   CONTINUE
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
      IF ( .NOT. OK ) THEN
        PRINT *,' Subroutine OFIFR.'
        PRINT *,' NH = ', NH,' NHMAX = ', NHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
 44   FORMAT(I2,I2,I2,I2)
      RETURN
      END
C*********************************************************************
