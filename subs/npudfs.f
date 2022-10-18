C*********************************************************************
C subroutine New Program Unit Dependencies Find Subroutine ***********
C            -   -       -    -            -    -          ***********
C Steve Gibbons Sun Mar 19 11:25:49 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Takes a filename, FNAME, and a free logical unit number, LU,       C
C and opens the file (probably with a .dep extension) reading in the C
C dependencies (each an 11 charater long string). If one is found    C
C which is not already in the list, PUDL, the string is added to the C
C list, checking that the length of the list has not been exceeded.  C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LU        : Logical file unit number.                          C
C     MNPU      : Maximum number of program units. (Defines PUDL)    C
C     NIPU      : Number of initial program units. This is the       C
C                 number of filled elements in PUDL when the sub. is C
C                 first called.                                      C
C     NFPU      : Number of final program units.                     C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C     PUDL      : *(11). Dimension (MNPU).                           C
C                 Program unit dependency list.                      C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NPUDFS( LU, MNPU, NIPU, NFPU, FNAME, PUDL )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LU, MNPU, NIPU, NFPU
      CHARACTER *(*) FNAME
      CHARACTER *(11) PUDL( MNPU )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IREAD, IPU
      PARAMETER ( IREAD = 1 )
      CHARACTER *(11) ENTRY
      CHARACTER *(80) LINE
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NFPU = NIPU
      OK = .TRUE.
      CALL FOPEN( LU, FNAME, IREAD )
 80   FORMAT(A)
 49   CONTINUE
      READ ( LU, 80, END=50 ) LINE
      IF ( LINE(1:3).EQ.'END' ) GOTO 50
      IF ( LINE(1:1).EQ.'*' ) GOTO 49
      IF ( LINE(1:5).EQ.'UNFIN' ) THEN
        PRINT *,' Subroutine NPUDFS.'
        PRINT *,' File ', FNAME
        PRINT *,' is not yet finished.'
        PRINT *,' Please complete and re-run.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      ENTRY = LINE(1:11)
C
      DO IPU = 1, NFPU
        IF ( ENTRY.EQ.PUDL(IPU) ) GOTO 49
      ENDDO
      CALL CNTRIC( NFPU, MNPU, OK )
      IF ( OK ) PUDL( NFPU ) = ENTRY
      IF ( .NOT. OK ) THEN
        PRINT *,' Subroutine NPUDFS.'
        PRINT *,' MNPU exceeded.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      GOTO 49
 50   CONTINUE
      CALL FCLOSE( LU, FNAME, 'Error' )
C
      RETURN
      END
C*********************************************************************
