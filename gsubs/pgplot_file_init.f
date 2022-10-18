C*********************************************************************
C subroutine PGPLOT FILE INITialisation ******************************
C            ------ ---- ----           ******************************
C Steve Gibbons Mon Jan 22 15:04:51 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Calls PGOPEN, but looks after some of the bookwork.                C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ISTAT     : Integer returned by function PGOPEN.               C
C     IDEV      : Key for selection of device.                       C
C                                                                    C
C            Possibilities are:-                                     C
C                                                                    C
C   idev = 1: gif file. (Landscape)                                  C
C   idev = 2: gif file. (Portrait)                                   C
C   idev = 3: ps file. (Landscape)                                   C
C   idev = 4: ps file. (Portrait)                                    C
C   idev = 5: colour ps file. (Landscape)                            C
C   idev = 6: colour ps file. (Portrait)                             C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     STEM      : Root (maybe including directory tree ) *(*)        C
C                 A space character, " ", terminates STEM.           C
C     CHDEV     : "G" if we have selected gif and "P" for postscript C
C                 (single character string).                         C
C                                                                    C
C*********************************************************************
      SUBROUTINE PGPLOT_FILE_INIT( ISTAT, IDEV, STEM, CHDEV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER        ISTAT, IDEV
      CHARACTER *(*) STEM
      CHARACTER *(1) CHDEV
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          PGOPEN, ILEN
      CHARACTER *(200) DEVICE
      EXTERNAL         PGOPEN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IDEV.LT.1 .OR. IDEV.GT.6 ) THEN
        PRINT *,' Subroutine PGPLOT_FILE_INIT.'
        PRINT *,' IDEV = ', IDEV
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Find length of STEM
C
      ILEN = 0
 50   CONTINUE
      ILEN = ILEN + 1
      IF ( ILEN.GT.190 ) THEN
        PRINT *,' Subroutine PGPLOT_FILE_INIT.'
        PRINT *,' No termination character found in STEM.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      IF ( STEM(ILEN:ILEN).EQ.' ' ) GOTO 60
      GOTO 50
 60   CONTINUE
      ILEN = ILEN - 1
C
C Fix the character string DEVICE.
C
      DEVICE(  1:  1) = '"'
      DEVICE(  2:ILEN+1) = STEM(1:ILEN)
C
      IF ( IDEV.EQ.1 ) THEN
        DEVICE(ILEN+2:ILEN+10) = '.gif"/GIF'
        DEVICE                 = DEVICE(1:ILEN+10)
        CHDEV                  = 'G'
      ENDIF
C
      IF ( IDEV.EQ.2 ) THEN
        DEVICE(ILEN+2:ILEN+11) = '.gif"/VGIF'
        DEVICE                 = DEVICE(1:ILEN+11)
        CHDEV                  = 'G'
      ENDIF
C
      IF ( IDEV.EQ.3 ) THEN
        DEVICE(ILEN+2:ILEN+8)  = '.ps"/PS'
        DEVICE                 = DEVICE(1:ILEN+8)
        CHDEV                  = 'P'
      ENDIF
C
      IF ( IDEV.EQ.4 ) THEN
        DEVICE(ILEN+2:ILEN+9)  = '.ps"/VPS'
        DEVICE                 = DEVICE(1:ILEN+9)
        CHDEV                  = 'P'
      ENDIF
C
      IF ( IDEV.EQ.5 ) THEN
        DEVICE(ILEN+2:ILEN+9)  = '.ps"/CPS'
        DEVICE                 = DEVICE(1:ILEN+9)
        CHDEV                  = 'P'
      ENDIF
C
      IF ( IDEV.EQ.6 ) THEN
        DEVICE(ILEN+2:ILEN+10) = '.ps"/VCPS'
        DEVICE                 = DEVICE(1:ILEN+10)
        CHDEV                  = 'P'
      ENDIF
C
      ISTAT = PGOPEN( DEVICE )
C
      IF ( ISTAT.LE.0 ) THEN
        PRINT *,' Subroutine PGPLOT_FILE_INIT.'
        PRINT *,' PGOPEN has given ISTAT = ', ISTAT
        PRINT *,' PGOPEN was supplied with DEVICE = '
        PRINT *,  DEVICE
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
