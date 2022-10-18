C*********************************************************************
C subroutine LATitude in DEGrees 2 THE (in radians) ******************
C            ---         ---     - ---              ******************
C Steve Gibbons Mon Mar 12 11:26:32 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     LATDEG    : Latitude in degrees (between 90 and -90)           C
C     THE       : Co-latitude in radians                             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LATDEG_2_THE( LATDEG, THE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      REAL    LATDEG, THE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      REAL    PI
      PARAMETER ( PI=3.1415926 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( LATDEG.GT.90.0 .OR. LATDEG.LT.-90.0 ) THEN
        PRINT *,' Subroutine LATDEG_2_THE.'
        PRINT *,' LATDEG = ', LATDEG
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      THE = PI*(90.0 - LATDEG)/180.0
C     .
      RETURN
      END
C*********************************************************************

