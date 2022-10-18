C*********************************************************************
C integer function Sixteen Bit Character String 2 Number *************
C                  -       -   -         -      - -      *************
C Steve Gibbons Fri Jul 28 09:45:49 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Takes a character string of length sixteen characters, each being  C
C either '0' or '1' and returns the integer value represented.       C
C                                                                    C
C e.g. '0000000000000000' = 0                                        C
C      '0000000000000010' = 2                                        C
C      '1111100011100010' = 63714                                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C   Character *(16)                                                  C
C   ---------------                                                  C
C     CHARST    : Representation of 16 bit number.                   C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION SBCS2N( CHARST )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER         SBCS2N
      CHARACTER *(16) CHARST
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IND, I, IPT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IND    = 17
      SBCS2N = 0
      DO I = 1, 16
        IND = IND - 1
        IF ( CHARST(I:I).NE.'0' .AND. CHARST(I:I).NE.'1' ) THEN
          PRINT *,' Function SBCS2N.'
          PRINT *,' Position ',I,' is not 0 or 1.'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        IF ( I.EQ.1 ) IPT = 1
        IF ( I.GT.1 ) IPT = IPT*2
        IF ( CHARST(IND:IND).EQ.'1' ) SBCS2N = SBCS2N + IPT
      ENDDO
C
      RETURN
      END
C*********************************************************************
