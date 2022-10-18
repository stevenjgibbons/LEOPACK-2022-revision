C********************************************************************
C subroutine QST coefficient DiSPlay ********************************
C            ---             - --    ********************************
C Steve Gibbons  2.8.99                                             C
C___________________________________________________________________C
C Displays non-zero values of the QST array.                        C
C___________________________________________________________________C
C Input Variables :-                                                C
C ===============                                                   C
C  Integers                                                         C
C  --------                                                         C
C     LH        : Maximum degree of spherical harmonics             C
C     LU        : Logical file unit number.                         C
C  Double Precision                                                 C
C  ----------------                                                 C
C     QST       : Array containing the QST coefficients             C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE QSTDSP ( LH, LU, QST )
      IMPLICIT NONE
C___________________________________________________________________C
C Variable declarations - Parameters ...............................C
      INTEGER LH, LU
      DOUBLE PRECISION QST ( LH*(LH + 2) , 3)
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      INTEGER NH, L, M, ICS, IH, IT
      DOUBLE PRECISION TOL
      PARAMETER ( TOL = 1.0d-10 )
      CHARACTER *(3) TYPE
      CHARACTER *(1) AQST
C___________________________________________________________________C
C START OF PROGRAM *************************************************C
C___________________________________________________________________C
C
      NH = LH * ( LH + 2 )
C
C First do scaloidal harmonics
C
      DO IT = 1, 3
       IF ( IT.EQ.1 ) AQST = 'Q'
       IF ( IT.EQ.2 ) AQST = 'S'
       IF ( IT.EQ.3 ) AQST = 'T'
       DO IH = 1, NH
        CALL LMFIND ( IH, L, M, ICS )
        IF ( ABS( QST( IH, IT ) ).GT.TOL ) THEN
          IF ( ICS.EQ.1 ) TYPE = 'COS'
          IF ( ICS.EQ.2 ) TYPE = 'SIN'
          WRITE ( LU  , 789 ) AQST, L, M, TYPE, QST( IH, IT )
        ENDIF
       ENDDO
      ENDDO
 789  FORMAT(' ',A1,'_',I2,'^',I2,' ',A3,' = ',D14.5 )
      RETURN
      END
C********************************************************************
