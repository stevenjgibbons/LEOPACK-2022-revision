C********************************************************************
C subroutine Poloidal / Toroidal coefficient DiSPlay ****************
C            -          -                    - --    ****************
C Steve Gibbons 19.7.99                                             C
C___________________________________________________________________C
C Reads in QST array and displays the poloidal / toroidal values.   C
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
C     RAD       : Value of radius.                                  C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE PTDISP ( LH, LU, QST, RAD )
      IMPLICIT NONE
C___________________________________________________________________C
C Variable declarations - Parameters ...............................C
      INTEGER LH, LU
      DOUBLE PRECISION RAD, QST ( LH*(LH + 2) , 3)
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      INTEGER NH, L, M, ICS, IH
      DOUBLE PRECISION RLL1, SQRLL1, FAC, TOL
      PARAMETER ( TOL = 1.0d-10 )
      CHARACTER *(3) TYPE
C___________________________________________________________________C
C START OF PROGRAM *************************************************C
C___________________________________________________________________C
C
      NH = LH * ( LH + 2 )
C
C First do poloidal harmonics
C
      DO IH = 1, NH
        CALL LMFIND ( IH, L, M, ICS )
        RLL1 = DBLE( L*L + L )
        IF ( ABS( QST( IH, 1) ).GT.TOL ) THEN
          FAC = QST( IH, 1) * RAD / RLL1
          IF ( ICS.EQ.1 ) TYPE = 'COS'
          IF ( ICS.EQ.2 ) TYPE = 'SIN'
          WRITE ( LU  , 789 ) L, M, TYPE, FAC
        ENDIF
 789    FORMAT(' Pol_',I2,'^',I2,' ',A3,' = ',D14.5 )
      ENDDO
C
C Now do toroidal harmonics
C
      DO IH = 1, NH
        CALL LMFIND ( IH, L, M, ICS )
        RLL1 = DBLE( L*L + L )
        SQRLL1 = DSQRT( RLL1 )
        IF ( ABS( QST( IH, 3) ).GT.TOL ) THEN
          FAC = ( -1.0d0 )*QST( IH, 3) / SQRLL1
          IF ( ICS.EQ.1 ) TYPE = 'COS'
          IF ( ICS.EQ.2 ) TYPE = 'SIN'
          WRITE ( LU  , 790 ) L, M, TYPE, FAC
        ENDIF
 790    FORMAT(' Tor_',I2,'^',I2,' ',A3,' = ',D14.5 )
      ENDDO
C
      RETURN
      END
C********************************************************************
