C*********************************************************************
C subroutine Two Dimensional Real Array ReaD *************************
C            -   -           -    -     -  - *************************
C Steve Gibbons Mon Mar 19 15:12:41 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C The file with logical number LU is already open!                   C
C TDRARD then reads in the real array F( N1, N2 ) from file LU with  C
C format flag IFORM.                                                 C
C Currently, IFORM must be 1, for format (5(1PE16.7)).               C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N1        : First dimension of array F.                        C
C     N2        : Second dimension of array F.                       C
C     LU        : Logical unit number of file.                       C
C     IFORM     : Format of data stored in file.                     C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Real array of dimensions ( N1, N2 ).               C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE TDRARD( N1, N2, LU, IFORM, F )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, LU, IFORM
      REAL    F( N1, N2 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I1, I2
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine SVFRD.'
        PRINT *,' IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C OK, so read SV values ...
C
      DO I2 = 1, N2
        IF ( IFORM.EQ.1 ) READ ( LU, 41 )
     1         ( F( I1, I2 ), I1 = 1, N1 )
      ENDDO
C
 41   FORMAT(5(1PE16.7))
      RETURN
      END
C*********************************************************************

