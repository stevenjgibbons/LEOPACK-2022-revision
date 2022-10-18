C*********************************************************************
C subroutine Spherical Harmonic Scalar File WriTe ********************
C            -         -        -      -    -  -  ********************
C Steve Gibbons Wed Oct 10 14:45:05 WEST 2001                        C
C____________________________________________________________________C
C                                                                    C
C Writes an 'nview' format file containing spherical harmonic        C
C coefficients.                                                      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input Variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum degree, l, of spherical harmonics.         C
C     LU        : Logical unit number.                               C
C     IFORM     : Format flag for data in file.                      C
C                                                                    C
C     Currently, only option is IFORM = 1 for a format (5d16.7)      C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : Filename - length not specified                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SHC       : Spherical harmonic coefficients.                   C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHSFWT( LH, LU, FNAME, SHC, IFORM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LU, LH, IFORM
      DOUBLE PRECISION SHC( * )
      CHARACTER *( * ) FNAME
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      CHARACTER *(80) LINE
      INTEGER IWR, NH, I
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IFORM.NE.1 ) THEN
         PRINT *,' Subroutine SHSFWT.'
         PRINT *,' IFORM = ', IFORM
         PRINT *,' Not a recognised option.'
         STOP
      ENDIF
C
      IWR = 3
C open file ...
      CALL FOPEN ( LU, FNAME, IWR )
C write a blank line
      LINE = '                '
      WRITE ( LU, 88 ) LINE
C read second line ind read LH ...
      WRITE ( LINE, * ) LH
      WRITE ( LU, 88 ) LINE
      NH = LH * ( LH + 2 )
      IF ( IFORM.EQ.1 ) WRITE ( LU, 89 ) ( SHC( I ), I = 1, NH )
      CALL FCLOSE ( LU, FNAME, FNAME )

 88   FORMAT (A80)
 89   FORMAT (5d16.7)
      RETURN
      END
C*********************************************************************
