C*********************************************************************
C subroutine LINE ReaD ***********************************************
C            ---- -  - ***********************************************
C Steve Gibbons Tue Jan  9 09:48:37 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Reads in a line from a file LU and if the first character is       C
C the 1 character string ESCPCH, then a subsequent line is read.     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LU        : Number of input file.                              C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     LINE      : String (*) to be read.                             C
C     ESCPCH    : String (1): Miss out line if ESCPCH found.         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LINERD( LU, LINE, ESCPCH )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER        LU
      CHARACTER *(*) LINE
      CHARACTER *(1) ESCPCH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
 80   FORMAT(A)
 50   CONTINUE
      READ ( LU, 80 ) LINE
      IF ( LINE(1:1).EQ.ESCPCH(1:1) ) GOTO 50
C
      RETURN
      END
C*********************************************************************
