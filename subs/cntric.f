C*********************************************************************
C subroutine CouNTeR Increment and Check *****************************
C            -  -- - -             -     *****************************
C Steve Gibbons Mon Jan 10 10:16:42 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Adds 1 to the number IC.                                           C
C If the resulting number is greater than ICMAX then the logical     C
C flag OK is set to .FALSE.                                          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IC       : Value of counter variable.                          C
C     ICMAX    : Maximum permitted value of IC.                      C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     OK       : Unaltered if IC.le.ICMAX.                           C
C                Set to .FALSE. if IC.gt.ICMAX.                      C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CNTRIC( IC, ICMAX, OK )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IC, ICMAX
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      IC = IC + 1
      IF ( IC.GT.ICMAX ) OK = .FALSE.
      RETURN
      END
C*********************************************************************
