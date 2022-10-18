C*********************************************************************
C subroutine POWer of TWO checker ************************************
C            ---      ---         ************************************
C Steve Gibbons 11.4.97        					     C
C____________________________________________________________________C
C Checks an integer INPUT and returns a .TRUE. logical variable POT  C
C if and only if INPUT is a power of 2.				     C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     INPUT	: Integer to be tested				     C
C____________________________________________________________________C
C Output :-                                                          C
C ======                                                             C
C  Logical 							     C
C  ------- 							     C
C     POT	: Power Of Two ? 				     C
C____________________________________________________________________C
C*********************************************************************
      SUBROUTINE POWTWO(INPUT,POT)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INPUT
      LOGICAL POT
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITWO,N,IREM
      PARAMETER (ITWO=2)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C First of all put N = INPUT so that INPUT won't be altered
      N=INPUT
C Then check that N is atleast equal to 2 - otherwise
C it is obviously not a power of 2.
      IF ( N.LT.2 ) THEN
         POT=.FALSE.
         RETURN
      ENDIF
C Begin our iteration ....
 500  CONTINUE
      IREM = MOD ( N , ITWO )
      IF ( IREM.NE.0 ) THEN
         POT=.FALSE.
         RETURN
      ENDIF
      IF ( N.EQ.2 ) THEN
         POT=.TRUE.
         RETURN
      ELSE
         N=N/2
         GOTO 500
      ENDIF
      END
C*********************************************************************
