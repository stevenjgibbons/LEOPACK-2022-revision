C*********************************************************************
C subroutine CHange of ReSolution Error Find *************************
C            --        - -        -     -    *************************
C Steve Gibbons Tue Jan 25 08:25:20 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VL        : Value obtained using the lower of two resolutions. C
C     VH        : Value obtained using the finer of two resolutions. C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Output :-                                                          C
C ======                                                             C
C  Integer                                                           C
C  -------                                                           C
C     IERR      : Returned = 0 if RESERR is (vl-vh)/vl               C
C                 Returned = 1 if vl = 0                             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RESERR    : If vl .ne. 0, RESERR = (vl-vh)/vl, else RESERR =vh C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CHRSEF( VL, VH, IERR, RESERR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
C
      DOUBLE PRECISION VL, VH, RESERR
      INTEGER IERR
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION QUOT, LOW
      PARAMETER ( LOW = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( DABS( VL ).LT.LOW ) THEN
        QUOT = 1.0d0
        IERR = 1
      ELSE
        QUOT = 1.0d0/VL
        IERR = 0
      ENDIF
C
      RESERR = QUOT*( VL - VH )
C
      RETURN
      END
C*********************************************************************

