C*********************************************************************
C subroutine VALMIN and VALMAX UPDATE ********************************
C            ------     ------ ------ ********************************
C Steve Gibbons Tue Feb 27 12:14:18 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C F is an array of dimensions ( N1, N2 ).                            C
C                                                                    C
C VALMIN and VALMAX are respectively a global minimum and a global   C
C maximum value. VALMIN_VALMAX_UPDATE tests every element in F and   C
C if such an element is greater than VALMAX, VALMAX is modified and  C
C similarly with VALMIN.                                             C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N1        : First dimension of F                               C
C     N2        : Second dimension of F                              C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Function Dim ( N1, N2 )                            C
C                                                                    C
C     VALMIN    : Global minimum to be tested.                       C
C     VALMAX    : Global maximum to be tested.                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VALMIN_VALMAX_UPDATE( N1, N2, F, VALMIN, VALMAX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2
      REAL    F( N1, N2 ), VALMIN, VALMAX
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I1, I2
      REAL    FBIG, FVAL
      PARAMETER ( FBIG = 1.0e8 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO I2 = 1, N2
        DO I1 = 1, N1
          FVAL = F( I1, I2 )
          IF ( FVAL.GT.VALMAX .AND. FVAL.LT.FBIG ) VALMAX = FVAL
          IF ( FVAL.LT.VALMIN ) VALMIN = FVAL
        ENDDO
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
