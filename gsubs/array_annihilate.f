C*********************************************************************
C subroutine ARRAY ANNIHILATE ****************************************
C            ----- ---------- ****************************************
C Steve Gibbons Tue Mar 20 10:24:45 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C F( N1, N2 ) is a real array.                                       C
C                                                                    C
C F( I1, I2 ) is a function evaluated at XXX, YYY                    C
C                                                                    C
C    where XXX = C1V1 + REAL( I1 - 1 )*DXXX                          C
C               with DXXX = (C1V2 - C1V1)/REAL( N1 - 1 )             C
C  and                                                               C
C          YYY = C2V1 + REAL( J2 - 1 )*DYYY                          C
C               with DYYY = (C2V2 - C2V1)/REAL( N2 - 1 )             C
C                                                                    C
C ARRAY_ANHILIATE loops around all the points in F and, if           C
C the XXX and YYY coordinates lie in the respective open intervals   C
C (C1B1,C1B2) and (C2B1,C2B2) then the value of F is set to          C
C 1.0e10.                                                            C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N1        : First dimension of the array F                     C
C     N2        : Second dimension of the array F                    C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Data array with dimensions ( N1, N2 )              C
C     C1V1      : See above.                                         C
C     C1V2      : See above.                                         C
C     C2V1      : See above.                                         C
C     C2V2      : See above.                                         C
C     C1B1      : See above.                                         C
C     C1B2      : See above.                                         C
C     C2B1      : See above.                                         C
C     C2B2      : See above.                                         C
C                                                                    C
C*********************************************************************
      SUBROUTINE ARRAY_ANNIHILATE( N1, N2, F, C1V1, C1V2, C2V1, C2V2,
     1                             C1B1, C1B2, C2B1, C2B2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2
      REAL    F( N1, N2 ), C1V1, C1V2, C2V1, C2V2,
     1        C1B1, C1B2, C2B1, C2B2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I1, I2
      REAL    XX, YY, FBIG, RDX, RDY
      PARAMETER ( FBIG = 1.0e10 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      RDX = (C1V2 - C1V1)/REAL( N1 - 1 )
      RDY = (C2V2 - C2V1)/REAL( N2 - 1 )
C
      DO I2 = 1, N2
        YY = C2V1 + REAL( I2 - 2 )*RDY
        DO I1 = 1, N1
          XX = C1V1 + REAL( I1 - 2 )*RDX
C         .
          IF (    XX.GT.C1B1 .AND. XX.LT.C1B2 .AND.
     1            YY.GT.C2B1 .AND. YY.LT.C2B2           )  THEN
            F( I1, I2 ) = FBIG
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C
