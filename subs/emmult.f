C*********************************************************************
C function Element of Matrix MULTiplication evaluate *****************
C          -          -      ----                    *****************
C Steve Gibbons Fri Oct 22 17:23:46 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Let the matrices B and C be respectively 'm by k' and 'k by n'     C
C double precision matrices. Then A = B C is an 'm by n' matrix      C
C with a_{ij} = \sum_{l=1,k} b_{il} c_{lj}                           C
C                                                                    C
C EMMULT returns the value of A( I, J)                               C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     I         : Row containing desired element of A.               C
C     J         : Column containing desired element of A.            C
C     LDB       : Leading dimension of matrix B.                     C
C     LDC       : Leading dimension of matrix C.                     C
C     K         : Number of columns of matrix B and                  C
C                 Number of rows of matrix C.                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     B         : First matrix. Dim ( LDB, * )                       C
C     C         : First matrix. Dim ( LDC, * )                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION EMMULT( I, J, LDB, LDC, K, B, C )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER I, J, LDB, LDC, K
      DOUBLE PRECISION EMMULT, B( LDB, * ), C( LDB, * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First check the validity of the integer inputs.
C K must not exceed LDC
C I must not exceed LDB
C     .
      IF ( K.GT.LDC .OR. I.GT.LDB ) THEN
        PRINT *,' Function EMMULT.'
        PRINT *,' K    = ', K,' I    = ', I
        PRINT *,' LDB  = ', LDB,' LDC  = ', LDC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      EMMULT = 0.0d0
C     . 
      DO L = 1, K
        EMMULT = EMMULT + B( I, L )*C( L, J )
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

