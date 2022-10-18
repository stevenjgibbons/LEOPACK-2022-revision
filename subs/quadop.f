C*********************************************************************
C subroutine QUAD dimensional matrix OPeration ***********************
C Steve Gibbons 27.10.99 Does operation on a four-dimensional array. C
C                                                Can set equal to a  C
C                       constant; multiply by a constant or have a   C
C                       constant added to it.                        C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IOP	: Type of operation to be done.                      C
C                  IOP=0  -->  Each element of the matrix = CONST    C
C                  IOP=1  -->  Each el. is multiplied by CONST       C
C                  IOP=2  -->  Each el. is added to CONST            C
C     N1     : First dimension of the matrix.	                     C
C     N2     : Second dimension of the matrix.	                     C
C     N3     : Third dimension of the matrix.	                     C
C     N4     : Fourth dimension of the matrix.	                     C
C  Double Precision                                                  C
C  ----------------                                                  C
C     QUAD	: Matrix with dimension ( N1, N2, N3, N4 )           C
C     CONST     : Double precision constant.                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE QUADOP( QUAD, CONST, N1, N2, N3, N4, IOP)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, N3, N4, IOP
      DOUBLE PRECISION QUAD( N1, N2, N3, N4 ), CONST
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, J, K, L
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First do case IOP=0 
      IF ( IOP.EQ.0 ) THEN
         DO L = 1, N4
           DO K = 1, N3
             DO J = 1, N2
                DO I = 1, N1
                   QUAD( I, J, K, L) = CONST
                ENDDO
             ENDDO
           ENDDO
         ENDDO
         RETURN
      ENDIF
C Now do case IOP=1
      IF ( IOP.EQ.1 ) THEN
         DO L = 1, N4
           DO K = 1, N3
             DO J = 1, N2
                DO I = 1, N1
                   QUAD( I, J, K, L) = QUAD( I, J, K, L)*CONST
                ENDDO
             ENDDO
           ENDDO
         ENDDO
         RETURN
      ENDIF
C Now do case IOP=2
      IF ( IOP.EQ.2 ) THEN
         DO L = 1, N4
           DO K = 1, N3
             DO J = 1, N2
                DO I = 1, N1
                   QUAD( I, J, K, L) = QUAD( I, J, K, L) + CONST
                ENDDO
             ENDDO
           ENDDO
         ENDDO
         RETURN
      ENDIF
      PRINT *,' Subroutine QUADOP. IOP must be 0,1 or 2.'
      STOP
      END
C*********************************************************************

