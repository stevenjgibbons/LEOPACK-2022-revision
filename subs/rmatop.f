C*********************************************************************
C subroutine Real MATrix OPeration ***********************************
C Steve Gibbons 23.4.97 Does operation on a two-dimensional array.   C
C                                                Can set equal to a  C
C                       constant; multiply by a constant or have a   C
C                       constant added to it.                        C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IOP	: Type of operation to be done.                      C
C                    WHOLE MATRIX OPERATIONS                         C
C                  IOP=0  -->  Each element of the matrix = CONST    C
C                  IOP=1  -->  Each el. is multiplied by CONST       C
C                  IOP=2  -->  Each el. is added to CONST            C
C     NDIM1     : First dimension of the matrix.	             C
C     NDIM2     : Second dimension of the matrix.	             C
C  Real                                                              C
C  ----                                                              C
C     RMAT	: Matrix with dimension ( NDIM1, NDIM2 )             C
C     CONST     : Real constant.                                     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RMATOP ( RMAT, CONST, NDIM1, NDIM2, IOP)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM1, NDIM2, IOP
      REAL    RMAT( NDIM1, NDIM2 ), CONST
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, J
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First do case IOP=0 
      IF ( IOP.EQ.0 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               RMAT ( I, J) = CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
C Now do case IOP=1
      IF ( IOP.EQ.1 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               RMAT ( I, J) = RMAT ( I, J)*CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
C Now do case IOP=2
      IF ( IOP.EQ.2 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               RMAT ( I, J) = RMAT ( I, J) + CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
      PRINT *,' Subroutine RMATOP. IOP must be 0,1 or 2.'
      STOP
      END
C*********************************************************************
