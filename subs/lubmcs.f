C*********************************************************************
C subroutine LU decomposed Banded Matrix + Column Solve **************
C            --            -      -        -      -     **************
C Steve Gibbons Wed Feb 23 14:21:49 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given a banded matrix, A, stored in LAPACK format (i.e. element    C
C A( i, j ) is stored in A( kle + ku + 1 + i - j , j ) - WHICH HAS   C
C ALREADY HAD AN LU DECOMPOSITION PERFORMED ON IT - then LUBMCS will C
C solve the linear matrix equation                                   C
C                                                                    C
C     ( A + C ) x = b                                                C
C                                                                    C
C where the matrix C is zero apart from a single column whose        C
C elements are contained in the vector CVEC.                         C
C                                                                    C
C b is provided in the array RHS: which also returns the solution, x C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. Must be equal to KL in this routine.      C
C                                                                    C
C     ICN       : Number of column which non-zero elements of C      C
C                   occupy.                                          C
C                                                                    C
C     IPIV      : Work array dim ( N2 ) (LAPACK pivotting info)      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dimensions ( N1, N2 ) format as above      C
C                  A is already LU decomposed.                       C
C     RHS       : Vector dim ( N2 ). Corresponds to (b) in above eq. C
C                 Solution X is returned in RHS.                     C
C     CVEC      : Vector which constitutes the non-zero elements     C
C                  of the matrix C. (Dimension (N2) )                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LUBMCS( N1, N2, KL, KU, KLE, ICN, IPIV,
     1                   A, RHS, CVEC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, ICN, IPIV( N2 )
      DOUBLE PRECISION A( N1, N2 ), RHS( N2 ), CVEC( N2 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NRHS, INFO, I
      DOUBLE PRECISION DLOW, FAC, VY, VZ, QUOT
      PARAMETER ( DLOW = 1.0d-8 )
      CHARACTER *(1) TRANS
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check inputs
C
      IF ( N1.NE.(KL+KU+KLE+1) ) THEN
        PRINT *,' Subroutine LUBMCS, N1 = ', N1
        PRINT *,' KL  = ', KL
        PRINT *,' KU  = ', KU
        PRINT *,' KLE = ', KLE
        STOP
      ENDIF
C
C Check size of KLE
C
      IF ( KL.NE.KLE ) THEN
        PRINT *,' Subroutine LUBMCS, KL = ', KL
        PRINT *,' Should be equal to KLE = ', KLE
        STOP
      ENDIF
C
      TRANS = 'N'
      NRHS  = 1
C
C -------------------------------------------------------- C
C The Sherman - Morrison formula says if                   C
C Ay = b and Az = u then x = y - [ vy/(1 + vz) ] z         C
C where  ( A + UV ) x = b                                  C
C and UV( i, j ) = u(i)v(j)                                C
C -------------------------------------------------------- C
C
C     . V is a unit vector with only V( ICN ) non zero
C     . So, simply solve the equations
C     . A.y = RHS and A.z = CVEC
C     .
      CALL DGBTRS( TRANS, N2, KL, KU, NRHS, A, N1,
     1             IPIV, RHS, N2, INFO )
C     .
C     . Check that solution has occured without error
C     .
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine LUBMCS.'
         PRINT *,' The LAPACK subroutine DGBTRS'
         PRINT *,' has returned ',INFO,' as a value of'
         PRINT *,' INFO solving A.y = RHS .'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
      CALL DGBTRS( TRANS, N2, KL, KU, NRHS, A, N1,
     1             IPIV, CVEC, N2, INFO )
C     .
C     . Check that solution has occured without error
C     .
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine LUBMCS.'
         PRINT *,' The LAPACK subroutine DGBTRS'
         PRINT *,' has returned ',INFO,' as a value of'
         PRINT *,' INFO solving A.z = CVEC .'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . VY, VZ are easy to evaluate ...
C     .
      VY = RHS( ICN )
      VZ = CVEC( ICN )
C     .
      QUOT = 1.0d0+VZ
      IF ( DABS( QUOT ).LT.DLOW ) THEN
        PRINT *,' Subroutine LUBMCS.'
        PRINT *,' 1 + vz = ', QUOT
        PRINT *,' Division by zero imminent.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      FAC = VY/QUOT
C     .
C     . Now calculate X from the Sherman Morrison
C     . formula ....
C     .
      DO I = 1, N2
        RHS( I ) = RHS( I ) - FAC*CVEC( I )
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
