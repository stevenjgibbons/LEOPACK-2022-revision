C*********************************************************************
C subroutine Banded Matrix Sherman-Morrison formula SoLve ************
C            -      -      -       -                - -   ************
C Steve Gibbons 30.8.99                                              C
C____________________________________________________________________C
C                                                                    C
C Given a banded matrix, A, stored in LAPACK format (i.e. element    C
C A( i, j ) is stored in A( kle + ku + 1 + i - j , j ) .... ),       C
C then BMSMSL solves the matrix equation                             C
C                                                                    C
C     ( A + RC ) x = b                                               C
C                                                                    C
C where the matrix RC is zero apart from a single row or a single    C
C column whose elements are contained in the vector VRC.             C
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
C     IRCF      : The row/column flag. If IRCF = 1,                  C
C                  the non zero elements of the RC matrix form       C
C                  a single row.       If IRCF = 2,                  C
C                  the non zero elements of the RC matrix form       C
C                  a single column.                                  C
C                                                                    C
C     IRCN      : Number of row/column which non-zero elements of    C
C                  RC occupy.                                        C
C                                                                    C
C                                                                    C
C     IPIV      : Work array dim ( N2 ) (LAPACK pivotting info)      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dimensions ( N1, N2 ) format as above      C
C     RHS       : Vector dim ( N2 ). Corresponds to (b) in above eq. C
C                 Solution X is returned in RHS.                     C
C     VRC       : Vector which constitutes the non-zero elements     C
C                  of the matrix RC. (Dimension (N2) )               C
C                                                                    C
C                 If IRCF is 1 then VRC( j ) contains matrix el.     C
C                  RC( ircn, j )                                     C
C                                                                    C
C                 If IRCF is 2 then VRC( i ) contains matrix el.     C
C                  RC( i, ircn )                                     C
C                                                                    C
C     VEC       : Vector dim. ( N2 ). Does not need to be set on     C
C                                       entry.                       C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE BMSMSL( N1, N2, KL, KU, KLE, IRCF, IRCN, IPIV,
     1                   A, RHS, VRC, VEC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, IRCF, IRCN, IPIV( N2 )
      DOUBLE PRECISION A( N1, N2 ), RHS( N2 ), VRC( N2 ), VEC( N2 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NRHS, INFO, I
      DOUBLE PRECISION ZERO, FAC, VY, VZ
      PARAMETER ( ZERO = 0.0d0 )
      CHARACTER *(1) TRANS
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check inputs
C
      IF ( N1.NE.(KL+KU+KLE+1) ) THEN
        PRINT *,' Subroutine BMSMSL, N1 = ', N1
        PRINT *,' KL  = ', KL
        PRINT *,' KU  = ', KU
        PRINT *,' KLE = ', KLE
        STOP
      ENDIF
C
C Check size of KLE
C
      IF ( KL.NE.KLE ) THEN
        PRINT *,' Subroutine BMSMSL, KL = ', KL
        PRINT *,' Should be equal to KLE = ', KLE
        STOP
      ENDIF
C
      TRANS = 'N'
      NRHS  = 1
C
C First perform LU decomposition upon A.
C
      CALL DGBTRF( N2, N2, KL, KU, A, N1, IPIV, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine BMSMSL.'
         PRINT *,' The LAPACK subroutine DGBTRF has been '
         PRINT *,' and has returned ',INFO,' as a value of '
         PRINT *,' INFO. '
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C -------------------------------------------------------- C
C The Sherman - Morrison formula says if                   C
C Ay = b and Az = u then x = y - [ vy/(1 + vz) ] z         C
C where  ( A + UV ) x = b                                  C
C and UV( i, j ) = u(i)v(j)                                C
C -------------------------------------------------------- C
C
C
C Now consider the case of IRCF = 1 (row addition)
C
      IF ( IRCF.EQ.1 ) THEN
C        .
C        .'construct' the unit vector
C        .
         DO I = 1, N2
           VEC( I ) = ZERO
         ENDDO
         VEC( IRCN ) = 1.0d0
C        .
C        . VEC now corresponds to the 'u' in box above
C        . Now solve for 'y' (using RHS) and 'z'
C        . (using VEC)
C        .
         CALL DGBTRS( TRANS, N2, KL, KU, NRHS, A, N1,
     1             IPIV, RHS, N2, INFO )
C        .
C        . Check that solution has occured without error
C        .
         IF ( INFO.NE.0 ) THEN
            PRINT *,' Subroutine BMSMSL.'
            PRINT *,' The LAPACK subroutine DGBTRS'
            PRINT *,' has returned ',INFO,' as a value of'
            PRINT *,' INFO solving A.y = RHS .'
            PRINT *,' Program aborted.'
            STOP
         ENDIF
C        .
         CALL DGBTRS( TRANS, N2, KL, KU, NRHS, A, N1,
     1             IPIV, VEC, N2, INFO )
C        .
C        . Check that solution has occured without error
C        .
         IF ( INFO.NE.0 ) THEN
            PRINT *,' Subroutine BMSMSL.'
            PRINT *,' The LAPACK subroutine DGBTRS'
            PRINT *,' has returned ',INFO,' as a value of'
            PRINT *,' INFO solving A.z = VEC .'
            PRINT *,' Program aborted.'
            STOP
         ENDIF
C        .
C        . Now calculate VY ( the scalar product of
C        . vectors RHS and VRC )
C        .
         VY = 0.0d0
         DO I = 1, N2
           VY = VY + VRC( I )*RHS( I )
         ENDDO
C        .
C        . Now calculate VZ ( the scalar product of
C        . vectors VEC and VRC )
C        .
         VZ = 0.0d0
         DO I = 1, N2
           VZ = VZ + VRC( I )*VEC( I )
         ENDDO
C        .
C        . Now calculate FAC = [ VY / ( 1 + VZ ) ]
C        .
         FAC = VY/(1.0d0+VZ)
C        .
C        . Now calculate X from the Sherman Morrison
C        . formula ....
C        .
         DO I = 1, N2
           RHS( I ) = RHS( I ) - FAC*VEC( I )
         ENDDO
C        .
         RETURN
      ENDIF
C
C Now consider the case of IRCF = 2 (column addition)
C The proceedure here is a little different ...
C
      IF ( IRCF.EQ.2 ) THEN
C        .
C        . V is a unit vector with only V( IRCN ) non zero
C        . Unlike in the previous case, we do not need
C        . to use this explicitly in any arithmetic.
C        . So, simply solve the equations
C        . A.y = RHS and A.z = VRC
C        .
         CALL DGBTRS( TRANS, N2, KL, KU, NRHS, A, N1,
     1             IPIV, RHS, N2, INFO )
C        .
C        . Check that solution has occured without error
C        .
         IF ( INFO.NE.0 ) THEN
            PRINT *,' Subroutine BMSMSL.'
            PRINT *,' The LAPACK subroutine DGBTRS'
            PRINT *,' has returned ',INFO,' as a value of'
            PRINT *,' INFO solving A.y = RHS .'
            PRINT *,' Program aborted.'
            STOP
         ENDIF
C        .
         CALL DGBTRS( TRANS, N2, KL, KU, NRHS, A, N1,
     1             IPIV, VRC, N2, INFO )
C        .
C        . Check that solution has occured without error
C        .
         IF ( INFO.NE.0 ) THEN
            PRINT *,' Subroutine BMSMSL.'
            PRINT *,' The LAPACK subroutine DGBTRS'
            PRINT *,' has returned ',INFO,' as a value of'
            PRINT *,' INFO solving A.z = VRC .'
            PRINT *,' Program aborted.'
            STOP
         ENDIF
C        .
C        . VY, VZ are easier to evaluate ...
C        .
         VY = RHS( IRCN )
         VZ = VRC( IRCN )
         FAC = VY/(1.0d0+VZ)
C        .
C        . Now calculate X from the Sherman Morrison
C        . formula ....
C        .
         DO I = 1, N2
           RHS( I ) = RHS( I ) - FAC*VRC( I )
         ENDDO
C        .
         RETURN
      ENDIF
C
      PRINT *,' Subroutine BMSMSL. IRCF = ', IRCF
      PRINT *,' Should be either 1 or 2.'
      STOP
      END
C*********************************************************************
