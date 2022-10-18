C*********************************************************************
C subroutine Banded Matrix WooDbury Formula Solve ********************
C            -      -      -  -     -       -     ********************
C Steve Gibbons Sat Mar  4 12:37:40 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given a banded matrix, A, stored in LAPACK format (i.e. element    C
C A( i, j ) is stored in A( kle + ku + 1 + i - j , j ) .... ),       C
C then BMWDFS solves the matrix equation                             C
C                                                                    C
C     ( A + [UVM] ) x = b                                            C
C                                                                    C
C where the matrix [UVM] is the product of matrices U ( N2, NP )     C
C and the transpose of V ( N2, NP ).                                 C
C                                                                    C
C NP is the number of additional rows AND columns with which the     C
C banded matrix A is to be altered. For example if the linear        C
C system above is to be modified from the banded formation by the    C
C addition of 2 columns ( say columns 5 and 7 ) and 3 rows           C
C ( say rows 1, 7 and 8) then NP would be set to  5 (=2+3).          C
C                                                                    C
C  (PROCEDURE FOR COLUMN ALTERATION .... )                           C
C                                                                    C
C Now the 1st column of U must contain the elements of the additnl.  C
C elements to column 5 and the 1st column of V must be all zero      C
C except for the number 1.0d0 in row 5. Similarly,                   C
C the 2nd column of U must contain the elements of the additional    C
C elements to column 7 and the 2nd column of V must be all zero      C
C except for the number 1.0d0 in row 7.                              C
C                                                                    C
C  (PROCEDURE FOR ROW ALTERATION .... )                              C
C                                                                    C
C However, the 3rd column of U must be zero except for a 1.0d0 in    C
C row 1. The element V( i, 3 ) will contain the value of element     C
C (1, i ) of the matrix [UVH].                                       C
C                                                                    C
C Similarly, column 4 ( 5 ) of U must be zero except for a 1.0d0     C
C in row 7 ( 8 ). Column 4 ( 5 ) of V must contain the elements of   C
C the new row 7 ( 8 ).                                               C
C                                                                    C
C Hope this is sufficient to illustrate the principle.               C
C                                                                    C
C This routine differs from the earlier version BMWDSL in that       C
C it does not automatically perform an LU decomposition.             C
C This is to allow for repeated solution of a matrix equation        C
C where the matrix has previously been LU decomposed.                C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     NP        : Number of additional operations to matrix.         C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. Must be equal to KL in this routine.      C
C                                                                    C
C     IPIV      : Work array dim ( N2 ) (LAPACK pivotting info)      C
C     IPIVH     : Work array dim ( NP ) (LAPACK pivotting info)      C
C                                                                    C
C     ILUDF     : LU decomposition flag.                             C
C                 =1 if matrix has just been formed and LU decomp.   C
C                    is required.                                    C
C                 =2 if LU decomp. has already been done and we      C
C                    simply require a further solution.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dimensions ( N1, N2 ) format as above      C
C     RHS       : Vector dim ( N2 ). Corresponds to (b) in above eq. C
C                 Solution X is returned in RHS.                     C
C     U         : Matrix dim ( N2, NP ). See above for guidance.     C
C     V         : Matrix dim ( N2, NP ). See above for guidance.     C
C                                                                    C
C     HMAT      : Matrix dim ( NP, NP ). Workspace.                  C
C     VTY       : Matrix dim ( NP, 1  ). Workspace.                  C
C     HVTY      : Matrix dim ( NP, 1  ). Workspace.                  C
C     UHVTY     : Matrix dim ( N2, 1  ). Workspace.                  C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE BMWDFS( N1, N2, NP, KL, KU, KLE, IPIV, IPIVH, A,
     1                   RHS, U, V, HMAT, VTY, HVTY, UHVTY, ILUDF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, NP, KL, KU, KLE, IPIV( N2 ), IPIVH( NP ),
     1        ILUDF
      DOUBLE PRECISION A( N1, N2 ), RHS( N2 ), U( N2, NP ),
     1                 V( N2, NP ), HMAT( NP, NP ), VTY( NP, 1 ),
     2                 HVTY( NP, 1 ), UHVTY( N2, 1 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NRHS, INFO, I, MM, NN, KK
      DOUBLE PRECISION ALPHA, BETA
      CHARACTER *(1) TRANSA, TRANSB
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      INFO = 0
C
C Check inputs
C
      IF ( N1.NE.(KL+KU+KLE+1) ) THEN
        PRINT *,' Subroutine BMWDFS, N1 = ', N1
        PRINT *,' KL  = ', KL
        PRINT *,' KU  = ', KU
        PRINT *,' KLE = ', KLE
        STOP
      ENDIF
C
C Check size of KLE
C
      IF ( KL.NE.KLE ) THEN
        PRINT *,' Subroutine BMWDFS, KL = ', KL
        PRINT *,' Should be equal to KLE = ', KLE
        STOP
      ENDIF
C     .
      IF ( ILUDF.NE.1 .AND. ILUDF.NE.2 ) THEN
        PRINT *,' Subroutine BMWDFS.'
        PRINT *,' ILUDF = ', ILUDF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C First perform LU decomposition upon A.
C
      IF ( ILUDF.EQ.1 )
     1           CALL DGBTRF( N2, N2, KL, KU, A, N1, IPIV, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine BMWDFS.'
         PRINT *,' The LAPACK subroutine DGBTRF has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in LU decomposition of (A) matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C We now follow the prescription of Numerical Recipes in C
C (hereafter referred to as NR)
C 2nd edition ISBN 0-521-43108-5 page 76 to implement
C Woodbury formula to solve system.
C
C First we solve the system A. y = RHS
C
C     .
      TRANSA = 'N'
      NRHS = 1
C     .
      CALL DGBTRS( TRANSA, N2, KL, KU, NRHS, A, N1,
     1             IPIV, RHS, N2, INFO )
C     .
C     . Check that solution has occured without error
C     .
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine BMWDFS.'
         PRINT *,' The LAPACK subroutine DGBTRS'
         PRINT *,' has returned ',INFO,' as a value of'
         PRINT *,' INFO solving A.y = RHS.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . RHS now contains the vector NR refers to as y
C     . If NP is 0 then we can make an early exit!
C     .
      IF ( NP.EQ.0 ) RETURN
C     .
C     . Now we solve the auxiliary problems A.Z = U
C     . and so NRHS must be set to NP
C     .
      NRHS = NP
C     .
      CALL DGBTRS( TRANSA, N2, KL, KU, NRHS, A, N1,
     1             IPIV, U, N2, INFO )
C     .
C     . Check that solution has occured without error
C     .
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine BMWDFS.'
         PRINT *,' The LAPACK subroutine DGBTRS'
         PRINT *,' has returned ',INFO,' as a value of'
         PRINT *,' INFO solving A.Z = U.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . U now contains the vectors NR refers to as Z
C     . We must now construct matrix HMAT in stages
C     . First calculate V^T.Z: for this we use DGEMM.
C     .
      TRANSA = 'T'
      TRANSB = 'N'
      ALPHA  = 1.0d0
      BETA   = 0.0d0
C
      MM = NP
      NN = NP
      KK = N2
C
      CALL DGEMM ( TRANSA, TRANSB, MM, NN, KK,
     1             ALPHA, V, N2, U, N2, BETA, HMAT, NP)
C
C     .
C     . Now add the identity matrix to HMAT
C     .
      DO I = 1, NP
        HMAT( I, I ) = HMAT( I, I ) + 1.0d0
      ENDDO
C     .
C     . Now invert the square matrix HMAT
C     . First perform LU decomposition,
C     .
      CALL DGETRF( NP, NP, HMAT, NP, IPIVH, INFO )
C     .
C     . Check that LU decomposition has gone without
C     . problem.
C     .
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine BMWDFS.'
         PRINT *,' The LAPACK subroutine DGETRF has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in LU decomposition of (H) matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Now compute the inverse with the LAPACK routine
C     . DGETRI. For the workspace, we will supply
C     . UHVTY as this has plenty of space ( N2 elements )
C     . and has yet to be used.
C     .
      CALL DGETRI( NP, HMAT, NP, IPIVH, UHVTY, N2, INFO )
C     .
C     . Check that inversion has gone without problem.
C     .
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine BMWDFS.'
         PRINT *,' The LAPACK subroutine DGETRI has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in inversion of (H) matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Our matrix HMAT now represents the matrix H in
C     . NR on page 76 eqn. (2.7.19)
C     . We now proceed with a suite of matrix multiplications
C     . to evaluate the solutions as in eqn (2.7.21)
C     .
C     . First multiply V transpose by y ( RHS )
C     .
      TRANSA = 'T'
      TRANSB = 'N'
      ALPHA  = 1.0d0
      BETA   = 0.0d0
C
      MM = NP
      NN = 1
      KK = N2
C
      CALL DGEMM ( TRANSA, TRANSB, MM, NN, KK,
     1             ALPHA, V, N2, RHS, N2, BETA, VTY, NP)
C      
C     .
C     . Now premultiply resulting matrix (VTY) by
C     . HMAT (calculated above)
C     .
      TRANSA = 'N'
      TRANSB = 'N'
      ALPHA  = 1.0d0
      BETA   = 0.0d0
C
      MM = NP
      NN = 1
      KK = NP
C
      CALL DGEMM ( TRANSA, TRANSB, MM, NN, KK,
     1             ALPHA, HMAT, NP, VTY, NP, BETA, HVTY, NP)
C
C     .
C     . Now premultiply resulting matrix (HVTY) by
C     . HMAT (calculated above)
C     .
      TRANSA = 'N'
      TRANSB = 'N'
      ALPHA  = 1.0d0
      BETA   = 0.0d0
C
      MM = N2
      NN = 1
      KK = NP
C
      CALL DGEMM ( TRANSA, TRANSB, MM, NN, KK,
     1             ALPHA, U, N2, HVTY, NP, BETA, UHVTY, N2)
C
C  OK - so now evaluate solution!!!
C
      DO I = 1, N2
        RHS( I ) = RHS( I ) - UHVTY( I, 1 )
      ENDDO
C
      RETURN
      END
C*********************************************************************
