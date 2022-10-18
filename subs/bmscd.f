C*********************************************************************
C subroutine Banded Matrix Singularity Condition Determine ***********
C            -      -      -           -         -         ***********
C Steve Gibbons Wed Jul  5 15:55:07 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C A matrix, A, stored in LAPACK banded form with KL lower diagonals, C
C KU upper diagonals and KLE lower diagonals is known to have NS     C
C singularities. (The dimension of A is N2).                         C
C                                                                    C
C This means that there are NS non-trivial vectors                   C
C                                                                    C
C  X = ( x_{1}, x_{2}, .... , x_{N2} )^T   such that                 C
C                                                                    C
C  A. X = ( 0, 0, ... , 0 )^T                                        C
C                                                                    C
C This means that the rows of the matrix are linearly dependent.     C
C There are therefore NS vectors                                     C
C                                                                    C
C C = ( c_{1}, c_{2}, .... , c_{N2} )^T                              C
C                                                                    C
C such that A^T. C = ( 0, 0, ... , 0 )^T                             C
C                                                                    C
C Now one the non-zero c_{i} can be arbitrarily set to 1.0d0.        C
C The user must supply the index IADC( is ) [index of arbitrarily    C
C specified coefficients] of this c_{i}. (The user should know       C
C which rows of the matrix will be linearly dependent).              C
C                                                                    C
C The array V, with dimensions (N2, NS) is set to zero except        C
C for a 1.0d0 in each element [ IADC( is ), is ] for is = 1, .. , ns C
C                                                                    C
C The columns IADC( is ) (is = 1, ... , ns) of A are then set to     C
C zero except for a 1.0d0 in the diagonal element.                   C
C                                                                    C
C Following this operation, the matrix should be non-singular        C
C and so we proceed with an LU decomposition using the LAPACK        C
C routine DGBTRF.                                                    C
C                                                                    C
C We then solve for the NS vectors C using the LAPACK routine        C
C DGBTRS - the vectors are returned in the columns of V.             C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     NS        : Number of known singularities.                     C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. Must be equal to KL in this routine.      C
C                                                                    C
C     IPIV      : Work array dim ( N2 ) (LAPACK pivotting info)      C
C     IADC      : Rows for singularity fixing. Dim (NS).             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dimensions ( N1, N2 ) format as above      C
C     V         : Matrix dim ( N2, NS ). Outputs coefficients for    C
C                  linear summation of rows of A matrix.             C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE BMSCD( N1, N2, NS, KL, KU, KLE, IPIV, IADC, A, V )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, NS, KL, KU, KLE, IPIV( N2 ), IADC( NS )
      DOUBLE PRECISION A( N1, N2 ), V( N2, NS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NRHS, I, J, IOP, INFO
      DOUBLE PRECISION DZERO, DONE
      CHARACTER *(1) TRANS
      PARAMETER ( DZERO = 0.0d0, DONE = 1.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
C     . Early exit?
C     .
      IF ( NS.EQ.0 ) RETURN
C     .
C     . Check inputs
C     .
      IF ( N1.NE.(KL+KU+KLE+1) ) THEN
        PRINT *,' Subroutine BMSCD, N1 = ', N1
        PRINT *,' KL  = ', KL
        PRINT *,' KU  = ', KU
        PRINT *,' KLE = ', KLE
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check size of KLE
C     .
      IF ( KL.NE.KLE ) THEN
        PRINT *,' Subroutine BMSCD, KL = ', KL
        PRINT *,' Should be equal to KLE = ', KLE
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check that all NS values in IADC are in range
C     . and distinct.
C     .
      DO I = 1, NS
        IF ( IADC( I ).LT.1 .OR. IADC( I ).GT.N2 ) THEN
          PRINT *,' Subroutine BMSCD.'
          PRINT *,' IADC(',I,') = ', IADC( I )
          PRINT *,' N2 = ', N2
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        DO J = 1, NS
          IF ( I.NE.J .AND. IADC( I ).EQ.IADC( J ) ) THEN
            PRINT *,' Subroutine BMSCD.'
            PRINT *,' IADC(',I,') = IADC(',J,') = ', IADC( I )
            PRINT *,' Program aborted.'
            STOP
          ENDIF
        ENDDO
      ENDDO
C     .
C     . So far, so good.
C     . Now zero the array V.
C     .
      IOP = 0
      CALL MATOP( V, DZERO, N2, NS, IOP )
C     .
C     . Loop around singularities
C     . IOP now tells BMRCOP to act on columns
C     .
      IOP = 2
      DO I = 1, NS
        J = IADC( I )
        V( J, I ) = DONE
        CALL BMRCOP( KL, KU, KLE, N2, IOP, J, A, DZERO, DZERO )
        A( KU + KLE + 1, J ) = DONE
      ENDDO
C     .
C     . We will now try to LU decompose the matrix.
C     . Provided there are no problems other than the
C     . singularities dealt with here, this should
C     . not present a problem here.
C     .
      CALL DGBTRF( N2, N2, KL, KU, A, N1, IPIV, INFO )
C     .
C     . Check for an error from LU decomp.
C     .
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine BMSCD.'
        PRINT *,' The LAPACK subroutine DGBTRF has been'
        PRINT *,' and has returned ',INFO,' as a value of INFO.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . We now need to solve A^T.C = V
C     .
      NRHS  = NS
      TRANS = 'T'
      CALL DGBTRS( TRANS, N2, KL, KU, NRHS, A, N1, IPIV, V, N2,
     1             INFO )
C     .
C     . Check solution occured without error.
C     .
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine BMSCD.'
        PRINT *,' The LAPACK subroutine DGBTRS'
        PRINT *,' has returned ',INFO,' as a value of INFO.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      RETURN
      END
C*********************************************************************
