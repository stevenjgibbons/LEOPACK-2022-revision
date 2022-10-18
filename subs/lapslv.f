C*********************************************************************
C subroutine LAPack SoLVer *******************************************
C            ---    - --   *******************************************
C Steve Gibbons 17.7.97                                              C
C____________________________________________________________________C
C Takes the banded matrix ABAND and the RHS vector of dimensions     C
C ( KLE + KL + KU + 1 , NCOLS )                     and              C
C        ( NCOLS )                 respectively  and solves for the  C
C resulting velocity by driving the LAPACK subroutines DGBTRF ( to   C
C do the LU decomposition ) and DGBTRS ( to solve the matrix eqn. )  C
C - the solution vector is returned in RHS.                          C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NCOLS	: Number of cols in matrix.                          C
C     IPIV      : Array needed by the LAPACK routines. Has           C
C                  dimension ( NCOLS )                               C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     ABAND     : The banded convection matrix.                      C
C                  Has dimensions ( NDIM1 , NDIM2 )                  C
C                      where NDIM1 = KLE + KU + KL + 1               C
C                    and NDIM2 = NCOLS                               C
C                                                                    C
C     RHS       : On entry, contains the right hand side to the      C
C                  matrix equation and on exit containus the soln.   C
C                   vector. Dimension ( NCOLS )                      C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LAPSLV ( ABAND, RHS, NCOLS, IPIV, KL, KU, KLE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C

      INTEGER NCOLS, IPIV ( NCOLS ), KL, KU, KLE
      DOUBLE PRECISION ABAND ( KL+KU+KLE+1, NCOLS ), RHS ( NCOLS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
 
      INTEGER NROWS, INFO, LABAND, NRHS
      CHARACTER *1 TRANS
 
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
 
      LABAND = KLE + KL + KU + 1
      NROWS = NCOLS
      NRHS = 1

      IF ( KLE.NE.KL ) THEN
         PRINT *,' Subroutine LAPSLV. '
         PRINT *,' Need KLE = KL to use LAPACK routines.'
         PRINT *,' Program Aborted.'
         STOP
      ENDIF
 
C .......... first perform LU decomposition on ABAND .............

      CALL DGBTRF( NROWS, NCOLS, KL, KU, ABAND, LABAND, IPIV, INFO)

      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine LAPSLV. '
         PRINT *,' The LAPACK subroutine DGBTRF has been '
         PRINT *,' and has returned ',INFO,' as a value of '
         PRINT *,' INFO. '
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C .......... now solve the A.X = B system ...........................
 
      TRANS = 'N'
 
      CALL DGBTRS( TRANS, NCOLS, KL, KU, NRHS, ABAND, LABAND,
     1             IPIV, RHS, NCOLS, INFO )
 
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine LAPSLV. '
         PRINT *,' The LAPACK subroutine DGBTRS has been '
         PRINT *,' and has returned ',INFO,' as a value of '
         PRINT *,' INFO. '
         PRINT *,' Program aborted.'
         STOP
      ENDIF
 
      RETURN
      END
C*********************************************************************
 



