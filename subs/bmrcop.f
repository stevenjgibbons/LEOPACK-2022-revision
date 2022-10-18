C*********************************************************************
C subroutine Banded Matrix Row or Column OPeration *******************
C            -      -      -      -      --        *******************
C Steve Gibbons 3.12.97                                              C
C____________________________________________________________________C
C                                                                    C
C If a double precision matrix ABAND is stored in LAPACK format      C
C i.e. with the element a(i,j) being stored in                       C
C  ABAND( KLE + KU + 1 + i - j, j ),                                 C
C then BMRCOP will take either a row or a column of A and            C
C multiply it by ALPHA and add BETA.                                 C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C     NDIM	: Dimension of matrix.                               C
C     IRC	: = 1 to act on ROW                                  C
C                 = 2 to act on COLUMN                               C
C     NRC	: Number of row or column.                           C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     ABAND     : Banded matrix. Dimensions ( NDIM1, NDIM2 ) with    C
C                  NDIM1 = KL + KU + KLE + 1                         C
C                  NDIM2 = Dimension of problem (NDIM ).             C
C     ALPHA	: First multiplies current matrix entry by ALPHA     C
C     BETA 	: Then adds on BETA.                                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE BMRCOP( KL, KU, KLE, NDIM, IRC, NRC, ABAND, 
     1                   ALPHA, BETA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KL, KU, KLE, NDIM, IRC, NRC
      DOUBLE PRECISION ABAND( KLE + KU + KL + 1, NDIM ) , ALPHA, BETA
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IC, IR, IROW, IND1, IND2
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
C     .                    First deal with the case acting on ROWS ...
      IF ( IRC.EQ.1 ) THEN
         IR = NRC
         IND1 = MAX( 1, IR - KL )
         IND2 = MIN( IR + KU, NDIM )
         DO IC = IND1, IND2
            IROW = KLE + KU + 1 + IR - IC
            ABAND( IROW, IC ) = ABAND( IROW, IC )*ALPHA + BETA
         ENDDO
         RETURN
      ENDIF
c
C     .                    First deal with the case acting on COLUMNS.
      IF ( IRC.EQ.2 ) THEN
         IC = NRC
         IND1 = MAX( 1, IC - KU )
         IND2 = MIN( IC + KL, NDIM )
         DO IR = IND1, IND2
            IROW = KLE + KU + 1 + IR - IC
            ABAND( IROW, IC ) = ABAND( IROW, IC )*ALPHA + BETA
         ENDDO
         RETURN
      ENDIF
C
      PRINT *,' Subroutine BMRCOP '
      PRINT *,' IRC must be set to either 1 or 2.'
      PRINT *,' Program aborted.'
      STOP
C
      END
C*********************************************************************
