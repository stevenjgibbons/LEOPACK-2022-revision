C*********************************************************************
C subroutine Fintie Difference Coefficient matrix INVert *************
C            -      -          -                  ---    *************
C Steve Gibbons Wed Sep 15 09:07:59 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C For a given value of H (the distance between adjacent and equally  C
C spaced radial grid nodes) fdcinv calculates the coefficients       C
C which multiply the values of f( x_j )                              C
C       [ i - nln .LE. j .LE. i - nrn ] in order to calculate the    C
C n^{th} derivative of f at x_i.                                     C
C                                                                    C
C Coefm is a square matrix with the dimensions ( NCFM, NCFM).        C
C                                                                    C
C If NLN and NRN are the numbers of grid nodes to the left and       C
C right of node j, where the derivative is required, then            C
C NCFM must be atleast ( NLN + NRN + 1 ).                            C
C                                                                    C
C This will return coefficients for calculating derivatives up to    C
C the ( NLN + NRN )^{th} deriv.                                      C
C                                                                    C
C If we are doing central derivatives, i.e. NBN = NLN = NRN,         C
C then the NBN necessary to give the (ND)^{th} derivative to an      C
C accuracy of H^{IACC} (with IACC an even integer) is                C
C                                                                    C
C NBN = ( ND + 1 )/2 + IACC/2 - 1                                    C
C                                                                    C
C Lower derivatives than ( ND + 1 )/2*2 will be accurate to          C
C a higher power of H than the ND^{th}. Subsequent careful calls     C
C must be made to fdcinv in order that all deriv.s are correct       C
C to the same order.                                                 C
C                                                                    C
C On returning from fdcinv, the coefficient are stored in            C
C COEFM such that if:                                                C
C                                                                    C
C f^{ ND }( x_i ) = \sum_{j = i - nln}^{j = i + nrn} \left(          C
C                   COEFM( ND + 1, j - i + nln + 1)*f( x_j )         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     H         : Distance between radial grid nodes.                C
C     COEFM     : Dimension ( NCFM, NCFM).                           C
C     WORK      : Workspace array for LAPACK inversion routine.      C
C                 Dimension ( NCFM )                                 C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NLN       : Number of left nodes. This is how many             C
C                 nodes to the left of node i can be used            C
C                 in order to evaluate a derivative of a             C
C                 function of f at x_i.                              C
C                 For instance if NLH = 2 then the values of         C
C                 f( x_{i-2} ), f( x_{i-1} ) and f( x_i )            C
C                 may be used to evaluate a deriv of f( x_i ).       C
C                                                                    C
C     NRN       : Number of right nodes. This is how many            C
C                 nodes to the right of node i can be used           C
C                 in order to evaluate a derivative of a             C
C                 function of f at x_i.                              C
C                 For instance if NRH = 2 then the values of         C
C                 f( x_{i+2} ), f( x_{i+1} ) and f( x_i )            C
C                 may be used to evaluate a deriv of f( x_i ).       C
C                                                                    C
C     NBN       : Number of banded nodes. This is a maximum for      C
C                 NLN and NRN.                                       C
C                                                                    C
C     NCFM      : Order of the square matrix, COEFM.                 C
C                 Must be atleast NLN + NRN + 1.                     C
C                                                                    C
C     IPCM      : Work array for LAPACK routines to perform          C
C                 pivotting in the matrix inversion.                 C
C                 Dimension ( NCVM )                                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FDCINV ( H, NLN, NRN, NBN, NCFM, COEFM, IPCM, WORK)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NLN, NRN, NBN, NCFM, IPCM( NCFM )
      DOUBLE PRECISION H, COEFM( NCFM, NCFM ), WORK( NCFM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NCF, I, NDER, INODE, MNLN, INFO, ICOL, IROW
      DOUBLE PRECISION DZERO, LOW, FAC
      PARAMETER ( DZERO = 0.0d0, LOW = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check that H is Non-zero
C
      IF ( H.LT.LOW ) THEN
         PRINT *,' Subroutine FDCINV: H = ', H
         PRINT *,' This will lead to a singular matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C Check bounds of integer parameters
C
      IF ( NLN.GT.NBN .OR. NRN.GT.NBN ) THEN
         PRINT *,' Subroutine FDCINV: NBN = ', NBN
         PRINT *,' NLN = ', NLN,'. NRN = ',NRN
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C Check bounds of COEFM
C
      NCF = NLN+NRN+1
      IF ( NCFM.LT.NCF ) THEN
         PRINT *,' Subroutine FDCINV: NCFM = ', NCFM
         PRINT *,' NCF = ',NCF
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C____________________________________________________________________C
C Parameters are ok so let's zero COEFM
C
      I = 0
      CALL MATOP( COEFM, DZERO, NCFM, NCFM, I )
C
C (nder+1) is the number of the matrix column being filled in.
C inode is the number of the matrix row being filled in.
C
      MNLN = (-1)*NLN
      DO NDER = 0, NLN + NRN
        ICOL = NDER + 1
        DO INODE = MNLN, NRN
          IROW = INODE + NLN + 1
          IF ( NDER.EQ.0 ) THEN
            COEFM( IROW, ICOL )  = 1.0d0
          ELSE
            FAC = DBLE( INODE )*H/DBLE( NDER )
            COEFM( IROW, ICOL ) = COEFM( IROW, ICOL-1 )*FAC
          ENDIF
        ENDDO
      ENDDO
C
C Ok - this matrix is now ready for inversion -
C For this we use the LAPACK routines DGETRF and DGETRI
C First perform LU decomposition
C
      CALL DGETRF( NCF, NCF, COEFM, NCFM, IPCM, INFO )
C
C     . Check that LU decomposition has gone without
C     . problem.
C     .
C
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine FDCINV.'
         PRINT *,' The LAPACK subroutine DGETRF has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in LU decomposition of COEFM matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Now compute the inverse with the LAPACK routine
C     . DGETRI.
C     .
      CALL DGETRI( NCF, COEFM, NCFM, IPCM, WORK, NCFM, INFO )
C     .
C     . Check that inversion has gone without problem.
C     .
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine FDCINV.'
         PRINT *,' The LAPACK subroutine DGETRI has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in inversion of COEFM matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
      RETURN
      END
C*********************************************************************
