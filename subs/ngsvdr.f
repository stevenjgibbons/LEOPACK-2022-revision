C*********************************************************************
C subroutine Non-uniform Grid Solution Vector DeRivative *************
C            -           -    -        -      - -        *************
C Steve Gibbons Thu Sep 23 08:24:08 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C V is a solution vector with NH (= inarr(3) ) harmonic radial       C
C functions and NR (= inarr(2) ) grid nodes in each function.        C
C The position of the j^{th} node of the i^{th} harmonic is given by C
C INDFUN( j, i, INARR) and the radial value is given by              C
C  XARR( j ) - an array which is passed to the routine FDCMBD        C
C in order to calculate the finite difference coefficients, FDCM.    C
C XARR itself is not referenced by NGSVDR.                           C
C NGSVDR returns the radial derivatives 1, ..., IHD of radial func.  C
C IH evaluated at node IR.                                           C
C                                                                    C
C IFORMF = INARR(1) should be either 3 or 4 since this is the        C
C arbitrarily spaced mesh version of the code.                       C
C                                                                    C
C NBN is the maximum number of nodes on                              C
C either side which maybe used in central differences.               C
C For instance, if to calculate the derivative of                    C
C f at r_j, you may use the values of f at r = r_{j-2}, r_{j-1},     C
C r_j, r_{j+1} and r_{j+2} then NBN = 2. The value of IHD is checked C
C only for being positive and no greater than NDVDS (the             C
C number of the highest derivative for which coefficients            C
C are stored by the array FDCM), as FDCM must be calculated in       C
C advance by a call to FDCMBD which checks NDVDS against             C
C the physical restrictions imposed by the value of NBN, ILNC and    C
C IRNC; which are respectively the furthest left and furthest right  C
C nodes which may be used to calculate derivatives and are the       C
C same values which must be supplied to FDCMBD as NLMC and NRMC.     C
C                                                                    C
C NBN must be as supplied to FDCMBD.                                 C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Number of radial grid node.                        C
C     IH        : Number of radial function (harmonic).              C
C     NBN       : Number of bounding nodes. See above.               C
C     IHD       : Highest derivative requested.                      C
C     NFDCM     : Leading dimension of FDCM. At least (2*NBN+1)      C
C     NR        : Number of radial grid nodes in each function.      C
C     NDVDS     : Number of highest derivative for which             C
C                  coefficients are stored by the array FDCM.        C
C     ILNC      : Number of left most node which can be used to      C
C                     calculate derivatives.                         C
C     IRNC      : Number of right most node which can be used to     C
C                     calculate derivatives.                         C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NRR. Must be consistent with NR.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V         : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C     DERV      : Derivatives. Dim ( * ) but length atleast IHD      C
C                  DERV( i ) is returned containing the value of     C
C                  the i^[th} derivative.                            C
C                                                                    C
C     FDCM      : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVS ).                   C
C                   Array is generated by the routine fdcmbd         C
C                 See documentation for FDCMBD for details.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NGSVDR ( V, IR, IH, NBN, IHD, NFDCM, NR, NDVDS,
     1                    DERV, ILNC, IRNC, INARR, FDCM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, IH, NBN, IHD, NFDCM, NR, NDVDS, ILNC, IRNC,
     1        INARR( * )
      DOUBLE PRECISION V( * ), DERV( * ), 
     1                 FDCM( NFDCM, NR, NDVDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION COEF
      INTEGER INODE, ILN, IRN, INDFUN, ID, IND, NRR, IFORMF, IK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameters ...
C
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
C
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine NGSVDR.'
         PRINT *,' INARR( 2 ) = ', NRR
         PRINT *,' NR = ', NR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
         PRINT *,' Subroutine NGSVDR.'
         PRINT *,' INARR( 1 ) = ', IFORMF
         PRINT *,' This is an irregular grid routine.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IR.LT.ILNC .OR. IR.GT.IRNC ) THEN
        PRINT *,' Subroutine NGSVDR.'
        PRINT *,' IR   = ', IR
        PRINT *,' ILNC = ', ILNC
        PRINT *,' IRNC = ', IRNC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IHD.LT.1 .OR. IHD.GT.NDVDS ) THEN
        PRINT *,' Subroutine NGSVDR.'
        PRINT *,' IHD   = ', IHD
        PRINT *,' NDVDS = ', NDVDS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Calculate furthest left and furthest right mode
C     . to be used to form derivative.
C     .
      ILN = MAX( ILNC, IR - NBN )
      IRN = MIN( IRNC, IR + NBN )
C
      DO ID = 1, IHD
        DERV( ID ) = 0.0d0
        DO INODE = ILN, IRN
          IK = INODE - IR + NBN + 1
          COEF = FDCM( IK, IR, ID )
          IND = INDFUN( INODE, IH, INARR )
          DERV( ID ) = DERV( ID ) + COEF*V( IND )
        ENDDO
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

