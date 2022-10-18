C*********************************************************************
C subroutine Finite Difference Coefficient Matrix BuilD **************
C            -      -          -           -      -   - **************
C Steve Gibbons Tue Sep 21 09:25:54 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  If XARR is an array of length ( NR ) such that the j^{th}         C
C  element is the value of x_j, then FDCMBD builds an array          C
C  FDCM of dimension ( NFDCM, NR, NDRVS ) such that for a given      C
C  node number, j, the ND^{th} derivative of a function f( x )       C
C  will be given by                                                  C
C                                                                    C
C  f^{ND}( x_j ) = \sum_{i=LN}^{RN} FDCM( IRAD, j, ND ) f( x_i )     C
C                                                                    C
C  where LN (the left node)  = MAX( NLMC, j - NBN ) and              C
C        RN (the right node) = MIN( NRMC, j + NBN )                  C
C                                                                    C
C  and IRAD = i - j + NBN + 1                                        C
C                                                                    C
C  NLMC and NRMC are respectively the left most and right most       C
C  nodes (columns) which may be used to obtain a difference formula. C
C  In most matrix applications NLMC = 1 and NRMC = NR, although      C
C  when differentiating a vector it may be necessary to omit an      C
C  extreme point; for instance when this would result in a division  C
C  by zero.                                                          C
C                                                                    C
C  The elements of this array are filled in from j = NLMN            C
C  to j = NRMN ( number of the left most node and number of the      C
C  right most node ) - other rows are left unreferred to.            C
C  This is incase a higher order derivative is required for          C
C  central nodes than boundary nodes; in which case FDCMBD must      C
C  be called for the remaining nodes with modified parameters.       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     FDCM      : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM ).                   C
C                                                                    C
C     COEFM     : Coefficient work array. Dimension ( NCFM, NCFM )   C
C     WORK1     : Coefficient work array. Dimension ( NCFM )         C
C     WORK2     : Coefficient work array. Dimension ( NCFM )         C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C     NLMN      : This is the lowest j for which the terms are       C
C                  calculated for FDCM( i, j, ND )                   C
C     NRMN      : This is the highest j for which the terms are      C
C                  calculated for FDCM( i, j, ND )                   C
C     NLMC      : This is the lowest i for which the terms are       C
C                  calculated for FDCM( i, j, ND )                   C
C     NRMC      : This is the highest i for which the terms are      C
C                  calculated for FDCM( i, j, ND )                   C
C                                                                    C
C     NCFM      : Leading order of working coefficient matrix.       C
C                 Must be atleast (2*NBN + 1) where NBN is the       C
C                 maximum number of nodes on either side of the      C
C                 central node.                                      C
C     NFDCM     : Leading order of the array FDCM.                   C
C                 This must be atleast (2*NBN + 1)                   C
C     NDRVS     : Number of derivatives required.                    C
C                  This will be limited by the available bandwidth.  C
C                                                                    C
C                  Let NLCS = NLMN - NLMC and let                    C
C                      NRCS = NRMC - NRMN                            C
C                                                                    C
C                  Now, let I = MIN( NLCS, NRCS) + NBN               C
C                                                                    C
C                  then NDRVS must be no greater than I.             C
C                  This is checked for.                              C
C     NDRVM     : Maximum number of derivatives required.            C
C                                                                    C
C     IWORK     : Integer work array. Dimension ( NCFM )             C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FDCMBD( NR, NBN, NLMN, NRMN, NLMC, NRMC, NCFM,
     1                   NFDCM, NDRVS, NDRVM, IWORK, XARR, FDCM,
     2                   COEFM, WORK1, WORK2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NBN, NLMN, NRMN, NCFM, NFDCM, NDRVS, NDRVM, 
     1        IWORK( NCFM ), NLMC, NRMC
      DOUBLE PRECISION XARR( NR ), FDCM( NFDCM, NR, NDRVM ),
     1                 COEFM( NCFM, NCFM ), WORK1( NCFM ),
     2                 WORK2( NCFM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NDER, IRAD, NLCS, NRCS, NLN, NRN,
     1        NNDS, INDS, I, INODE
      DOUBLE PRECISION DZERO, X0
      PARAMETER ( DZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     .
C     . Check the values of integers ...
C     .
      IF ( NDRVS.GT.NDRVM ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NDRVS  = ',NDRVS
         PRINT *,' NDRVM  = ',NDRVM
         STOP
      ENDIF
C     .
      NLCS = NLMN - NLMC
      NRCS = NRMC - NRMN
C     . 
C     . Check that sufficient points are allowed
C     . for the derivatives ...
C     . 
      IF ( (NRMC-NLMC).LT.(NBN+1) ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NBN  = ', NBN
         PRINT *,' NLMC = ', NLMC
         PRINT *,' NRMC = ', NRMC
         PRINT *,' Insufficient nodes for differencing.'
         STOP
      ENDIF
C     . 
      IF ( NLCS.LT.0 ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NLMN = ', NLMN
         PRINT *,' NLMC = ', NLMC
         STOP
      ENDIF
C     . 
      IF ( NRCS.LT.0 ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NRMN = ', NRMN
         PRINT *,' NRMC = ', NRMC
         STOP
      ENDIF
C     .
      I = MIN( NLCS, NRCS) + NBN
C     .
      IF ( NDRVS.GT.I ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' You have requested deriv.s to order ',NDRVS
         PRINT *,' At one node, you have only', I
         PRINT *,' side nodes to use for differencing.'
         STOP
      ENDIF
C     .
      I = 2*NBN + 1
C     .
      IF ( NCFM.LT.I ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NCFM = ', NCFM
         PRINT *,' NBN  = ', NBN
         STOP
      ENDIF
C     .
      IF ( NFDCM.LT.I ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NFDCM = ', NFDCM
         PRINT *,' NBN  = ', NBN 
         STOP
      ENDIF
C     .
      IF ( NLMN.GT.NRMN ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NLMN = ', NLMN
         PRINT *,' NRMN = ', NRMN
         STOP
      ENDIF
C     .
C     . Whew ... all input parameters seem to be o.k.
C     . Now loop around the requested rows
C     .
      DO IRAD = NLMN, NRMN
C
        DO NDER = 1, NDRVS
          DO I = 1, NFDCM
            FDCM( I, IRAD, NDER ) = DZERO
          ENDDO
        ENDDO
C
C irad is the node for which we want to calculate
C our coefficients
C
        NLCS = IRAD - NLMC
        NRCS = NRMC - IRAD
C
C we wish to calculate NLN (number of left nodes)
C and NRN ( number of right nodes )
C NNDS ( total number of nodes) is then NLN + NRN + 1 ...
C
        NLN = MIN( NBN, NLCS )
        NRN = MIN( NBN, NRCS )
C
        NNDS = NLN + NRN + 1
C
        X0 = XARR( IRAD )
        DO INDS = 1, NNDS
          INODE = IRAD - NLN - 1 + INDS
          WORK1( INDS ) = XARR( INODE )
        ENDDO
C
C Now ready to calculate the coefficients
C
        CALL GFDCFD( X0, WORK1, NNDS, COEFM, NCFM, 
     1               IWORK, WORK2 )
C
C coefm matrix should now contain the coeff.s
C
        DO NDER = 1, NDRVS
          DO INDS = 1, NNDS
            INODE = INDS - NLN + NBN
            FDCM( INODE, IRAD, NDER ) = COEFM( NDER+1, INDS )
          ENDDO
        ENDDO
C
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
