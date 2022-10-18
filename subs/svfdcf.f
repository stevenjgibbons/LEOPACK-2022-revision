C*********************************************************************
C subroutine Solution Vector Finite Difference Coefficients Form *****
C            -        -      -      -          -            -    *****
C Steve Gibbons Fri Oct 22 09:33:36 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  If XARR is an array of length ( NR ) such that the j^{th}         C
C  element is the value of x_j, then SVFDCF builds an array          C
C  SVFDC of dimension ( NFDCM, NR, NDRVM+1, NDCS ) such that for a   C
C  node number, j, the ND^{th} derivative of radial function         C
C  given f_{IH} ( x ) will be given by                               C
C                                                                    C
C  f_{IH}^{ND}( x_j ) =                                              C
C         \sum_{i=LN}^{RN} SVFDC ( IRAD, j, ND+1, K ) f_{IH} ( x_i ) C
C                                                                    C
C  where LN (the left node)  = MAX( NLMR, j - NBN ) and              C
C        RN (the right node) = MIN( NRMC, j + NBN ),                 C
C                                                                    C
C  IRAD = i - j + NBN + 1 and K = MHP( ih ).                         C
C                                                                    C
C  NDCS is the number of distinct sets of coefficients required for  C
C  different types of harmonics (in generally will be considerably   C
C  smaller than the number of harmonics, NH).                        C
C                                                                    C
C  The arrays MHIBC and MHOBC instruct SVFDCF how to manipulate      C
C  the finite difference coefficients at the boundaries.             C
C                                                                    C
C  MHIBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( ih ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                        where L = LARR( ih )                        C
C                                                                    C
C  Similarly, at the outer boundary:-                                C
C                                                                    C
C  MHOBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( ih ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                        where L = LARR( ih )                        C
C                                                                    C
C  The elements of this array are filled in from j = NLMR            C
C  to j = NRMR ( number of the left most node and number of the      C
C  right most node ) - other rows are left unreferred to.            C
C  This is incase a higher order derivative is required for          C
C  central nodes than boundary nodes; in which case SVFDCF must      C
C  be called for the remaining nodes with modified parameters.       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
C                                                                    C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C     NLMR      : This is the lowest j for which the terms are       C
C                  calculated for SVFDC( i, j, ND+1, K )             C
C     NRMR      : This is the highest j for which the terms are      C
C                  calculated for SVFDC( i, j, ND+1, K )             C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). See above.                     C
C     MHOBC     : Dimension ( NDCS ). See above.                     C
C     LARR      : Spherical harmonic degree, L. Dim. ( NDCS ).       C
C                 This value is only useful for harmonics            C
C                 when calculating derivatives for magnetic fields.  C
C                 If LARR( ih ) = -1, the harmonic is ignored        C
C                 completely.                                        C
C                                                                    C
C     NCFM      : Leading order of working coefficient matrix.       C
C                 Must be atleast (2*NBN + 1) where NBN is the       C
C                 maximum number of nodes on either side of the      C
C                 central node.                                      C
C     NFDCM     : Leading order of the array SVFDC.                  C
C                 This must be atleast (2*NBN + 1)                   C
C     NDRVS     : Number of derivatives required.                    C
C                  This will be limited by the available bandwidth.  C
C                                                                    C
C                  Let NLCS = NLMR - 1    and let                    C
C                      NRCS = NR - NRMR                              C
C                                                                    C
C                  Now, let I = MIN( NLCS, NRCS) + NBN               C
C                                                                    C
C                  then NDRVS must be no greater than I.             C
C                  This is checked for.                              C
C     NDRVM     : Maximum number of derivatives allowed.             C
C                                                                    C
C     IWORK     : Integer work array. Dimension ( NCFM )             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                                                                    C
C     COEFM1    : Coefficient work array. Dimension ( NCFM, NCFM )   C
C     COEFM2    : Coefficient work array. Dimension ( NCFM, NCFM )   C
C     WORK1     : Coefficient work array. Dimension ( NCFM )         C
C     WORK2     : Coefficient work array. Dimension ( NCFM )         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVFDCF( NR, NDCS, NBN, NLMR, NRMR, MHIBC, MHOBC,
     1                   LARR, NCFM, NFDCM, NDRVS, NDRVM, XARR,
     2                   IWORK, SVFDC, COEFM1, COEFM2, WORK1, WORK2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NDCS, NBN, NLMR, NRMR, MHIBC( NDCS ), MHOBC( NDCS ),
     1        LARR( NDCS ), NCFM, NFDCM, NDRVS, NDRVM, 
     1        IWORK( NCFM )
      DOUBLE PRECISION XARR( NR ), SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                 COEFM1( NCFM, NCFM ), COEFM2( NCFM, NCFM ),
     2                 WORK1( NCFM ), WORK2( NCFM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NDER, IRAD, NLCS, NRCS, NLN, NRN, IDCS, L,
     1        NNDS, INDS, I, INODE, NSNIB, NSNOB, NALF, NARF,
     2        IIBC, IOBC, ND1
C
C nsnib is the number of special nodes on the inner boundary
C nsnob is the number of special nodes on the outer boundary
C
      DOUBLE PRECISION DZERO, X0, EMMULT, FAC
      PARAMETER ( DZERO = 0.0d0 )
C
      LOGICAL OCHNGE
C
C ochnge is .TRUE. when the boundary conditions
C play a part in the finite difference coefficients
C and .FALSE. otherwise
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     .
C     . Check the values of integers ...
C     .
      IF ( NDRVS.GT.NDRVM ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NDRVS  = ',NDRVS
         PRINT *,' NDRVM  = ',NDRVM
         STOP
      ENDIF
C     .
      NLCS = NLMR - 1
      NRCS = NR - NRMR
C     . 
C     . Check that sufficient points are allowed
C     . for the derivatives ...
C     . 
      IF ( (NR-1).LT.(NBN+1) ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NBN  = ', NBN
         PRINT *,' NLMR = ', NLMR
         PRINT *,' NRMR = ', NRMR
         PRINT *,' Insufficient nodes for differencing.'
         STOP
      ENDIF
C     . 
      IF ( NLCS.LT.0 ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NLMR = ', NLMR
         STOP
      ENDIF
C     . 
      IF ( NRCS.LT.0 ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NRMR = ', NRMR
         STOP
      ENDIF
C     .
      I = MIN( NLCS, NRCS) + NBN
C     .
      IF ( NDRVS.GT.I ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' You have requested deriv.s to order ',NDRVS
         PRINT *,' At one node, you have only', I
         PRINT *,' side nodes to use for differencing.'
         STOP
      ENDIF
C     .
      I = 2*NBN + 1
C     .
      IF ( NCFM.LT.I ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NCFM = ', NCFM
         PRINT *,' NBN  = ', NBN
         STOP
      ENDIF
C     .
      IF ( NFDCM.LT.I ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NFDCM = ', NFDCM
         PRINT *,' NBN  = ', NBN 
         STOP
      ENDIF
C     .
      IF ( NLMR.GT.NRMR ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NLMR = ', NLMR
         PRINT *,' NRMR = ', NRMR
         STOP
      ENDIF
C     .
C     . Whew ... all input parameters seem to be o.k.
C     . Loop around the requested nodes
C     .
      DO IRAD = NLMR, NRMR
C      .
C      . Loop around the different 'harmonic forms'
C      .
       DO IDCS = 1, NDCS
C
C Check for flag to ignore this harmonic
C
        IF ( LARR( IDCS ).EQ.-1 ) GOTO 50
C
C We now need to check which boundary condition is
C required. Make sure that it is valid.
C
        IF ( MHIBC( IDCS ).LT.1 .AND. MHIBC( IDCS ).GT.7 ) THEN
          PRINT *,' Subroutine SVFDCF.'
          PRINT *,' MHIBC(',IDCS,') = ', MHIBC( IDCS )
          STOP
        ENDIF
C
C O.k. inner b.c. is fine. Now need to 
C see how many points this effects.
C
        IF ( MHIBC( IDCS ).EQ.1 ) THEN
          NSNIB = 0
        ENDIF
C
        IF ( MHIBC( IDCS ).EQ.2 .OR. MHIBC( IDCS ).EQ.3 .OR.
     1       MHIBC( IDCS ).EQ.6 .OR. MHIBC( IDCS ).EQ.7     ) THEN
          NSNIB = 1
        ENDIF
C
        IF ( MHIBC( IDCS ).EQ.4 .OR. MHIBC( IDCS ).EQ.5 ) THEN
          NSNIB = 2
        ENDIF
C
        IF ( MHOBC( IDCS ).LT.1 .AND. MHOBC( IDCS ).GT.7 ) THEN
          PRINT *,' Subroutine SVFDCF.'
          PRINT *,' MHOBC(',IDCS,') = ', MHOBC( IDCS )
          STOP
        ENDIF
C
C O.k. outer b.c. is fine. Now need to
C see how many points this effects.
C
        IF ( MHOBC( IDCS ).EQ.1 ) THEN
          NSNOB = 0
        ENDIF
C
        IF ( MHOBC( IDCS ).EQ.2 .OR. MHOBC( IDCS ).EQ.3 .OR.
     1       MHOBC( IDCS ).EQ.6 .OR. MHOBC( IDCS ).EQ.7     ) THEN
          NSNOB = 1 
        ENDIF
C
        IF ( MHOBC( IDCS ).EQ.4 .OR. MHOBC( IDCS ).EQ.5 ) THEN
          NSNOB = 2 
        ENDIF
C
        DO NDER = 0, NDRVS
          ND1 = NDER + 1
          DO I = 1, NFDCM
            SVFDC( I, IRAD, ND1, IDCS ) = DZERO
          ENDDO
        ENDDO
C
C irad is the node for which we want to calculate
C our coefficients
C
        NLCS = IRAD - 1
        NRCS = NR - IRAD
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
C We must work out how many nodes are affected
C to the left and the right. NALF and NARF
C are respectively the number of affected nodes
C to the left and right.
C
        IF ( (IRAD-NLN).GT.NSNIB ) NALF = 0
        IF ( (IRAD-NLN).EQ.NSNIB ) NALF = 1
        IF ( (IRAD-NLN).LT.NSNIB ) NALF = 2
C
        IF ( (IRAD+NRN).LT.(NR+1-NSNOB) ) NARF = 0
        IF ( (IRAD+NRN).EQ.(NR+1-NSNOB) ) NARF = 1
        IF ( (IRAD+NRN).GT.(NR+1-NSNOB) ) NARF = 2
C
        IF ( NALF.EQ.0 .AND. NARF.EQ.0 ) THEN
          OCHNGE = .FALSE.
        ELSE
          OCHNGE = .TRUE.
        ENDIF
        IF ( .NOT. OCHNGE ) GOTO 51
C       .
C       . OK - we need to form a matrix COEFM2 such that
C       . the correct coeffcients are given when
C       . COEFM1 is multiplied by COEFM2
C       .
        L    = LARR( IDCS )
        IIBC = MHIBC( IDCS )
        IOBC = MHOBC( IDCS )
C       .
        CALL LDGNMF( NR, NNDS, NALF, NARF, L, IIBC, IOBC, NCFM,
     1             XARR, COEFM2, COEFM1, WORK1, WORK2, IWORK )
C       .
 51     CONTINUE
        X0 = XARR( IRAD )
        DO INDS = 1, NNDS
          INODE = IRAD - NLN - 1 + INDS
          WORK1( INDS ) = XARR( INODE )
        ENDDO
C
C Now ready to calculate the coefficients
C
        CALL GFDCFD( X0, WORK1, NNDS, COEFM1, NCFM, 
     1               IWORK, WORK2 )
C
C coefm matrix should now contain the coeff.s
C
        IF ( OCHNGE ) THEN
C        .
C        . Our coefficients are modified
C        . by the boundary conditions
C        .
         DO NDER = 0, NDRVS
          ND1 = NDER + 1
          DO INDS = 1, NNDS
            INODE = INDS - NLN + NBN
            FAC = EMMULT( ND1, INDS, NCFM, NCFM, NNDS,
     1                      COEFM1, COEFM2 )
            SVFDC( INODE, IRAD, ND1, IDCS ) = FAC
          ENDDO
         ENDDO
        ELSE
C        .
C        . Our coefficients are not modified
C        . by the boundary conditions
C        .
         DO NDER = 0, NDRVS
          ND1 = NDER + 1
          DO INDS = 1, NNDS
            INODE = INDS - NLN + NBN
            SVFDC( INODE, IRAD, ND1, IDCS ) = 
     1                      COEFM1( ND1, INDS )
          ENDDO
         ENDDO
C        .
        ENDIF
C
 50    CONTINUE
       ENDDO
C      .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
