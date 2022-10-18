C*********************************************************************
C subroutine Adapted Solution Vector DeRivative **********************
C            -       -        -      - -        **********************
C Steve Gibbons Mon Oct 25 15:17:59 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C V is a solution vector with NH (= inarr(3) ) harmonic radial       C
C functions and NR (= inarr(2) ) grid nodes in each function.        C
C The position of the j^{th} node of the i^{th} harmonic is given by C
C INDFUN( j, i, INARR) and the radial value is given by              C
C  XARR( j ) - an array which is passed to the routine FDCMBD        C
C in order to calculate the finite difference coefficients, FDCM.    C
C XARR itself is not referenced by ASVDR.                            C
C ASVDR returns the radial derivatives 0, 1, ..., IHD of radial      C
C function IH evaluated at node IR.                                  C
C                                                                    C
C IFORMF = INARR(1) should be either 3 or 4 since this is the        C
C arbitrarily spaced mesh version of the code.                       C
C                                                                    C
C NBN is the maximum number of nodes on                              C
C either side which maybe used in central differences.               C
C For instance, if to calculate the derivative of                    C
C f at r_j, you may use the values of f at r = r_{j-2}, r_{j-1},     C
C r_j, r_{j+1} and r_{j+2} then NBN = 2. The value of IHD is checked C
C only for being positive and no greater than NDRVS (the             C
C number of the highest derivative for which coefficients            C
C are stored by the array SVFDC), as SVFDC must be calculated in     C
C advance by a call to SVFDCF which checks NDRVS against             C
C the physical restrictions imposed by the value of NBN.             C
C                                                                    C
C NDRVM restricts the size of NDRVS and is a defining parameter      C
C of the array SVFDCF.                                               C
C                                                                    C
C NBN must be as supplied to SVFDCF.                                 C
C                                                                    C
C IS must be supplied to indicate the finite difference scheme       C
C being employed.                                                    C
C                                                                    C
C ASVDR will use whichever coefficients SVFDC( k, j, nd + 1, K ),    C
C that SVFDCF formed with K = IS.                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Number of radial grid node.                        C
C     IS        : Number of finite difference scheme used.           C
C     IH        : Number of radial function (harmonic).              C
C     NBN       : Number of bounding nodes. See above.               C
C     IHD       : Highest derivative requested.                      C
C     NFDCM     : Leading dimension of SVFDC. At least (2*NBN+1)     C
C     NR        : Number of radial grid nodes in each function.      C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NRR. Must be consistent with NR.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
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
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ASVDR ( V, IR, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                   NDRVM, DERV, INARR, SVFDC, NDCS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, IS, IH, NBN, IHD, NFDCM, NR, NDRVS, NDRVM,
     1        INARR( * ), NDCS
      DOUBLE PRECISION V( * ), DERV( * ), 
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION COEF
      INTEGER INODE, ILN, IRN, INDFUN, ID, IND, NRR, IFORMF, IK, ID1
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
         PRINT *,' Subroutine ASVDR.'
         PRINT *,' INARR( 2 ) = ', NRR
         PRINT *,' NR = ', NR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
         PRINT *,' Subroutine ASVDR.'
         PRINT *,' INARR( 1 ) = ', IFORMF
         PRINT *,' This is an irregular grid routine.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IR.LT.1 .OR. IR.GT.NR ) THEN
        PRINT *,' Subroutine ASVDR.'
        PRINT *,' IR   = ', IR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IHD.LT.0 .OR. IHD.GT.NDRVS .OR. NDRVS.GT.NDRVM ) THEN
        PRINT *,' Subroutine ASVDR.'
        PRINT *,' IHD   = ', IHD
        PRINT *,' NDRVS = ', NDRVS
        PRINT *,' NDRVM = ', NDRVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IS.LT.1 .OR. IS.GT.NDCS ) THEN
        PRINT *,' Subroutine ASVDR.'
        PRINT *,' IS    = ', IS
        PRINT *,' NDCS  = ', NDCS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Calculate furthest left and furthest right mode
C     . to be used to form derivative.
C     .
      ILN = MAX(  1, IR - NBN )
      IRN = MIN( NR, IR + NBN )
C
      DO ID = 0, IHD
        ID1 = ID + 1
        DERV( ID1 ) = 0.0d0
        DO INODE = ILN, IRN
          IK = INODE - IR + NBN + 1
          COEF = SVFDC( IK, IR, ID1, IS )
          IND = INDFUN( INODE, IH, INARR )
          DERV( ID1 ) = DERV( ID1 ) + COEF*V( IND )
        ENDDO
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

