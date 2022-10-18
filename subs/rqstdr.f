C*********************************************************************
C subroutine Radial QST DeRivative ***********************************
C            -      --- - -        ***********************************
C Steve Gibbons Tue Sep 28 11:27:29 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Calculates a derivative of a radial function to a harmonic         C
C in the RQST array.                                                 C
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
C                 Given as a function of l, m etc. by INDSHC         C
C                                                                    C
C     ICOMP     : 1 for Q, 2 for S and 3 for T.                      C
C                                                                    C
C     NBN       : Number of bounding nodes. See above.               C
C     IHD       : Highest derivative requested.                      C
C     NFDCM     : Leading dimension of FDCM. At least (2*NBN+1)      C
C     NR        : Number of radial grid nodes in each function.      C
C     LH        : Maximum value of degree, l.                        C
C     NDVDS     : Number of highest derivative for which             C
C                  coefficients are stored by the array FDCM.        C
C     ILNC      : Number of left most node which can be used to      C
C                     calculate derivatives.                         C
C     IRNC      : Number of right most node which can be used to     C
C                     calculate derivatives.                         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RQST      : Dim ( LH*(LH+2), 3, NR )                           C
C                                                                    C
C                 Input array containing scaloidal/spheroidal        C
C                  decomposition of vector.                          C
C                  Has dimensions (  LH*(LH+2) ,3, NR ).             C
C              RQST (l*l+2m,1,I) = q_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,1,I) = q_l^mc(r_i)                       C
C              RQST (l*l+2m,2,I) = s_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,2,I) = s_l^mc(r_i)                       C
C              RQST (l*l+2m,3,I) = t_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,3,I) = t_l^mc(r_i)                       C
C                                                                    C
C                                                                    C
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
      SUBROUTINE RQSTDR ( RQST, IR, IH, NBN, IHD, NFDCM, NR, NDVDS,
     1                    DERV, ILNC, IRNC, FDCM, LH, ICOMP )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, IH, NBN, IHD, NFDCM, NR, NDVDS, ILNC, IRNC, LH,
     1        ICOMP
      DOUBLE PRECISION RQST( LH*(LH+2), 3, NR ), DERV( * ), 
     1                 FDCM( NFDCM, NR, NDVDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION COEF, FAC
      INTEGER INODE, ILN, IRN, ID, IK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameters ...
C
      IF ( IH.LT.1 .OR. IH.GT.(LH*LH + 2*LH) ) THEN
        PRINT *,' Subroutine RQSTDR.'
        PRINT *,' IH   = ', IH
        PRINT *,' LH   = ', LH
        PRINT *,' NH   = ', LH*(LH+2)
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( IR.LT.ILNC .OR. IR.GT.IRNC ) THEN
        PRINT *,' Subroutine RQSTDR.'
        PRINT *,' IR   = ', IR
        PRINT *,' ILNC = ', ILNC
        PRINT *,' IRNC = ', IRNC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IHD.LT.1 .OR. IHD.GT.NDVDS ) THEN
        PRINT *,' Subroutine RQSTDR.'
        PRINT *,' IHD   = ', IHD
        PRINT *,' NDVDS = ', NDVDS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( ICOMP.NE.1 .AND. ICOMP.NE.2 .AND. ICOMP.NE.3 ) THEN
        PRINT *,' Subroutine RQSTDR.'
        PRINT *,' ICOMP = ', ICOMP
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
          FAC  = RQST( IH, ICOMP, INODE )
          DERV( ID ) = DERV( ID ) + COEF*FAC
        ENDDO
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

