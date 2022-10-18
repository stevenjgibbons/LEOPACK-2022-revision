C*********************************************************************
C subroutine Radial QST Curl Add *************************************
C            -      --- -    -   *************************************
C Steve Gibbons Tue Feb 15 08:17:08 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C If RQST contains a QST decomposition of a vector function          C
C then CRQST is returned containing                                  C
C                                                                    C
C  A * CRQST + B * curl ( RQST )                                     C
C                                                                    C
C There must be NR grid nodes for both functions and both functions  C
C must be stored with harmonics of degree l up to LH.                C
C                                                                    C
C The radial value of the node, j is stored in XARR( j )             C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Maximum value of degree, l.                        C
C     NR        : Number of radial grid nodes in each function.      C
C     M0        : Basic wavenumber. Modes M which do not satisfy     C
C                 MOD( M, M0 ) = 0 are bypassed. This increases      C
C                 the speed in 2.5D calculations.                    C
C                                                                    C
C     MMAX      : Maximum wavenumber. Modes M which do not satisfy   C
C                 M .le. MMAX are ignored.                           C
C                                                                    C
C     NBN       : Number of bounding nodes. See RQSTDR.              C
C                                                                    C
C     NFDCM     : Leading dimension of FDCM. At least (2*NBN+1)      C
C     NDVDS     : Number of highest derivative for which             C
C                  coefficients are stored by the array FDCM.        C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C     ILNC      : Number of left most node which can be used to      C
C                     calculate derivatives.                         C
C     IRNC      : Number of right most node which can be used to     C
C                     calculate derivatives.                         C
C                                                                    C
C  ILNC and IRNC correspond to NLMC and NRMC in the call to FDCMBD.  C
C                                                                    C
C     ICLS      : Set to 1 to zero the array. No action for other    C
C                                                     values.        C
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
C     CRQST     : Dim ( LH*(LH+2), 3, NR ) - Is returned containing  C
C                 the curl of the function in RQST.                  C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     FDCM      : Dim ( NFDCM, NR, NDRVS ). Finite diff. coeff.s     C
C                Must be formed in advance by a call to FDCMBD.      C
C                                                                    C
C     A         : Factor multiplying original values of CRQST.       C
C                                                                    C
C     B         : Factor multiplying values of curl of RQST.         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RQSTCA( LH, NR, M0, MMAX, NBN, NFDCM, NDVDS, ILNR,
     1                   IRNR, ILNC, IRNC, RQST, CRQST, XARR, FDCM,
     2                   ICLS, A, B )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, LH, M0, MMAX, NBN, NFDCM, NDVDS, ILNR, IRNR,
     1        ILNC, IRNC, ICLS
      DOUBLE PRECISION XARR( NR ), A, B,
     1                 FDCM( NFDCM, NR, NDVDS ),
     2                 RQST( LH*(LH+2), 3, NR ),
     3                 CRQST( LH*(LH+2), 3, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IH, L, M, ICS, IR, IQ, IS, IT, IHD, ICOMP
      DOUBLE PRECISION D0F, D1F, DERV(1), ZERO, RAD, SQRLL1
      PARAMETER ( ZERO = 0.0d0, IQ = 1, IS = 2, IT = 3, IHD = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      IF ( M0.EQ.0 ) THEN
        PRINT *,' Subroutine RQSTCA.'
        PRINT *,' M0 = ', M0
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . First zero all elements of CRQST which correspond
C     . to IR with IR.GE.ILNR and IR.LE.IRNC if required.
C     .
      IF ( ICLS.EQ.1 ) THEN
       DO IR = ILNR, IRNR
         DO ICOMP = 1, 3
           DO IH = 1, LH*(LH+2)
             CRQST( IH, ICOMP, IR ) = ZERO
           ENDDO
         ENDDO
       ENDDO
      ENDIF
C     .
      DO IH = 1, LH*(LH+2)
        CALL LMFIND( IH, L, M, ICS )
        IF ( M.GT.MMAX ) GOTO 50
        IF ( MOD( M, M0).NE.0 ) GOTO 50
C       .
C       . curl of scaloidal function q(r)
C       . subtracts -sqrll1( l )q(r)/r from toroidal part
C       .
        DO IR = ILNR, IRNR
          RAD = XARR( IR )
          D0F = RQST( IH, IQ, IR )
          CRQST( IH, IT, IR ) = A*CRQST( IH, IT, IR ) -
     1                          B*SQRLL1( L )*D0F/RAD
        ENDDO
C       .
C       . curl of spheroidal function s(r)
C       . adds ( ds(r)/dr + s(r)/r ) to toroidal part
C       .
        DO IR = ILNR, IRNR
          RAD = XARR( IR )
          D0F = RQST( IH, IS, IR )
          CALL RQSTDR( RQST, IR, IH, NBN, IHD, NFDCM, NR, NDVDS,
     1                 DERV, ILNC, IRNC, FDCM, LH, IS )
          D1F = DERV( 1 )
          CRQST( IH, IT, IR ) = CRQST( IH, IT, IR ) +
     1                                   B*( D1F + D0F/RAD )
        ENDDO
C       .
C       . curl of toroidal function t(r)
C       . subtracts sqrll1( l ) t(r)/r from the scaloidal
C       . part and 
C       . subtracts ( dt(r)/dr + t(r)/r ) from the spheroidal
C       . part
C       .
        DO IR = ILNR, IRNR
          RAD = XARR( IR )
          D0F = RQST( IH, IT, IR )
          CALL RQSTDR( RQST, IR, IH, NBN, IHD, NFDCM, NR, NDVDS,
     1                 DERV, ILNC, IRNC, FDCM, LH, IT )
          D1F = DERV( 1 )
          CRQST( IH, IQ, IR ) = A*CRQST( IH, IQ, IR ) -
     1                          B*SQRLL1( L )*D0F/RAD
          CRQST( IH, IS, IR ) = A*CRQST( IH, IS, IR ) -
     1                          B*( D1F + D0F/RAD )
        ENDDO
C       .
 50   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
