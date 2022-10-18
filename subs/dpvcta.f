C*********************************************************************
C subroutine Double Precision Vector Coriolis Term Add ***************
C            -      -         -      -        -    -   ***************
C Steve Gibbons 17.6.99                                              C
C____________________________________________________________________C
C Adds the term to the double precision vector resulting from        C
C the curl ( k cross v ) term in the momentum equation.              C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Level of harmonics.                                C
C     NDIM      : Length of vector, SV.                              C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     NR        : Number of radial grid nodes                        C
C     NH        : Number of harmonics (all types)                    C
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C     MHM       : Dimension ( NH )                                   C
C     MHC       : Dimension ( NH )                                   C
C                                                                    C
C  mht, mhl, mhm and mhc define the spherical harmonics present      C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     ORD       : *(2). Determines the order of accuracy. Options :  C
C             SS - Strictly second order                             C
C             SF - Strictly fourth order                             C
C             O5 - Optimum accuracy for bandwidth 5; this gives      C
C                  Fourth order accuracy for 1st and 2nd derivatives C
C                  and second order accuracy for 3rd and 4th der.s   C
C             O7 - Optimum accuracy for bandwidth 7; this gives      C
C                  Sixth order accuracy for 1st and 2nd derivatives  C
C                  and fourth order accuracy for 3rd and 4th der.s   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RVEC      : Receiving vector for output.                       C
C                  Dim ( NDIM ) with NDIM = NH * NR                  C
C     OLDSV     : Old Solution Vector Dimensions                     C
C                  Dim ( NDIM ) with NDIM = NH * NR                  C
C     FAC       : Multiplication factor for term.                    C
C     RI        : Radius of inner core.                              C
C     RO        : Radius of outer core.                              C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Gaussian weights evaluated by the routine          C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LHMAX + 1 )*( LHMAX + 2 )/2 , NTHMAX }    C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C____________________________________________________________________C
C Working variables :-                                               C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     QST       : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector                           C
C                  Has dimensions (  LH*(LH+2) , 3).                 C
C     TMPSV1    : Old Solution Vector Dimensions                     C
C     TMPSV2    : Old Solution Vector Dimensions                     C
C                  Dim ( NDIM ) with NDIM = NH * NR                  C
C     VF        : Vector Function. An array of dimensions            C
C                  ( NPHPTS, NTHPTS, 3) which contain the R, THETA   C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IPHI, ITHETA, 2 ) is the Theta      C
C                  compontent of the vector at (iphi, itheta).       C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPVCTA ( LH, NDIM, NTHPTS, NPHPTS, NR, NH, MHT, MHL,
     1                    MHM, MHC, ORD, RVEC, OLDSV, RI, RO, FAC,
     2                    QST, TMPSV1, TMPSV2, VF, FTF1, FTF2, FTF3,
     3                    GAUX, GAUW, PA, DPA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NDIM, NTHPTS, NPHPTS, NR, NH, MHT( NH ), MHL( NH ),
     1        MHM( NH ), MHC( NH )
      CHARACTER *(2) ORD
      DOUBLE PRECISION RVEC( NDIM ), OLDSV( NDIM ), TMPSV1( NDIM ),
     1                 RI, RO, GAUX( NTHPTS ), FTF1( 2*NPHPTS ),
     2                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ), FAC,
     3                 VF( NPHPTS, NTHPTS, 3), TMPSV2( NDIM ),
     4                 GAUW( NTHPTS )
      DOUBLE PRECISION QST( LH*(LH + 2), 3 ),
     1                 PA( (LH+1)*(LH+2)/2 , NTHPTS ),
     2                 DPA( (LH+1)*(LH+2)/2 , NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRN, IH, L, M, ICS, NOHARM, INDSHC, IND, IPTF, IVLMF
      DOUBLE PRECISION RAD, D0F, D1F, D2F, D3F, D4F, SQRLL1, H, TOL
      PARAMETER ( TOL = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of NDIM against NH and NR etc.
C
      IF ( NDIM.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPVCTA, bad array size'
        PRINT *,' NDIM = ',NDIM
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FAC ).LT.TOL ) RETURN
C
C Set value of H
C
      H = ( RO - RI )/DBLE( NR - 1 )
C
      DO IRN = 1, NR
         RAD = RI + H*DBLE( IRN - 1 )
C        .
C        . Find qst coefficients for the velocity
C        . at this node.  
C        .
         IPTF = 1
         IVLMF = 1
         CALL SV2QST ( LH, NDIM, NR, NH, MHT, MHL, MHM, MHC, IRN,
     1                 IVLMF, IPTF, ORD, OLDSV, RI, RO, QST )
C        .
C        . evaluate velocity in space at this node
C        .
         CALL QST2VA ( QST, VF, GAUX, PA, DPA, FTF1, FTF2, FTF3,
     1                 LH, LH, NTHPTS, NTHPTS, NPHPTS, NPHPTS )
C        .
C        . evaluate Coriolis force in space
C        .
         CALL VFCOR ( NTHPTS, NTHPTS, NPHPTS, NPHPTS, VF, GAUX )
C        .
C        . re-convert vector back to qst values
C        .
         CALL VF2QSA ( QST, VF, GAUX, GAUW, PA, DPA, FTF1,
     1                 FTF2, FTF3, LH, LH, NTHPTS, NTHPTS,
     2                 NPHPTS, NPHPTS )
C        .
C        . We now have the QST values corresponding to
C        . K x V. We store them in the appropriate 
C        . places in the TMPSV1 and TMPSV2 arrays.
C        .
         DO IH = 1, NH
            L   = MHL( IH )
            M   = MHM( IH )
            ICS = MHC( IH )
            NOHARM = INDSHC( L, M, ICS )
            IND = ( IRN - 1 )*NH + IH
C           .
C           . Put scaloidal parts into poloidal TMPSV1
C           . Put spheroidal parts into poloidal TMPSV2
C           .
            IF ( MHT(IH).EQ.1 ) THEN
               TMPSV1( IND ) = QST( NOHARM, 1 )
               TMPSV2( IND ) = QST( NOHARM, 2 )
            ENDIF
C           .
C           . Put toroidal parts into toroidal TMPSV1
C           .
            IF ( MHT(IH).EQ.2 ) THEN
               TMPSV1( IND ) = QST( NOHARM, 3 )
            ENDIF
         ENDDO
C        .
      ENDDO
C
C Now K x V for all grid nodes is stored in TMPSV1 and TMPSV2
C Now take the curl and put resulting curl ( k x v )
C into RVEC.
C
      DO IRN = 1, NR
        RAD = RI + H*DBLE( IRN - 1 )
        DO IH = 1, NH
           L = MHL( IH )
           IF ( MHT( IH ).EQ.1 ) THEN
             IND = ( IRN - 1 )*NH + IH
             RVEC( IND ) = RVEC( IND ) + FAC*TMPSV1( IND )/RAD
             CALL SVDERF ( NDIM, NH, IH, NR, IRN, TMPSV2, ORD,
     1                     D0F, D1F, D2F, D3F, D4F, RI, RO )
             RVEC( IND ) = RVEC( IND ) - FAC*(D1F+D0F/RAD)/SQRLL1(L)
           ENDIF
           IF ( MHT( IH ).EQ.2 ) THEN
             IND = ( IRN - 1 )*NH + IH
             RVEC( IND ) = RVEC( IND ) - FAC*TMPSV1( IND )/SQRLL1(L)
           ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
