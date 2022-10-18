C*********************************************************************
C subroutine Double Precision Vector Velocity dot gradient of Theta  *
C            -      -         -      -                        -      *
C Steve Gibbons 26.6.99                                              C
C____________________________________________________________________C
C Adds the term to the double precision matrix resulting from        C
C vel_0 . Grad ( theta_0 ) in the heat equation.                     C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Level of harmonics.                                C
C     NDIM      : Length of vector, OLDSV.                           C
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
C                  Dim { ( LH+ 1 )*( LH+ 2 )/2 , NTHPTS }            C
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
C     VF1       : Vector Function. An array of dimensions            C
C                  ( NPHPTS, NTHPTS, 3) which contain the R, THETA   C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IPHI, ITHETA, 2 ) is the Theta      C
C                  compontent of the vector at (iphi, itheta).       C
C     VF2       : Second vector Function. See above                  C
C     SF        : Scalar Function. Dimensions ( NTHPTS, NPHPTS )     C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C     SHC       : Spherical Harmonic Coefficients.                   C
C                  Has the usual ordering.                           C
C                   SHC( l*l ) = P_l^0                               C
C                   SHC( l*l + 2*m - 1 ) = P_l^mc ( m non zero )     C
C                   SHC( l*l + 2*m ) = P_l^ms ( m non zero )         C
C                  Dimension is ( LH*( LH + 2) )                     C
C     DSHC      : The derivatives of the above at the radial nodes.  C
C                  Dimension is ( LH*( LH + 2) )                     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPVVTA ( LH, NDIM, NTHPTS, NPHPTS, NR, NH, MHT, MHL,
     1                    MHM, MHC, ORD, RVEC, OLDSV, RI, RO, FAC, SF,
     2                    QST, VF1, VF2, FTF1, FTF2, FTF3,
     3                    GAUX, GAUW, PA, DPA, SHC, DSHC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NDIM, NTHPTS, NPHPTS, NR, NH, MHT( NH ), MHL( NH ),
     1        MHM( NH ), MHC( NH )
      CHARACTER *(2) ORD
      DOUBLE PRECISION RVEC( NDIM ), OLDSV( NDIM ),
     1                 RI, RO, GAUX( NTHPTS ), FTF1( 2*NPHPTS ),
     2                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ), FAC,
     3                 VF1( NPHPTS, NTHPTS, 3), GAUW( NTHPTS ),
     4                 VF2( NPHPTS, NTHPTS, 3)
      DOUBLE PRECISION QST( LH*(LH + 2), 3 ),
     1                 PA( (LH+1)*(LH+2)/2 , NTHPTS ),
     2                 DPA( (LH+1)*(LH+2)/2 , NTHPTS ),
     2                 SF( NTHPTS, NPHPTS ),
     4                 SHC( LH*(LH+2) ), DSHC( LH*(LH+2) )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRN, IH, L, M, ICS, NOHARM, INDSHC, IND, IPTF, IVLMF
      DOUBLE PRECISION ZCOEF, TOL
      PARAMETER ( TOL = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of NDIM against NH and NR etc.
C
      IF ( NDIM.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPVVTA, bad array size'
        PRINT *,' NDIM = ',NDIM
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FAC ).LT.TOL ) RETURN
C
C Begin looping around radial grid nodes ...
C
      DO IRN = 1, NR
C       .
C       . First, evaluate grad(theta) over sphere at this point
C       .
        CALL SVGRTH ( LH, NDIM, NR, NH, NPHPTS, NTHPTS, MHT, MHL,
     1                MHM, MHC, IRN, GAUX, PA, DPA, RI, RO,
     2                OLDSV, ORD, VF1, SHC, DSHC, FTF1, FTF2,
     3                FTF3 )
C
C       grad( theta ) evaluated at node IRN is now stored
C       in the 3-d array VF1
C
        IPTF = 1
        IVLMF = 1
        CALL SV2QST ( LH, NDIM, NR, NH, MHT, MHL, MHM, MHC, IRN,
     1                IVLMF, IPTF, ORD, OLDSV, RI, RO, QST )
C
        CALL QST2VA ( QST, VF2, GAUX, PA, DPA, FTF1, FTF2, FTF3,
     1                LH, LH, NTHPTS, NTHPTS, NPHPTS, NPHPTS )
C
C The velocity at node IRN is now stored in VF2
C Now take dot product and put in scalar function SF
C
        CALL VFDP ( VF1, VF2, SF, NPHPTS, NTHPTS )
C
C let's transform SF back into spherical harmonic components
C
        CALL FORSSA ( SHC, SF, GAUW, PA, FTF1, LH,
     1            LH, NTHPTS, NTHPTS, NPHPTS, NPHPTS, ZCOEF )
C
C now add the transformed coeff.s at node IRN to RHS vector
C
C Loop around harmonics and select those with
C temperature allocation ...
C
        DO IH = 1, NH
           IF ( MHT( IH ).EQ.3 ) THEN
             L   = MHL( IH )
             M   = MHM( IH )
             ICS = MHC( IH )
             IND =  ( IRN - 1 )*NH + IH
             IF ( L.EQ.0 .AND. M.EQ.0 .AND. ICS.EQ.1 ) THEN
               RVEC( IND ) = RVEC( IND ) + FAC*ZCOEF
             ELSE
               NOHARM = INDSHC( L, M, ICS )
               RVEC( IND ) = RVEC( IND ) + FAC*SHC( NOHARM )
             ENDIF
           ENDIF
        ENDDO
C
      ENDDO
C
C End looping around radial grid nodes ...
C
      RETURN
      END
C*********************************************************************
