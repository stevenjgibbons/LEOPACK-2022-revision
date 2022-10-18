C*********************************************************************
C subroutine Double Precision Vector Velocity dot grad ( Velocity )  *
C            -      -         -      -                   -           *
C Steve Gibbons 26.6.99                                              C
C____________________________________________________________________C
C Adds the term to the double precision matrix resulting from        C
C curl [ vel_0 . Grad ( vel_0 ) ] in the momentum equation.          C
C In practice, this routine SUBTRACTS curl[ u \times curl [ u ] ]    C
C from RVEC.                ^^^^^^^^^                                C
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
C     SV1       : Workspace of dimension ( NDIM )                    C
C     SV2       : Workspace of dimension ( NDIM )                    C
C     SV3       : Workspace of dimension ( NDIM )                    C
C     QST       : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector                           C
C                  Has dimensions (  LH*(LH+2) , 3).                 C
C     VF1       : Vector Function. An array of dimensions            C
C                  ( NPHPTS, NTHPTS, 3) which contain the R, THETA   C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IPHI, ITHETA, 2 ) is the Theta      C
C                  compontent of the vector at (iphi, itheta).       C
C     VF2       : Second vector Function. See above                  C
C     VF3       : Third vector Function. See above                   C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPVVVA ( LH, NDIM, NTHPTS, NPHPTS, NR, NH, MHT, MHL,
     1                    MHM, MHC, ORD, RVEC, OLDSV, RI, RO, FAC, 
     2                    SV1, SV2, SV3, QST, VF1, VF2, VF3, FTF1,
     3                    FTF2, FTF3, GAUX, GAUW, PA, DPA )
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
     4                 VF2( NPHPTS, NTHPTS, 3),
     5                 VF3( NPHPTS, NTHPTS, 3)
      DOUBLE PRECISION QST( LH*(LH + 2), 3 ),
     1                 PA( (LH+1)*(LH+2)/2 , NTHPTS ),
     2                 DPA( (LH+1)*(LH+2)/2 , NTHPTS )
      DOUBLE PRECISION SV1( NDIM ), SV2( NDIM ), SV3( NDIM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRN, IH, L, M, ICS, NOHARM, INDSHC, IND, IPTF, IVLMF
      DOUBLE PRECISION RAD, TOL, TEMFAC, D0F, D1F, D2F, D3F, D4F, H,
     1                 SQRLL1
      PARAMETER ( TOL = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      H = (RO - RI)/DBLE( NR - 1 )
C
C Check value of NDIM against NH and NR etc.
C
      IF ( NDIM.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPVVVA, bad array size'
        PRINT *,' NDIM = ',NDIM
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FAC ).LT.TOL ) RETURN
C
C Need to put the curl of OLDSV into SV1
C Put the curl of the toroidal harmonics of OLDSV
C into the toroidal parts of SV1 - note however
C that here, they actually represent Poloidal harmonics
C Put the curl of the poloidal harms of OLDSV into the
C poloidal parts of SV1. Here though, they represent
C toroidal harmonics.
C
      DO IH = 1, NH
        IF ( MHT( IH ).EQ.1 ) THEN
          L = MHL( IH )
          DO IRN = 1, NR
            IND = ( IRN - 1 )*NH + IH
            RAD = RI + H*DBLE( IRN - 1 )
            CALL SVDERF ( NDIM, NH, IH, NR, IRN, OLDSV, ORD,
     1                    D0F, D1F, D2F, D3F, D4F, RI, RO )
            TEMFAC = DBLE( L*L + L )*D0F/(RAD*RAD)
            TEMFAC = TEMFAC - D2F - 2.0d0*D1F/RAD
            SV1( IND ) = TEMFAC
          ENDDO
        ENDIF
        IF ( MHT( IH ).EQ.2 ) THEN
          DO IRN = 1, NR
            IND = ( IRN - 1 )*NH + IH
            SV1( IND ) = OLDSV( IND )
          ENDDO
        ENDIF
      ENDDO
C
C We now have the poloidal/toroidal decompositions of u
C and curl ( u ) in OLDSV and SV1 respectively.
C
C
C Begin looping around radial grid nodes ...
C
      DO IRN = 1, NR
c        .
c        . Find qst coeff.s for the velocity at this node
c        .
         IPTF = 1
         IVLMF = 1
         CALL SV2QST ( LH, NDIM, NR, NH, MHT, MHL, MHM, MHC, IRN,
     1                 IVLMF, IPTF, ORD, OLDSV, RI, RO, QST )
c        .
c        . Evaluate the vector in space
c        .
         CALL QST2VA ( QST, VF1, GAUX, PA, DPA, FTF1, FTF2, FTF3,
     1                 LH, LH, NTHPTS, NTHPTS, NPHPTS, NPHPTS )
c        .
c        . Find qst coeff.s for ( curl [u] ) at this node
c        .
         IPTF = 2
         IVLMF = 1
         CALL SV2QST ( LH, NDIM, NR, NH, MHT, MHL, MHM, MHC, IRN,
     1                 IVLMF, IPTF, ORD, SV1, RI, RO, QST )
c        .
c        . Evaluate the vector in space
c        .
         CALL QST2VA ( QST, VF2, GAUX, PA, DPA, FTF1, FTF2, FTF3,
     1                 LH, LH, NTHPTS, NTHPTS, NPHPTS, NPHPTS )
c        .
c        . VF1 now contains [u] in space and 
c        . VF2 now contains curl [u] in space
c        . so let's take the cross product, the
c        . result going into VF3.
c        .
         CALL VFCP ( VF1, VF2, VF3, NPHPTS, NTHPTS )
c        .
c        . Now we reconvert this vector back to qst values
c        .
         CALL VF2QSA ( QST, VF3, GAUX, GAUW, PA, DPA, FTF1,
     1                 FTF2, FTF3, LH, LH, NTHPTS, NTHPTS,
     2                 NPHPTS, NPHPTS )
c        .
c        . So now, at node IRN, we have the QST values for
c        . u x curl u ... these are the same harmonics as
c        . the initial velocity and so we will put
c        . the Q harmonics into the poloidal part of SV2,
c        . the S harmonics into the poloidal part of SV3 and
c        . the T harmonics into the toroidal part of SV2.
c        .
         DO IH = 1, NH
            L   = MHL( IH )
            M   = MHM( IH )
            ICS = MHC( IH )
            NOHARM = INDSHC( L, M, ICS )
            IND = ( IRN - 1 )*NH + IH
            IF ( MHT( IH ).EQ.1 ) THEN
              SV2( IND ) = QST( NOHARM, 1 )
              SV3( IND ) = QST( NOHARM, 2 )
            ENDIF
            IF ( MHT( IH ).EQ.2 ) THEN
              SV2( IND ) = QST( NOHARM, 3 )
            ENDIF
         ENDDO
c        .
      ENDDO
C
C End looping around radial grid nodes ...
C
C We now need to take the curl of the qst vectors
C stored in SV2 and SV3, and add to RVEC ...
C
      DO IRN = 1, NR
        RAD = RI + H*DBLE( IRN - 1 )
        DO IH = 1, NH
           L = MHL( IH )
           IND = ( IRN - 1 )*NH + IH
           IF ( MHT( IH ).EQ.1 ) THEN
             RVEC( IND ) = RVEC( IND ) - FAC*SV2( IND )/RAD
             CALL SVDERF ( NDIM, NH, IH, NR, IRN, SV3, ORD,
     1                     D0F, D1F, D2F, D3F, D4F, RI, RO )
             RVEC( IND ) = RVEC( IND ) + FAC*(D1F+D0F/RAD)/SQRLL1(L)
           ENDIF
           IF ( MHT( IH ).EQ.2 ) THEN
             RVEC( IND ) = RVEC( IND ) + FAC*SV2( IND )/SQRLL1(L)
           ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
