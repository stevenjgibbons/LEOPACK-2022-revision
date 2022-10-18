C*********************************************************************
C subroutine Solution Vector GRadient of THeta find ******************
C            -        -      --          --         ******************
C Steve Gibbons 25.6.99 (adapted from code of 8.9.97                 C
C____________________________________________________________________C
C For a given radial grid node IRN, SVGRTH calculates the gradient   C
C of the temperature function in space.                              C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Maximum degree of spherical harmonics.             C
C     NDIM      : Length of vector, SV.                              C
C     NR        : Number of radial grid nodes                        C
C     NH        : Number of harmonics (all types)                    C
C     NPHPTS    : The number of phi points.                          C
C     NTHPTS    : The number of theta points.                        C
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
C     IRN       : Number of radial node at which values to be        C
C                                                    taken.          C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2,NTHPTS)       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C     RI        : Radius of inner core.                              C
C     RO        : Radius of outer core.                              C
C     SV        : Solution Vector Dimensions ( NDIM = NR * NH )      C
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
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SVGRVF    : Solution Vector GRadient Vector Function           C
C                 has dimensions ( NPHPTS, NTHPTS , 3 )              C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SHC       : Spherical Harmonic Coefficients.                   C
C                  Has the usual ordering.                           C
C                   SHC( l*l ) = P_l^0                               C
C                   SHC( l*l + 2*m - 1 ) = P_l^mc ( m non zero )     C
C                   SHC( l*l + 2*m ) = P_l^ms ( m non zero )         C
C                  Dimension is ( LH*( LH + 2) )                     C
C     DSHC      : The derivatives of the above at the radial nodes.  C
C                  Dimension is ( LH*( LH + 2) )                     C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Have dimensions ( 2*NPHPTS )                      C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVGRTH ( LH, NDIM, NR, NH, NPHPTS, NTHPTS, MHT, MHL,
     1                    MHM, MHC, IRN, GAUX, PA, DPA, RI, RO,
     2                    SV, ORD, SVGRVF, SHC, DSHC, FTF1, FTF2,
     3                    FTF3 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NDIM, NR, NH, NPHPTS, NTHPTS, MHT( NH ),
     1        MHL( NH ), MHM( NH ), MHC( NH ), IRN
      DOUBLE PRECISION SHC( LH*(LH+2) ), DSHC( LH*(LH+2) ),
     1                 SVGRVF( NPHPTS, NTHPTS, 3 ),
     2                 FTF1( 2*NPHPTS ), FTF2( 2*NPHPTS ),
     3                 FTF3( 2*NPHPTS )
      DOUBLE PRECISION GAUX( NTHPTS ), SV( NDIM ), RI, RO
      DOUBLE PRECISION
     1                 PA( (LH+1)*(LH+2)/2 , NTHPTS ),
     2                 DPA( (LH+1)*(LH+2)/2 , NTHPTS )
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L, M, ICS, ILEN, IOP, IH, NOHARM, INDSHC
      DOUBLE PRECISION ZERO, DZCOEF, D0F, D1F, D2F, D3F, D4F, H, RAD
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      H = ( RO - RI ) / DBLE ( NR - 1 )
      RAD = RI + H*DBLE( IRN - 1 )
c
c Zero the arrays SHC and DSHC
c
      ILEN =  LH*(LH+2)
      IOP = 0
      CALL VECOP ( SHC, ZERO, ILEN, IOP )
      CALL VECOP ( DSHC, ZERO, ILEN, IOP )
      DZCOEF = ZERO
c
c Loop around harmonics and select temperature ones
c
      DO IH = 1, NH
         IF ( MHT(IH).EQ.3 ) THEN
            L   = MHL( IH )
            M   = MHM( IH )
            ICS = MHC( IH )
            CALL SVDERF ( NDIM, NH, IH, NR, IRN, SV, ORD,
     1                    D0F, D1F, D2F, D3F, D4F, RI, RO )
C check for the monopole term which is treated
c differently ...
            IF ( L.EQ.0 .AND. M.EQ.0 .AND. ICS.EQ.1 ) THEN
              DZCOEF = D1F
            ELSE
              NOHARM = INDSHC( L, M, ICS )
              SHC( NOHARM ) = D0F
              DSHC( NOHARM ) = D1F
            ENDIF
         ENDIF
      ENDDO
c
c Now transform
c
      CALL GRINVA ( SHC, DSHC, GAUX, RAD, PA, DPA, FTF1, FTF2,
     1              FTF3, SVGRVF, LH, LH, NTHPTS, NTHPTS,
     2              NPHPTS, NPHPTS, DZCOEF )
c
      RETURN
      END
C*********************************************************************
