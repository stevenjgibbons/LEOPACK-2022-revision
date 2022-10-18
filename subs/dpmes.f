C*********************************************************************
C subroutine Double Precision Matrix Entry Subroutine ****************
C            -      -         -      -     -          ****************
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C Enters contributions from one harmonic to another into the         C
C convection/dynamo matrix.                                          C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of harmonics (all types)                    C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C                                                                    C
C     IMF       : Matrix format flag.                                C
C                                                                    C
C          imf = 1; Matrix is in LAPACK banded format                C
C                   ie element a_{i,j} is stored in                  C
C                   A( kle + ku + 1 + i - j , j )                    C
C                                                                    C
C          imf = 2; Matrix is banded but with element a_{i,j}        C
C                   stored in A( kl + 1 + j - i , i ).               C
C                                                                    C
C          imf = 3; Matrix is square - ie a_{i,j} is stored          C
C                   in A( i, j ).                                    C
C                                                                    C
C      IIH      : Number of the input (column) harmonic              C
C      IOH      : Number of the output (row) harmonic                C
C                                                                    C
C     INTPAR    : Integer parameters required for subroutine SUB1.   C
C                   ( Dimension (*) )                                C
C  Subroutines                                                       C
C  -----------                                                       C
C     SUB1      : Determines what multiplies each derivative in      C
C                  the matrix. Must have calling sequence ...        C
C                                                                    C
C     CALL SUB1 ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD, RI, RO,    C
C                 INTPAR, DPPARS, VEC, ORD, IRN )                    C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
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
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     AMAT      : Matrix. Dimensions ( N1, N2 )                      C
C              Will generally be banded due to the nature of the     C
C              numerical scheme. KL, KU and KLE parameterise this.   C
C                                                                    C
C     DPPARS    : Arbitrary (*) array.                               C
C     VEC       : Arbitrary (*) array.                               C
C     FAC       : Multiplying factor.                                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH, 
     1                   IOH, INTPAR, SUB1, ORD, RI, RO, AMAT,
     2                   DPPARS, VEC, FAC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, NR, NH, KL, KU, KLE, IMF, IIH, IOH, INTPAR( * )
      EXTERNAL SUB1
      CHARACTER *(2) ORD
      DOUBLE PRECISION AMAT( N1, N2 ), RI, RO, DPPARS( * ), VEC( * ),
     1                 FAC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION FM4, FM3, FM2, FM1, F00, FP1, FP2, FP3, FP4,
     1                 CM4( 4 ), CM3( 4 ), CM2( 4 ), CM1( 4 ),
     2                 C00( 4 ), CP1( 4 ), CP2( 4 ), CP3( 4 ),
     3                 CP4( 4 ), H
      DOUBLE PRECISION RAD, D0FAC, D1FAC, D2FAC, D3FAC, D4FAC,
     1                 TOL
      INTEGER IRN, IR, ICOL, IROW, IM4, IM3, IM2, IM1, I00,
     1        IP1, IP2, IP3, IP4
      LOGICAL LB
      PARAMETER ( TOL = 1.0d-8 )
C The flag LB stands for Long Band width - this is set to true
C in the cases of ORD = 'SF' and 'O7' and false otherwise
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of N2 against NH and NR etc.
C
      IF ( N2.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPMES, bad array size'
        PRINT *,' N2 (second array dimension) = ',N2
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FAC ).LT.TOL ) RETURN
C
      IF ( ORD.EQ.'SF' .OR. ORD.EQ.'O7' ) THEN
         LB = .TRUE.
      ELSE
         LB = .FALSE.
      ENDIF
C
      H = ( RO - RI ) / DBLE( NR - 1 )
C     ........................ consider IRN = 1 ...................
      IRN = 1
      RAD = RI
      CALL GSLDCF ( H, C00, CP1, CP2, CP3, CP4 )
      IR  = ( IRN - 1 )*NH + IOH
      I00 = ( IRN - 1 )*NH + IIH
      IP1 = ( IRN     )*NH + IIH
      IP2 = ( IRN + 1 )*NH + IIH
      IP3 = ( IRN + 2 )*NH + IIH
      IP4 = ( IRN + 3 )*NH + IIH
      CALL SUB1 ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD, RI, RO,
     1            INTPAR, DPPARS, VEC, ORD, IRN)
      F00 = C00( 1 )*D1FAC + C00( 2 )*D2FAC +
     1      C00( 3 )*D3FAC + C00( 4 )*D4FAC     + D0FAC
      FP1 = CP1( 1 )*D1FAC + CP1( 2 )*D2FAC +
     1      CP1( 3 )*D3FAC + CP1( 4 )*D4FAC
      FP2 = CP2( 1 )*D1FAC + CP2( 2 )*D2FAC +
     1      CP2( 3 )*D3FAC + CP2( 4 )*D4FAC
      FP3 = CP3( 1 )*D1FAC + CP3( 2 )*D2FAC +
     1      CP3( 3 )*D3FAC + CP3( 4 )*D4FAC
      FP4 = CP4( 1 )*D1FAC + CP4( 2 )*D2FAC +
     1      CP4( 3 )*D3FAC + CP4( 4 )*D4FAC
C
      CALL MATIND (IR,I00,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * F00
      CALL MATIND (IR,IP1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP1
      CALL MATIND (IR,IP2,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP2
      CALL MATIND (IR,IP3,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP3
      CALL MATIND (IR,IP4,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP4
c     .
C     ........................ consider IRN = 2 ...................
      IRN = 2
      RAD = RI + H
      CALL GFLDCF ( H, CM1, C00, CP1, CP2, CP3 )
      IR  = ( IRN - 1 )*NH + IOH
      IM1 = ( IRN - 2 )*NH + IIH
      I00 = ( IRN - 1 )*NH + IIH
      IP1 = ( IRN     )*NH + IIH
      IP2 = ( IRN + 1 )*NH + IIH
      IP3 = ( IRN + 2 )*NH + IIH
      CALL SUB1 ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD, RI, RO,
     1            INTPAR, DPPARS, VEC, ORD, IRN)
      FM1 = CM1( 1 )*D1FAC + CM1( 2 )*D2FAC +
     1      CM1( 3 )*D3FAC + CM1( 4 )*D4FAC
      F00 = C00( 1 )*D1FAC + C00( 2 )*D2FAC +
     1      C00( 3 )*D3FAC + C00( 4 )*D4FAC     + D0FAC
      FP1 = CP1( 1 )*D1FAC + CP1( 2 )*D2FAC +
     1      CP1( 3 )*D3FAC + CP1( 4 )*D4FAC
      FP2 = CP2( 1 )*D1FAC + CP2( 2 )*D2FAC +
     1      CP2( 3 )*D3FAC + CP2( 4 )*D4FAC
      FP3 = CP3( 1 )*D1FAC + CP3( 2 )*D2FAC +
     1      CP3( 3 )*D3FAC + CP3( 4 )*D4FAC
C
      CALL MATIND (IR,IM1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM1
      CALL MATIND (IR,I00,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * F00
      CALL MATIND (IR,IP1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP1
      CALL MATIND (IR,IP2,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP2
      CALL MATIND (IR,IP3,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP3
c     .
C     ........................ consider IRN = 3 ...................
      IRN = 3
      RAD = RI + 2.0d0*H
      CALL CENDIF( H, 'O5', CM3, CM2, CM1, C00, CP1, CP2, CP3 )
      IR  = ( IRN - 1 )*NH + IOH
      IM2 = ( IRN - 3 )*NH + IIH
      IM1 = ( IRN - 2 )*NH + IIH
      I00 = ( IRN - 1 )*NH + IIH
      IP1 = ( IRN     )*NH + IIH
      IP2 = ( IRN + 1 )*NH + IIH
      CALL SUB1 ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD, RI, RO,
     1            INTPAR, DPPARS, VEC, ORD, IRN)
      FM2 = CM2( 1 )*D1FAC + CM2( 2 )*D2FAC +
     1      CM2( 3 )*D3FAC + CM2( 4 )*D4FAC
      FM1 = CM1( 1 )*D1FAC + CM1( 2 )*D2FAC +
     1      CM1( 3 )*D3FAC + CM1( 4 )*D4FAC
      F00 = C00( 1 )*D1FAC + C00( 2 )*D2FAC +
     1      C00( 3 )*D3FAC + C00( 4 )*D4FAC     + D0FAC
      FP1 = CP1( 1 )*D1FAC + CP1( 2 )*D2FAC +
     1      CP1( 3 )*D3FAC + CP1( 4 )*D4FAC
      FP2 = CP2( 1 )*D1FAC + CP2( 2 )*D2FAC +
     1      CP2( 3 )*D3FAC + CP2( 4 )*D4FAC
C
      CALL MATIND (IR,IM2,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM2
      CALL MATIND (IR,IM1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM1
      CALL MATIND (IR,I00,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * F00
      CALL MATIND (IR,IP1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP1
      CALL MATIND (IR,IP2,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP2
c     .
C     ................ now loop around IRN = 4 to NR - 3 ........
      DO IRN = 4, NR - 3
         RAD = RI + H*DBLE( IRN - 1 )
         IR = ( IRN - 1 )*NH + IOH
         IF ( LB ) IM3 = ( IRN - 4 )*NH + IIH
         IM2 = ( IRN - 3 )*NH + IIH
         IM1 = ( IRN - 2 )*NH + IIH
         I00 = ( IRN - 1 )*NH + IIH
         IP1 = ( IRN     )*NH + IIH
         IP2 = ( IRN + 1 )*NH + IIH
         IF ( LB ) IP3 = ( IRN + 2 )*NH + IIH
         CALL CENDIF( H, ORD, CM3, CM2, CM1, C00, CP1, CP2, CP3 )
         CALL SUB1 ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD, RI, RO,
     1               INTPAR, DPPARS, VEC, ORD, IRN)
         IF ( LB ) FM3 = CM3( 1 )*D1FAC + CM3( 2 )*D2FAC +
     1         CM3( 3 )*D3FAC + CM3( 4 )*D4FAC
         FM2 = CM2( 1 )*D1FAC + CM2( 2 )*D2FAC +
     1         CM2( 3 )*D3FAC + CM2( 4 )*D4FAC
         FM1 = CM1( 1 )*D1FAC + CM1( 2 )*D2FAC +
     1         CM1( 3 )*D3FAC + CM1( 4 )*D4FAC
         F00 = C00( 1 )*D1FAC + C00( 2 )*D2FAC +
     1         C00( 3 )*D3FAC + C00( 4 )*D4FAC     + D0FAC
         FP1 = CP1( 1 )*D1FAC + CP1( 2 )*D2FAC +
     1         CP1( 3 )*D3FAC + CP1( 4 )*D4FAC
         FP2 = CP2( 1 )*D1FAC + CP2( 2 )*D2FAC +
     1         CP2( 3 )*D3FAC + CP2( 4 )*D4FAC
         IF ( LB ) FP3 = CP3( 1 )*D1FAC + CP3( 2 )*D2FAC +
     1         CP3( 3 )*D3FAC + CP3( 4 )*D4FAC
         IF ( LB ) CALL MATIND (IR,IM3,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
         IF ( LB ) AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM3
         CALL MATIND (IR,IM2,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
         AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM2
         CALL MATIND (IR,IM1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
         AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM1
         CALL MATIND (IR,I00,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
         AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * F00
         CALL MATIND (IR,IP1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
         AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP1
         CALL MATIND (IR,IP2,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
         AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP2
         IF ( LB ) CALL MATIND (IR,IP3,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
         IF ( LB ) AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP3
      ENDDO
C     ................ end loop around IRN = 4 to NR - 3 ........
c     .
C     ........................ consider IRN = NR - 2 ............
      IRN = NR - 2
      RAD = RO - 2.0d0*H
      CALL CENDIF( H, 'O5', CM3, CM2, CM1, C00, CP1, CP2, CP3 )
      IR  = ( IRN - 1 )*NH + IOH
      IM2 = ( IRN - 3 )*NH + IIH
      IM1 = ( IRN - 2 )*NH + IIH
      I00 = ( IRN - 1 )*NH + IIH
      IP1 = ( IRN     )*NH + IIH
      IP2 = ( IRN + 1 )*NH + IIH
      CALL SUB1 ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD, RI, RO,
     1            INTPAR, DPPARS, VEC, ORD, IRN)
      FM2 = CM2( 1 )*D1FAC + CM2( 2 )*D2FAC +
     1      CM2( 3 )*D3FAC + CM2( 4 )*D4FAC
      FM1 = CM1( 1 )*D1FAC + CM1( 2 )*D2FAC +
     1      CM1( 3 )*D3FAC + CM1( 4 )*D4FAC
      F00 = C00( 1 )*D1FAC + C00( 2 )*D2FAC +
     1      C00( 3 )*D3FAC + C00( 4 )*D4FAC     + D0FAC
      FP1 = CP1( 1 )*D1FAC + CP1( 2 )*D2FAC +
     1      CP1( 3 )*D3FAC + CP1( 4 )*D4FAC
      FP2 = CP2( 1 )*D1FAC + CP2( 2 )*D2FAC +
     1      CP2( 3 )*D3FAC + CP2( 4 )*D4FAC
C
      CALL MATIND (IR,IM2,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM2
      CALL MATIND (IR,IM1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM1
      CALL MATIND (IR,I00,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * F00
      CALL MATIND (IR,IP1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP1
      CALL MATIND (IR,IP2,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP2
c     .
C           ................. consider IRN = NR - 1  ............
      IRN = NR - 1
      RAD = RO - H
      CALL GFRDCF ( H, CM3, CM2, CM1, C00, CP1 )
      IR  = ( IRN - 1 )*NH + IOH
      IM3 = ( IRN - 4 )*NH + IIH
      IM2 = ( IRN - 3 )*NH + IIH
      IM1 = ( IRN - 2 )*NH + IIH
      I00 = ( IRN - 1 )*NH + IIH
      IP1 = ( IRN     )*NH + IIH
      CALL SUB1 ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD, RI, RO,
     1            INTPAR, DPPARS, VEC, ORD, IRN)
      FM3 = CM3( 1 )*D1FAC + CM3( 2 )*D2FAC +
     1      CM3( 3 )*D3FAC + CM3( 4 )*D4FAC
      FM2 = CM2( 1 )*D1FAC + CM2( 2 )*D2FAC +
     1      CM2( 3 )*D3FAC + CM2( 4 )*D4FAC
      FM1 = CM1( 1 )*D1FAC + CM1( 2 )*D2FAC +
     1      CM1( 3 )*D3FAC + CM1( 4 )*D4FAC
      F00 = C00( 1 )*D1FAC + C00( 2 )*D2FAC +
     1      C00( 3 )*D3FAC + C00( 4 )*D4FAC     + D0FAC
      FP1 = CP1( 1 )*D1FAC + CP1( 2 )*D2FAC +
     1      CP1( 3 )*D3FAC + CP1( 4 )*D4FAC
C
      CALL MATIND (IR,IM3,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM3
      CALL MATIND (IR,IM2,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM2
      CALL MATIND (IR,IM1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM1
      CALL MATIND (IR,I00,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * F00
      CALL MATIND (IR,IP1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FP1
C     .
C           ................. consider IRN = IR .................
      IRN = NR
      RAD = RO
      CALL GSRDCF ( H, CM4, CM3, CM2, CM1, C00 )
      IR  = ( IRN - 1 )*NH + IOH
      IM4 = ( IRN - 5 )*NH + IIH
      IM3 = ( IRN - 4 )*NH + IIH
      IM2 = ( IRN - 3 )*NH + IIH
      IM1 = ( IRN - 2 )*NH + IIH
      I00 = ( IRN - 1 )*NH + IIH
      CALL SUB1 ( D0FAC, D1FAC, D2FAC, D3FAC, D4FAC, RAD, RI, RO,
     1            INTPAR, DPPARS, VEC, ORD, IRN)
      FM4 = CM4( 1 )*D1FAC + CM4( 2 )*D2FAC +
     1      CM4( 3 )*D3FAC + CM4( 4 )*D4FAC
      FM3 = CM3( 1 )*D1FAC + CM3( 2 )*D2FAC +
     1      CM3( 3 )*D3FAC + CM3( 4 )*D4FAC
      FM2 = CM2( 1 )*D1FAC + CM2( 2 )*D2FAC +
     1      CM2( 3 )*D3FAC + CM2( 4 )*D4FAC
      FM1 = CM1( 1 )*D1FAC + CM1( 2 )*D2FAC +
     1      CM1( 3 )*D3FAC + CM1( 4 )*D4FAC
      F00 = C00( 1 )*D1FAC + C00( 2 )*D2FAC +
     1      C00( 3 )*D3FAC + C00( 4 )*D4FAC     + D0FAC
C
      CALL MATIND (IR,IM4,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM4
      CALL MATIND (IR,IM3,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM3
      CALL MATIND (IR,IM2,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM2
      CALL MATIND (IR,IM1,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * FM1
      CALL MATIND (IR,I00,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + FAC * F00
C
      RETURN
      END
C*********************************************************************

