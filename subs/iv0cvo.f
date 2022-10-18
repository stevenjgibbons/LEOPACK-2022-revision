C*********************************************************************
C subroutine Identical node Velocity_0 cross Curl Velocity add (Opt) *
C            -              -        -       -    -             -    *
C Steve Gibbons Tue Oct  9 10:05:16 WEST 2001                        C
C____________________________________________________________________C
C                                                                    C
C Subtracts the v_0 x curl v  terms from the matrix for solving      C
C the vorticity equation. The v_0 terms are contained in a vector    C
C VEC0 which has the same grid spacing as that in the matrix         C
C (NR grid points defined by the array XARR).                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NR        : Number of radial grid nodes.                       C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
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
C     ILNRP     : First radial node to act on.                       C
C                  when contributing to the poloidal vorticity eqn.s C
C     IRNRP     : Last radial node to act on.                        C
C                  when contributing to the poloidal vorticity eqn.s C
C     ILNRT     : First radial node to act on.                       C
C                  when contributing to the toroidal vorticity eqn.s C
C     IRNRT     : Last radial node to act on.                        C
C                  when contributing to the toroidal vorticity eqn.s C
C                                                                    C
C     INI       : Int. parameter array corresponding to VI.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INI( 1 ) = IFORMF                                  C
C                 INI( 2 ) = NRR      See INDFUN for details         C
C                 INI( 3 ) = NH                                      C
C                                                                    C
C     MTI      : Array length ( * ) - atleast length NHI             C
C                  See above for key. (corresponds to input vec.)    C
C     MLI      : Array length ( * ) - atleast length NHI             C
C                  Sph. harm. degree, l.                             C
C     MMI      : Array length ( * ) - atleast length NHI             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MPI      : Array length ( * ) - atleast length NHI             C
C                  Pointer array to finite difference coefficients.  C
C                  MPI( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     MTO      : Array length ( * ) - atleast length NHI             C
C                  See above for key. (corresponds to output vec.)   C
C     MLO      : Array length ( * ) - atleast length NHI             C
C                  Sph. harm. degree, l.                             C
C     MMO      : Array length ( * ) - atleast length NHI             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     IN0       : Int. parameter array corresponding to VEC0.        C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 IN0( 1 ) = IFORM0                                  C
C                 IN0( 2 ) = NR0      See INDFUN for details         C
C                 IN0( 3 ) = NH0                                     C
C                                                                    C
C     MT0      : Array length ( * ) - atleast length NH0             C
C                  See above for key. (corresponds to VEC0).         C
C     ML0      : Array length ( * ) - atleast length NH0             C
C                  Sph. harm. degree, l.                             C
C     MM0      : Array length ( * ) - atleast length NH0             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MP0      : Array length ( * ) - atleast length NH0             C
C                  Pointer array to finite difference coefficients.  C
C                  MP0( ih ) = is, which is the 4th index of         C
C                  array SVFDC0 - indicates f.d. scheme used.        C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     NDRVS     : Highest derivative stored in SVFDC.                C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     NBN0      : Number of nodes on each side of point for          C
C                  central differences in VEC0.                      C
C                                                                    C
C     NDCS0     : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC0.                       C
C     NDRVS0    : Highest derivative stored in SVFDC0.               C
C     NDRVM0    : Limit on NDRVS0. Array bound for SVFDC0.           C
C                                                                    C
C     NFDCM0    : Leading dim of SVFDC0. See SVFDCF0.                C
C                  (Must be atleast 2*NBN0 + 1 )                     C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NTHP      : The number of theta points.                        C
C     NPHP      : The number of phi points.                          C
C     MMAX      : Maximum sph. harmonic order, m.                    C
C     M0        : Minimum non-zero wavenumber.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 )                             C
C     FAC       : Multiplier of v0. grad Theta to be added.          C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     VEC0      : Dim ( * ). Vector containing v0.                   C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM , NR, NDRVM +1, NDCS  ).        C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C     SVFDC0    : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM0, NR, NDRVM0+1, NDCS0 ).        C  
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     FF1       : Work array - dim. (2*NPHP)                         C
C     FF2       : Work array - dim. (2*NPHP)                         C
C     FF3       : Work array - dim. (2*NPHP)                         C
C     VF1       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF2       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF3       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     QST       : Work array - dim. ( LH*(LH+2) , 3)                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE IV0CVO( NR, N1, N2, KL, KU, KLE, IMF, ILNRP, IRNRP,
     1    ILNRT, IRNRT, INI, MTI, MLI, MMI, MPI, MTO, MLO, MMO, IN0,
     2    MT0, ML0, MM0, MP0, NBN, NDCS, NDRVS, NDRVM, NFDCM, NBN0,
     3    NDCS0, NDRVS0, NDRVM0, NFDCM0, LH, NTHP, NPHP, MMAX, A, FAC,
     4    XARR, VEC0, SVFDC, SVFDC0, GAUX, GAUW, PA, DPA, FF1, FF2,
     5    FF3, VF1, VF2, VF3, QST, M0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, N1, N2, KL, KU, KLE, IMF, ILNRP, IRNRP, ILNRT,
     1        IRNRT, INI( * ), MTI( * ), MLI( * ), MMI( * ), MPI( * ),
     2        MTO( * ), MLO( * ), MMO( * ), IN0( * ), MT0( * ),
     3        ML0( * ), MM0( * ), MP0( * ), M0
      INTEGER NBN, NDCS, NDRVS, NDRVM, NFDCM, NBN0, NDCS0, NDRVS0,
     1        NDRVM0, NFDCM0, LH, NTHP, NPHP, MMAX
      DOUBLE PRECISION A( N1, N2 ), FAC, XARR( NR ), VEC0( * ),
     1                 SVFDC( NFDCM , NR, NDRVM +1, NDCS  ),
     2                 SVFDC0( NFDCM0, NR, NDRVM0+1, NDCS0 ),
     3                 GAUX( NTHP ), GAUW( NTHP )
      DOUBLE PRECISION PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     1                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     2                 FF1( 2*NPHP ), FF2( 2*NPHP ), FF3( 2*NPHP )
      DOUBLE PRECISION VF1( NPHP, NTHP, 3), VF2( NPHP, NTHP, 3),
     1                 VF3( NPHP, NTHP, 3 ), QST( LH*( LH + 2), 3 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IHA, IHB, IHG, LA, LB, LG, MA, MB, MG, ICSA, ICSB,
     1        ICSG, IQSTA, IQSTB, IQSTG, IOLDF, IHD, IHD0,
     2        IPARS( 7 ), NHA, NHB, NHG, INDG, INDSHC, IS, IS0
C     .
C     . ioldf = 1 --> existing alpha function, new beta function
C     .
      PARAMETER ( IOLDF = 1 )
C     .
      DOUBLE PRECISION DPARS( 1 ), ZCOEF, LOW, WORK( 6 ),
     1                 CQTS, CQTT, CQSS, CQST, CSTQ, CSQS,
     2                 CSQT, CSSQ, CTTQ, CTQS, CTQT, CTSQ
      PARAMETER ( LOW = 1.0d-9 )
      EXTERNAL INVCVT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NHA = IN0( 3 )
      NHB = INI( 3 )
      NHG = INI( 3 )
C     .
C     . Early exit with zero value of multiplier
C     .
      IF ( DABS( FAC ).LT.LOW ) RETURN
C     .
C     . IHA is the 'alpha' harmonic.
C     . This is the velocity.
C     . In this case the velocity is in the VEC0
C     . vector and so loop around IHA from 1 to NH0
C     .
      DO IHA = 1, NHA
C       .
C       . Consider the case of alpha poloidal
C       .
        IF ( MT0( IHA ).EQ.1 ) THEN
          LA = ML0( IHA )
          IF ( MM0( IHA ).LT.0 ) THEN
            MA   = -MM0( IHA )
            ICSA = 2
          ELSE
            MA   = MM0( IHA )
            ICSA = 1
          ENDIF
C         .
          IS0    = MP0( IHA )
C         .
C         . Now evaluate q_{alpha} in VF1
C         .
          IQSTA = 1
C         .
          CALL SHVECO( LA, MA, ICSA, IQSTA, VF1, GAUX, PA, DPA,
     1                 NTHP, NPHP, LH, M0 )
C         .
C         . Now loop IHB around harmonics in solution vector
C         .
          DO IHB = 1, NHB
C           .
C           . Case of beta poloidal
C           .
            IF ( MTI( IHB ).EQ.1 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLI( IHB )
              IF ( MMI( IHB ).LT.0 ) THEN
                MB   = -MMI( IHB )
                ICSB = 2
              ELSE
                MB   = MMI( IHB )
                ICSB = 1
              ENDIF
C             .
              IS     = MPI( IHB )
C             .
C             . Evaluate t_{beta} in VF2
C             .
              IQSTB = 3
C             .
              CALL SHVECO( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                     PA, DPA, NTHP, NPHP, LH, M0 )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QSO( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX, M0 )
C             .
C             . QST now contains C_{qts}^{abg} in comp. 2
C             . QST now contains C_{qtt}^{abg} in comp. 3
C             .
C             . Loop around output harmonics
C             .
              DO IHG = 1, NHG
                LG = MLO( IHG )
                IF ( MMO( IHG ).LT.0 ) THEN
                  MG   = -MMO( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMO( IHG )
                  ICSG = 1
                ENDIF
C               .
                INDG   = INDSHC( LG, MG, ICSG )
C               .
                IF ( INDG.GT.0 ) CQTS   = QST( INDG, 2 )
                IF ( INDG.GT.0 ) CQTT   = QST( INDG, 3 )
C               .
C               . Case of Gamma toroidal
C               . See if we have C_{qts}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 .AND. DABS( CQTS ).GT.LOW ) THEN
                  IQSTG = 2
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 3
                  IHD0       = 1
C                 .
                  DPARS( 1 ) = CQTS
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{qtt}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.1 .AND. DABS( CQTT ).GT.LOW ) THEN
                  IQSTG = 3
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 2
                  IHD0       = 0
C                 .
                  DPARS( 1 ) = CQTT
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRP, IRNRP, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
            ENDIF
C           .
C           .
C           . Case of beta toroidal
C           .
            IF ( MTI( IHB ).EQ.2 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLI( IHB )
              IF ( MMI( IHB ).LT.0 ) THEN
                MB   = -MMI( IHB )
                ICSB = 2
              ELSE
                MB   = MMI( IHB )
                ICSB = 1
              ENDIF
C             .
              IS     = MPI( IHB )
C             .
C             . Evaluate s_{beta} in VF2
C             . (No need to evaluate q_{beta}
C             . - there is no Q x Q )
C             .
              IQSTB = 2
C             .
              CALL SHVECO( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                     PA, DPA, NTHP, NPHP, LH, M0 )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QSO( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX, M0 )
C             .
C             . QST now contains C_{qss}^{abg} in comp. 2
C             . QST now contains C_{qst}^{abg} in comp. 3
C             .
C             . Loop around output harmonics
C             .
              DO IHG = 1, NHG
                LG = MLO( IHG )
                IF ( MMO( IHG ).LT.0 ) THEN
                  MG   = -MMO( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMO( IHG )
                  ICSG = 1
                ENDIF
C               .
                INDG   = INDSHC( LG, MG, ICSG )
C               .
                IF ( INDG.GT.0 ) CQSS   = QST( INDG, 2 )
                IF ( INDG.GT.0 ) CQST   = QST( INDG, 3 )
C               .
C               . Case of Gamma toroidal
C               . See if we have C_{qss}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 .AND. DABS( CQSS ).GT.LOW ) THEN
                  IQSTG = 2
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 2
                  IHD0       = 1
C                 .
                  DPARS( 1 ) = CQSS
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{qst}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.1 .AND. DABS( CQST ).GT.LOW ) THEN
                  IQSTG = 3
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 1
                  IHD0       = 0
C                 .
                  DPARS( 1 ) = CQST
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRP, IRNRP, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
            ENDIF
C           .
          ENDDO
C         (end ihb = 1, nhb)
C         .
C         . Now evaluate s_{alpha} in VF1
C         .
          IQSTA = 2
C         .
          CALL SHVECO( LA, MA, ICSA, IQSTA, VF1, GAUX, PA, DPA,
     1                 NTHP, NPHP, LH, M0 )
C         .
C         . Now loop IHB around harmonics in solution vector
C         .
          DO IHB = 1, NHB
C           .
C           . Case of beta poloidal
C           .
            IF ( MTI( IHB ).EQ.1 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLI( IHB )
              IF ( MMI( IHB ).LT.0 ) THEN
                MB   = -MMI( IHB )
                ICSB = 2
              ELSE
                MB   = MMI( IHB )
                ICSB = 1
              ENDIF
C             .
              IS     = MPI( IHB )
C             .
C             . Evaluate t_{beta} in VF2
C             .
              IQSTB = 3
C             .
              CALL SHVECO( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                     PA, DPA, NTHP, NPHP, LH, M0 )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QSO( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX, M0 )
C             .
C             . QST now contains C_{stq}^{abg} in comp. 1
C             .
C             . Loop around output harmonics
C             .
              DO IHG = 1, NHG
                LG = MLO( IHG )
                IF ( MMO( IHG ).LT.0 ) THEN
                  MG   = -MMO( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMO( IHG )
                  ICSG = 1
                ENDIF
C               .
                INDG   = INDSHC( LG, MG, ICSG )
C               .
                IF ( INDG.GT.0 ) CSTQ   = QST( INDG, 1 )
C               .
C               . Case of Gamma toroidal
C               . See if we have C_{stq}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 .AND. DABS( CSTQ ).GT.LOW ) THEN
                  IQSTG = 1
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 2
                  IHD0       = 1
C                 .
                  DPARS( 1 ) = CSTQ
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
            ENDIF
C           .
C           . Case of beta toroidal
C           .
            IF ( MTI( IHB ).EQ.2 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLI( IHB )
              IF ( MMI( IHB ).LT.0 ) THEN
                MB   = -MMI( IHB )
                ICSB = 2
              ELSE
                MB   = MMI( IHB )
                ICSB = 1
              ENDIF
C             .
              IS     = MPI( IHB )
C             .
C             . Evaluate q_{beta} in VF2
C             .
              IQSTB = 1
C             .
              CALL SHVECO( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                     PA, DPA, NTHP, NPHP, LH, M0 )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QSO( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX, M0 )
C             .
C             . QST now contains C_{sqs}^{abg} in comp. 2
C             . QST now contains C_{sqt}^{abg} in comp. 3
C             .
C             . Loop around output harmonics
C             .
              DO IHG = 1, NHG
                LG = MLO( IHG )
                IF ( MMO( IHG ).LT.0 ) THEN
                  MG   = -MMO( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMO( IHG )
                  ICSG = 1
                ENDIF
C               .
                INDG   = INDSHC( LG, MG, ICSG )
C               .
                IF ( INDG.GT.0 ) CSQS   = QST( INDG, 2 )
                IF ( INDG.GT.0 ) CSQT   = QST( INDG, 3 )
C               .
C               . Case of Gamma toroidal
C               . See if we have C_{sqs}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 .AND. DABS( CSQS ).GT.LOW ) THEN
                  IQSTG = 2
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 1
                  IHD0       = 2
C                 .
                  DPARS( 1 ) = CSQS
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{sqt}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.1 .AND. DABS( CSQT ).GT.LOW ) THEN
                  IQSTG = 3
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 0
                  IHD0       = 1
C                 .
                  DPARS( 1 ) = CSQT
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRP, IRNRP, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
C             . Evaluate s_{beta} in VF2
C             .
              IQSTB = 2
C             .
              CALL SHVECO( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                     PA, DPA, NTHP, NPHP, LH, M0 )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QSO( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX, M0 )
C             .
C             . QST now contains C_{ssq}^{abg} in comp. 1
C             .
C             . Loop around output harmonics
C             .
              DO IHG = 1, NHG
                LG = MLO( IHG )
                IF ( MMO( IHG ).LT.0 ) THEN
                  MG   = -MMO( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMO( IHG )
                  ICSG = 1
                ENDIF
C               .
                INDG   = INDSHC( LG, MG, ICSG )
C               .
                IF ( INDG.GT.0 ) CSSQ   = QST( INDG, 1 )
C               .
C               . Case of Gamma toroidal
C               . See if we have C_{ssq}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 .AND. DABS( CSSQ ).GT.LOW ) THEN
                  IQSTG = 1
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 1
                  IHD0       = 1
C                 .
                  DPARS( 1 ) = CSSQ
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
            ENDIF
C           .
          ENDDO
C         (end loop ihb = 1, nhb)
C         .
        ENDIF
C       .
C       . Consider the case of alpha toroidal
C       .
        IF ( MT0( IHA ).EQ.2 ) THEN
          LA = ML0( IHA )
          IF ( MM0( IHA ).LT.0 ) THEN
            MA   = -MM0( IHA )
            ICSA = 2
          ELSE
            MA   = MM0( IHA )
            ICSA = 1
          ENDIF
C         .
          IS0    = MP0( IHA )
C         .
C         . Now evaluate t_{alpha} in VF1
C         .
          IQSTA = 3
C         .
          CALL SHVECO( LA, MA, ICSA, IQSTA, VF1, GAUX, PA, DPA,
     1                 NTHP, NPHP, LH, M0 )
C         .
C         . Now loop IHB around harmonics in solution vector
C         .
          DO IHB = 1, NHB
C           .
C           . Case of beta poloidal
C           .
            IF ( MTI( IHB ).EQ.1 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLI( IHB )
              IF ( MMI( IHB ).LT.0 ) THEN
                MB   = -MMI( IHB )
                ICSB = 2
              ELSE
                MB   = MMI( IHB )
                ICSB = 1
              ENDIF
C             .
              IS     = MPI( IHB )
C             .
C             . Evaluate t_{beta} in VF2
C             .
              IQSTB = 3
C             .
              CALL SHVECO( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                     PA, DPA, NTHP, NPHP, LH, M0 )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QSO( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX, M0 )
C             .
C             . QST now contains C_{ttq}^{abg} in comp. 1
C             .
C             . Loop around output harmonics
C             .
              DO IHG = 1, NHG
                LG = MLO( IHG )
                IF ( MMO( IHG ).LT.0 ) THEN
                  MG   = -MMO( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMO( IHG )
                  ICSG = 1
                ENDIF
C               .
                INDG   = INDSHC( LG, MG, ICSG )
C               .
                IF ( INDG.GT.0 ) CTTQ   = QST( INDG, 1 )
C               .
C               . Case of Gamma toroidal
C               . See if we have C_{ttq}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 .AND. DABS( CTTQ ).GT.LOW ) THEN
                  IQSTG = 1
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 2
                  IHD0       = 0
C                 .
                  DPARS( 1 ) = CTTQ
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
            ENDIF
C           .
C           .
C           . Case of beta toroidal
C           .
            IF ( MTI( IHB ).EQ.2 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLI( IHB )
              IF ( MMI( IHB ).LT.0 ) THEN
                MB   = -MMI( IHB )
                ICSB = 2
              ELSE
                MB   = MMI( IHB )
                ICSB = 1
              ENDIF
C             .
              IS     = MPI( IHB )
C             .
C             . Evaluate q_{beta} in VF2
C             .
              IQSTB = 1
C             .
              CALL SHVECO( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                     PA, DPA, NTHP, NPHP, LH, M0 )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QSO( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX, M0 )
C             .
C             . QST now contains C_{tqs}^{abg} in comp. 2
C             . QST now contains C_{tqt}^{abg} in comp. 3
C             .
C             . Loop around output harmonics
C             .
              DO IHG = 1, NHG
                LG = MLO( IHG )
                IF ( MMO( IHG ).LT.0 ) THEN
                  MG   = -MMO( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMO( IHG )
                  ICSG = 1
                ENDIF
C               .
                INDG   = INDSHC( LG, MG, ICSG )
C               .
                IF ( INDG.GT.0 ) CTQS   = QST( INDG, 2 )
                IF ( INDG.GT.0 ) CTQT   = QST( INDG, 3 )
C               .
C               . Case of Gamma toroidal
C               . See if we have C_{tqs}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 .AND. DABS( CTQS ).GT.LOW ) THEN
                  IQSTG = 2
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 1
                  IHD0       = 1
C                 .
                  DPARS( 1 ) = CTQS
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{tqt}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.1 .AND. DABS( CTQT ).GT.LOW ) THEN
                  IQSTG = 3
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 0
                  IHD0       = 0
C                 .
                  DPARS( 1 ) = CTQT
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRP, IRNRP, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
C             .
C             . Evaluate s_{beta} in VF2
C             .
              IQSTB = 2
C             .
              CALL SHVECO( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                     PA, DPA, NTHP, NPHP, LH, M0 )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QSO( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX, M0 )
C             .
C             . QST now contains C_{tsq}^{abg} in comp. 1
C             .
C             . Loop around output harmonics
C             .
              DO IHG = 1, NHG
                LG = MLO( IHG )
                IF ( MMO( IHG ).LT.0 ) THEN
                  MG   = -MMO( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMO( IHG )
                  ICSG = 1
                ENDIF
C               .
                INDG   = INDSHC( LG, MG, ICSG )
C               .
                IF ( INDG.GT.0 ) CTSQ   = QST( INDG, 1 )
C               .
C               . Case of Gamma toroidal
C               . See if we have C_{tsq}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 .AND. DABS( CTSQ ).GT.LOW ) THEN
                  IQSTG = 1
C                 .
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 1
                  IHD0       = 0
C                 .
                  DPARS( 1 ) = CTSQ
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHA, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
            ENDIF
C           .
          ENDDO
C         (end ihb = 1, nhb)
C         .
        ENDIF
C       .
      ENDDO
C     (end loop iha = 1, nha)
      RETURN
      END
C*********************************************************************

