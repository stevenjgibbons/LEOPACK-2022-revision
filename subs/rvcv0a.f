C*********************************************************************
C subroutine Recalled coefficient Velocity cross Curl Velocity_0 Add *
C            -                    -        -     -    -        - -   *
C Steve Gibbons Tue Jan 18 09:22:22 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Subtracts the v x curl v_0  terms from the matrix for solving      C
C the vorticity equation. The v_0 terms are contained in a vector    C
C VEC0 which has the same grid spacing as that in the matrix         C
C (NR grid points defined by the array XARR).                        C
C                                                                    C
C Unlike IVCV0A which evaluates the vector interactions in situ,     C
C RVCV0A searches through the arrays IHNALP, IHNBET, IHNGAM and TVHI C
C (as pre-calculated by VCCPCC) with the help of VICEXR for          C
C the non-zero interactions.                                         C
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
C     NVI      : Number of vector interactions.                      C
C     IHNALP   : Number of alpha harmonics. Dim ( * )                C
C     IHNBET   : Number of beta harmonics. Dim ( * )                 C
C     IHNGAM   : Number of gamma harmonics. Dim ( * )                C
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
C     CVI       : Coefficients of vector interaction. Dim ( * )      C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. ( * )        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RVCV0A( NR, N1, N2, KL, KU, KLE, IMF, ILNRP, IRNRP,
     1    ILNRT, IRNRT, INI, MTI, MLI, MMI, MPI, MTO, MLO, MMO, IN0,
     2    MT0, ML0, MM0, MP0, NBN, NDCS, NDRVS, NDRVM, NFDCM, NBN0,
     3    NDCS0, NDRVS0, NDRVM0, NFDCM0, A, FAC, XARR, VEC0, SVFDC,
     4    SVFDC0, NVI, IHNALP, IHNBET, IHNGAM, TVHI, CVI )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, N1, N2, KL, KU, KLE, IMF, ILNRP, IRNRP, ILNRT,
     1        IRNRT, INI( * ), MTI( * ), MLI( * ), MMI( * ),
     2        MPI( * ), MTO( * ), MLO( * ), MMO( * ), IN0( * ),
     3        MT0( * ), ML0( * ), MM0( * ), MP0( * )
      INTEGER NBN, NDCS, NDRVS, NDRVM, NFDCM, NBN0, NDCS0, NDRVS0,
     1        NDRVM0, NFDCM0, NVI, IHNALP( * ), IHNBET( * ),
     2        IHNGAM( * )
      CHARACTER *(4) TVHI( * )
      DOUBLE PRECISION A( N1, N2 ), FAC, XARR( NR ), VEC0( * ),
     1                 SVFDC( NFDCM , NR, NDRVM +1, NDCS  ),
     2                 SVFDC0( NFDCM0, NR, NDRVM0+1, NDCS0 ),
     3                 CVI( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IHA, IHB, IHG, LA, LB, LG, MA, MB, MG, ICSA, ICSB,
     1        ICSG, IQSTA, IQSTB, IQSTG, IOLDF, IHD, IHD0,
     2        IPARS( 7 ), NHA, NHB, NHG, INDG, INDSHC, IS, IS0,
     3        INDA, INDB
C     .
C     . ioldf = 2 --> new alpha function, existing beta function
C     .
      PARAMETER ( IOLDF = 2 )
C     .
      DOUBLE PRECISION DPARS( 1 ), LOW, WORK( 6 ),
     1                 CQTS, CQTT, CQSS, CQST, CSTQ, CSQS,
     2                 CSQT, CSSQ, CTTQ, CTQS, CTQT, CTSQ
      PARAMETER ( LOW = 1.0d-9 )
      EXTERNAL INVCVT
      CHARACTER *(4) REQVIT
      LOGICAL ONZIC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NHA = INI( 3 )
      NHB = IN0( 3 )
      NHG = INI( 3 )
C     .
C     . Early exit with zero value of multiplier
C     .
      IF ( DABS( FAC ).LT.LOW ) RETURN
C     .
C     . IHA is the 'alpha' harmonic.
C     . This is the velocity.
C     . In this case the velocity is in matrix solution
C     . vector and so loop around IHA from 1 to NHI
C     .
      DO IHA = 1, NHA
C       .
C       . Consider the case of alpha poloidal
C       .
        IF ( MTI( IHA ).EQ.1 ) THEN
          LA = MLI( IHA )
          IF ( MMI( IHA ).LT.0 ) THEN
            MA   = -MMI( IHA )
            ICSA = 2
          ELSE
            MA   = MMI( IHA )
            ICSA = 1
          ENDIF
          INDA   = INDSHC( LA, MA, ICSA )
C         .
          IS     = MPI( IHA )
C         .
          IQSTA = 1
C         .
C         . Now loop IHB around harmonics in VEC0
C         .
          DO IHB = 1, NHB
C           .
C           . Case of beta poloidal
C           .
            IF ( MT0( IHB ).EQ.1 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = ML0( IHB )
              IF ( MM0( IHB ).LT.0 ) THEN
                MB   = -MM0( IHB )
                ICSB = 2
              ELSE
                MB   = MM0( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
              IS0    = MP0( IHB )
C             .
              IQSTB = 3
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
C               . Case of Gamma toroidal
C               . See if we have C_{qts}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 ) THEN
                  REQVIT = 'CQTS'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CQTS, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 2
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 1
                  IHD0       = 3
C                 .
                  DPARS( 1 ) = CQTS
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{qtt}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.1 ) THEN
                  REQVIT = 'CQTT'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CQTT, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 3
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 0
                  IHD0       = 2
C                 .
                  DPARS( 1 ) = CQTT
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRP, IRNRP, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
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
            IF ( MT0( IHB ).EQ.2 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = ML0( IHB )
              IF ( MM0( IHB ).LT.0 ) THEN
                MB   = -MM0( IHB )
                ICSB = 2
              ELSE
                MB   = MM0( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
              IS0    = MP0( IHB )
C             .
              IQSTB = 2
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
C               . Case of Gamma toroidal
C               . See if we have C_{qss}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 ) THEN
                  REQVIT = 'CQSS'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CQSS, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 2
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
                  DPARS( 1 ) = CQSS
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{qst}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.1 ) THEN
                  REQVIT = 'CQST'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CQST, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 3
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
                  DPARS( 1 ) = CQST
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRP, IRNRP, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
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
          IQSTA = 2
C         .
C         . Now loop IHB around harmonics in VEC0
C         .
          DO IHB = 1, NHB
C           .
C           . Case of beta poloidal
C           .
            IF ( MT0( IHB ).EQ.1 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = ML0( IHB )
              IF ( MM0( IHB ).LT.0 ) THEN
                MB   = -MM0( IHB )
                ICSB = 2
              ELSE
                MB   = MM0( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
              IS0    = MP0( IHB )
C             .
              IQSTB = 3
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
C               . Case of Gamma toroidal
C               . See if we have C_{stq}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 ) THEN
                  REQVIT = 'CSTQ'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CSTQ, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
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
                  IHD0       = 2
C                 .
                  DPARS( 1 ) = CSTQ
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
            ENDIF
C           .
C           . Case of beta toroidal
C           .
            IF ( MT0( IHB ).EQ.2 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = ML0( IHB )
              IF ( MM0( IHB ).LT.0 ) THEN
                MB   = -MM0( IHB )
                ICSB = 2
              ELSE
                MB   = MM0( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
              IS0    = MP0( IHB )
C             .
              IQSTB = 1
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
C               . Case of Gamma toroidal
C               . See if we have C_{sqs}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 ) THEN
                  REQVIT = 'CSQS'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CSQS, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 2
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
                  DPARS( 1 ) = CSQS
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{sqt}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.1 ) THEN
                  REQVIT = 'CSQT'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CSQT, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 3
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
                  DPARS( 1 ) = CSQT
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRP, IRNRP, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
              IQSTB = 2
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
C               . Case of Gamma toroidal
C               . See if we have C_{ssq}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 ) THEN
                  REQVIT = 'CSSQ'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CSSQ, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 1
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
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
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
        IF ( MTI( IHA ).EQ.2 ) THEN
          LA = MLI( IHA )
          IF ( MMI( IHA ).LT.0 ) THEN
            MA   = -MMI( IHA )
            ICSA = 2
          ELSE
            MA   = MMI( IHA )
            ICSA = 1
          ENDIF
          INDA   = INDSHC( LA, MA, ICSA )
C         .
          IS     = MPI( IHA )
C         .
          IQSTA = 3
C         .
C         . Now loop IHB around harmonics in VEC0
C         .
          DO IHB = 1, NHB
C           .
C           . Case of beta poloidal
C           .
            IF ( MT0( IHB ).EQ.1 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = ML0( IHB )
              IF ( MM0( IHB ).LT.0 ) THEN
                MB   = -MM0( IHB )
                ICSB = 2
              ELSE
                MB   = MM0( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
              IS0    = MP0( IHB )
C             .
              IQSTB = 3
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
C               . Case of Gamma toroidal
C               . See if we have C_{ttq}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 ) THEN
                  REQVIT = 'CTTQ'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CTTQ, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 1
                  IPARS( 1 ) = IQSTA
                  IPARS( 2 ) = IQSTB
                  IPARS( 3 ) = IQSTG
                  IPARS( 4 ) = IOLDF
                  IPARS( 5 ) = LA
                  IPARS( 6 ) = LB
                  IPARS( 7 ) = LG
C                 .
                  IHD        = 0
                  IHD0       = 2
C                 .
                  DPARS( 1 ) = CTTQ
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
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
            IF ( MT0( IHB ).EQ.2 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = ML0( IHB )
              IF ( MM0( IHB ).LT.0 ) THEN
                MB   = -MM0( IHB )
                ICSB = 2
              ELSE
                MB   = MM0( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
              IS0    = MP0( IHB )
C             .
C             . Evaluate q_{beta} in VF2
C             .
              IQSTB = 1
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
C               . Case of Gamma toroidal
C               . See if we have C_{tqs}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 ) THEN
                  REQVIT = 'CTQS'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CTQS, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 2
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
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{tqt}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.1 ) THEN
                  REQVIT = 'CTQT'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CTQT, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 3
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
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRP, IRNRP, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
              IQSTB = 2
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
C               . Case of Gamma toroidal
C               . See if we have C_{tsq}^{abg} interaction
C               .
                IF ( MTO( IHG ).EQ.2 ) THEN
                  REQVIT = 'CTSQ'
                  CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                   INDB, INDG, TVHI, REQVIT, CVI, CTSQ, ONZIC )
C                 .
                  IF ( ONZIC ) THEN
C                 .
                  IQSTG = 1
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
                  DPARS( 1 ) = CTSQ
C                 .
                  CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHA, IHG,
     1      INI, IHD, NBN, ILNRT, IRNRT, NFDCM, NR, NDCS, IS, NDRVS,
     2      NDRVM, IPARS, INVCVT, A, FAC, XARR, WORK, DPARS, SVFDC,
     3      IHB, IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0,
     4      VEC0, NFDCM0, SVFDC0 )
C                 .
                  ENDIF
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
C     .
      RETURN
      END
C*********************************************************************

