C*********************************************************************
C subroutine Recalled coefficient Velocity_0 . Gradient of Theta Add *
C            -                    -        -   -           -     -   *
C Steve Gibbons Tue Nov 23 07:51:01 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Adds the v_0 . \nabla Theta terms to the matrix for solving        C
C the heat equation. The v_0 terms are contained in a vector         C
C VEC0 which has the same grid spacing as that in the matrix         C
C (NR grid points defined by the array XARR).                        C
C                                                                    C
C Unlike IV0GTA which evaluates the vector interactions in situ,     C
C RV0GTA searches through the arrays IHNALP, IHNBET, IHNGAM and TVHI C
C (as pre-calculated by VSPCC) with the help of VICEXR for           C
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
C     ILN       : Left-most node at which to start adding rows to    C
C                   matrix.                                          C
C     IRN       : Right-most node at which to start adding rows to   C
C                   matrix.                                          C
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
      SUBROUTINE RV0GTA( NR, N1, N2, KL, KU, KLE, IMF, ILN, IRN, INI,
     1    MTI, MLI, MMI, MPI, MTO, MLO, MMO, IN0, MT0, ML0, MM0, MP0,
     2    NBN, NDCS, NDRVS, NDRVM, NFDCM, NBN0, NDCS0, NDRVS0, NDRVM0, 
     3    NFDCM0, A, FAC, XARR, VEC0, SVFDC, SVFDC0, 
     4    NVI, IHNALP, IHNBET, IHNGAM, TVHI, CVI )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, N1, N2, KL, KU, KLE, IMF, ILN, IRN, INI( * ),
     1        MTI( * ), MLI( * ), MMI( * ), MPI( * ), MTO( * ), 
     2        MLO( * ), MMO( * ), IN0( * ), MT0( * ), ML0( * ),
     3        MM0( * ), MP0( * )
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
     1        ICSG, IQSTA, IQSTB, IOLDF, IHD, IHD0, IPARS( 5 ),
     2        NHA, NHB, NHG, INDG, INDSHC, IS, IS0, INDA, INDB
C     .
C     . ioldf = 1 --> existing velocity function: new temp. func.
C     .
      PARAMETER ( IOLDF = 1 )
C     .
      DOUBLE PRECISION DPARS( 1 ), LOW, SQQ, SSS, STS,
     1                 WORK( 6 )
      PARAMETER ( LOW = 1.0d-9 )
      EXTERNAL INVGTT
      CHARACTER *(4) REQVIT
      LOGICAL ONZIC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
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
C     . In this case the velocity is in VEC0 and so
C     . loop around IHA from 1 to NH0
C     .
      DO IHA = 1, NHA
C       .
C       . Go onto next alpha harmonic if it is not a velocity term
C       .
        IF ( MT0( IHA ).NE.1 .AND. MT0( IHA ).NE.2 ) GOTO 50
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
          INDA   = INDSHC( LA, MA, ICSA )
C         .
          IS0    = MP0( IHA )
C         .
          IQSTA = 1
C         .
C         . Now loop IHB around harmonics in solution
C         . vector and for temperature
C         .
          DO IHB = 1, NHB
C           .
C           . Go onto next beta harmonic if it is not temperature
C           .
            IF ( MTI( IHB ).NE.3 ) GOTO 51
C           .
C           . Calculate indices of temperature harm. beta
C           .
            LB = MLI( IHB )
            IF ( MMI( IHB ).LT.0 ) THEN
              MB   = -MMI( IHB )
              ICSB = 2
            ELSE
              MB   = MMI( IHB )
              ICSB = 1
            ENDIF
            INDB   = INDSHC( LB, MB, ICSB )
C           .
            IS     = MPI( IHB )
C           .
            IQSTB = 1
C           .
C           . Now loop 'gamma' harmonic around the
C           . output harmonics in our matrix.
C           .
            DO IHG = 1, NHG
C             .
C             . Calculate indices of temperature harm. gamma
C             .
              IF ( MTO( IHG ).NE.3 ) GOTO 52
C             .
              LG = MLO( IHG )
              IF ( MMO( IHG ).LT.0 ) THEN
                MG   = -MMO( IHG )
                ICSG = 2
              ELSE
                MG   = MMO( IHG )
                ICSG = 1
              ENDIF
C             .
              INDG   = INDSHC( LG, MG, ICSG )
C             .
              REQVIT = 'SPQQ'
              CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, SQQ, ONZIC )
C             .
              IF ( .NOT. ONZIC ) GOTO 52
C             .
C             . OK - we have an S^{abg}_{qq} interaction
C             .
              IHD  = 1
              IHD0 = 0
C             .
              DPARS( 1 ) = SQQ
C             .
              IPARS( 1 ) = IQSTA
              IPARS( 2 ) = IQSTB
              IPARS( 3 ) = IOLDF
              IPARS( 4 ) = LA
              IPARS( 5 ) = LB
C             .
              CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG, INI,
     1      IHD, NBN, ILN, IRN, NFDCM, NR, NDCS, IS, NDRVS, NDRVM,
     2      IPARS, INVGTT, A, FAC, XARR, WORK, DPARS, SVFDC, IHA,
     3      IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0, VEC0,
     4      NFDCM0, SVFDC0 )
C             .
 52         CONTINUE
C           (if mto( ihg ).ne.3 ) goto 52)
C           (if dabs( sqq ).lt.low ) goto 52)
            ENDDO
C           .
 51       CONTINUE
C         (if mti( ihb ).ne.3 ) goto 51)
          ENDDO
C         .
          IQSTA = 2
C         .
C         . Now loop IHB around harmonics in solution
C         . vector and for temperature
C         .
          DO IHB = 1, NHB
C           .
C           . Go onto next beta harmonic if it is not temperature
C           .
            IF ( MTI( IHB ).NE.3 ) GOTO 53
C           .
C           . Also go on to next harmonic if this is the
C           . monopole. This is because we will evaluate
C           . s_{beta} in VF2 which is ofcourse zero for
C           . l = 0.
C           .
            LB = MLI( IHB )
            IF ( LB.EQ.0 ) GOTO 53
C           .
C           . Calculate indices of temperature harm. beta
C           .
            IF ( MMI( IHB ).LT.0 ) THEN
              MB   = -MMI( IHB )
              ICSB = 2
            ELSE
              MB   = MMI( IHB )
              ICSB = 1
            ENDIF
            INDB   = INDSHC( LB, MB, ICSB )
C           .
            IS     = MPI( IHB )
C           .
            IQSTB = 2
C           .
C           . Now loop 'gamma' harmonic around the
C           . output harmonics in our matrix.
C           .
            DO IHG = 1, NHG
C             .
C             . Calculate indices of temperature harm. gamma
C             .
              IF ( MTO( IHG ).NE.3 ) GOTO 54
C             .
              LG = MLO( IHG )
              IF ( MMO( IHG ).LT.0 ) THEN
                MG   = -MMO( IHG )
                ICSG = 2
              ELSE
                MG   = MMO( IHG )
                ICSG = 1
              ENDIF
C             .
              INDG   = INDSHC( LG, MG, ICSG )
C             .
              REQVIT = 'SPSS'
              CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, SSS, ONZIC )
C             .
              IF ( .NOT. ONZIC ) GOTO 54
C             .
C             . OK - we have an S^{abg}_{ss} interaction
C             .
              IHD  = 0
              IHD0 = 1
C             .
              DPARS( 1 ) = SSS
C             .
              IPARS( 1 ) = IQSTA
              IPARS( 2 ) = IQSTB
              IPARS( 3 ) = IOLDF
              IPARS( 4 ) = LA
              IPARS( 5 ) = LB
C             .
              CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG, INI,
     1      IHD, NBN, ILN, IRN, NFDCM, NR, NDCS, IS, NDRVS, NDRVM,
     2      IPARS, INVGTT, A, FAC, XARR, WORK, DPARS, SVFDC, IHA,
     3      IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0, VEC0,
     4      NFDCM0, SVFDC0 )
C             .
 54         CONTINUE
C           (if mto( ihg ).ne.3 ) goto 54)
C           (if dabs( sss ).lt.low ) goto 54)
            ENDDO
C           .
 53       CONTINUE
C         (if mti( ihb ).ne.3 ) goto 53)
C         (if lb.eq.0 ) goto 53)
          ENDDO
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
          INDA   = INDSHC( LA, MA, ICSA )
C         .
          IS0    = MP0( IHA )
C         .
          IQSTA = 3
C         .
C         . Now loop IHB around harmonics in solution
C         . vector and for temperature
C         .
          DO IHB = 1, NHB
C           .
C           . Go onto next beta harmonic if it is not temperature
C           .
            IF ( MTI( IHB ).NE.3 ) GOTO 55
C           .
C           . Also go on to next harmonic if this is the
C           . monopole. This is because we will evaluate
C           . s_{beta} in VF2 which is ofcourse zero for
C           . l = 0.
C           .
            LB = MLI( IHB )
            IF ( LB.EQ.0 ) GOTO 55
C           .
C           . Calculate indices of temperature harm. beta
C           .
            IF ( MMI( IHB ).LT.0 ) THEN
              MB   = -MMI( IHB )
              ICSB = 2
            ELSE
              MB   = MMI( IHB )
              ICSB = 1
            ENDIF
            INDB   = INDSHC( LB, MB, ICSB )
C           .
            IS     = MPI( IHB )
C           .
            IQSTB = 2
C           .
C           . Now loop 'gamma' harmonic around the
C           . output harmonics in our matrix.
C           .
            DO IHG = 1, NHG
C             .
C             . Calculate indices of temperature harm. gamma
C             .
              IF ( MTO( IHG ).NE.3 ) GOTO 56
C             .
              LG = MLO( IHG )
              IF ( MMO( IHG ).LT.0 ) THEN
                MG   = -MMO( IHG )
                ICSG = 2
              ELSE
                MG   = MMO( IHG )
                ICSG = 1
              ENDIF
C             .
              INDG   = INDSHC( LG, MG, ICSG )
C             .
              REQVIT = 'SPTS'
              CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, STS, ONZIC )
C             .
              IF ( .NOT. ONZIC ) GOTO 56
C             .
C             . OK - we have an S^{abg}_{ts} interaction
C             .
              IHD  = 0
              IHD0 = 0
C             .
              DPARS( 1 ) = STS
C             .
              IPARS( 1 ) = IQSTA
              IPARS( 2 ) = IQSTB
              IPARS( 3 ) = IOLDF
              IPARS( 4 ) = LA
              IPARS( 5 ) = LB
C             .
              CALL INNLCA( N1, N2, KL, KU, KLE, IMF, IHB, IHG, INI,
     1      IHD, NBN, ILN, IRN, NFDCM, NR, NDCS, IS, NDRVS, NDRVM,
     2      IPARS, INVGTT, A, FAC, XARR, WORK, DPARS, SVFDC, IHA,
     3      IHD0, NDCS0, IS0, NBN0, IN0, NDRVS0, NDRVM0, VEC0,
     4      NFDCM0, SVFDC0 )
C             .
 56         CONTINUE
C           (if mto( ihg ).ne.3 ) goto 56)
C           (if dabs( sts ).lt.low ) goto 56)
            ENDDO
C           .
 55       CONTINUE
C         (if mti( ihb ).ne.3 ) goto 55)
C         (if lb.eq.0 ) goto 55)
          ENDDO
C         .
        ENDIF
C       .
 50   CONTINUE
C     ( if mt0( iha ).ne.1 .and. mt0( iha ).ne.2 ) goto 50 )
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

