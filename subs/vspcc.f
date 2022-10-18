C*********************************************************************
C subroutine Vector Scalar Product Coefficients Calculate ************
C            -      -      -       -            -         ************
C Steve Gibbons Mon Jan 17 18:53:30 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Calculates the coefficients for the non-linear interaction         C
C                                                                    C
C v_{alpha} . Theta_{beta}                                           C
C                                                                    C
C The output takes the form of 5 arrays; IHNALP, IHNBET, IHNGAM,     C
C CVI and TVHI.                                                      C
C                                                                    C
C  IHNALP( in ) and IHNBET( in ) [in = interaction number]           C
C refer to the harmonic number (as defined by the function INDSHC)   C
C of the alpha and beta harmonics respectively.                      C
C                                                                    C
C  If a given reaction forming a third vector harmonic, gamma,       C
C is non-zero then the coefficient C^{abg}_{ABG} is stored in        C
C CVI( in ) [Coefficient of Vector Interaction].                     C
C  IHNGAM( in ) is the number defined by INDSHC for this gamma.      C
C                                                                    C
C  The corresponding element of TVHI [Type of Vector Harmonic        C
C Interaction] is given as 'SPAB' where each of A and B is           C
C replaced by Q, S or T depending upon the nature of the vector      C
C harmonics.                                                         C
C                                                                    C
C  The explicit formulae for the interactions 'SPAB' is given by     C
C the equations (B.39) through to (B.41) on page 188 of my thesis.   C
C                                                                    C
C NVI is the Number of Vector Interactions and, as an input          C
C parameter, indicates how many elements are already stored in the   C
C five arrays. (NVI=0 means that no such routine has yet been        C
C called; NVI .gt. 0 means that the coefficients found by this       C
C routine will begin at NVI + 1).                                    C
C                                                                    C
C MAXNVI is the maximum permitted number of vector interactions.     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NVI      : Number of vector interactions.                      C
C     MAXNVI   : Maximum number of vector interactions.              C
C                This is an upper limit for NVI and the dimension    C
C                of the arrays IHNALP, IHNBET, IHNGAM, CVI and       C
C                TVHI.                                               C
C                                                                    C
C     IHNALP   : Number of alpha harmonics. Dim ( MAXNVI )           C
C     IHNBET   : Number of beta harmonics. Dim ( MAXNVI )            C
C     IHNGAM   : Number of gamma harmonics. Dim ( MAXNVI )           C
C                                                                    C
C   [Key for MTA, MTB, MTG :- 1 = poloidal velocity, 2 = toroidal    C
C  velocity: 3, 4 and 5 are temp, pol mag. field and tor mag. f.]    C
C                                                                    C
C     NHA      : Number of alpha harmonics.                          C
C     MTA      : Array length ( * ) - atleast length NHA             C
C                  See above for key. (corresponds to alpha vec.)    C
C     MLA      : Array length ( * ) - atleast length NHA             C
C                  Sph. harm. degree, l.                             C
C     MMA      : Array length ( * ) - atleast length NHA             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     NHB      : Number of beta harmonics.                           C
C     MTB      : Array length ( * ) - atleast length NHB             C
C                  See above for key. (corresponds to beta vector)   C
C     MLB      : Array length ( * ) - atleast length NHB             C
C                  Sph. harm. degree, l.                             C
C     MMB      : Array length ( * ) - atleast length NHB             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     NHG      : Number of gamma harmonics.                          C
C     MTG      : Array length ( * ) - atleast length NHG             C
C                  See above for key. (corresponds to gamma vec.)    C
C     MLG      : Array length ( * ) - atleast length NHG             C
C                  Sph. harm. degree, l.                             C
C     MMG      : Array length ( * ) - atleast length NHG             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NTHP      : The number of theta points.                        C
C     NPHP      : The number of phi points.                          C
C     MMAX      : Maximum sph. harmonic order, m.                    C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. (MAXNVI).    C
C                 = 'CQSS', 'CQST' etc. according to the corresp.    C
C                 vector interaction.                                C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     CVI       : Coefficient of vector interaction. Dim ( MAXNVI )  C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     FF1       : Work array - dim. (2*NPHP)                         C
C     VF1       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF2       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     SF        : Work array - dim. ( NPHP, NTHP )                   C
C     SHC       : Work array - dim. ( LH * ( LH + 2 ) )              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VSPCC( NVI, MAXNVI, IHNALP, IHNBET, IHNGAM, NHA,
     1    MTA, MLA, MMA, NHB, MTB, MLB, MMB, NHG, MTG, MLG, MMG,
     2    LH, NTHP, NPHP, MMAX, TVHI, CVI, GAUX, GAUW, PA, DPA,
     3    FF1, VF1, VF2, SF, SHC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NVI, MAXNVI, IHNALP( MAXNVI ), IHNBET( MAXNVI ),
     1        IHNGAM( MAXNVI ), NHA, MTA( * ), MLA( * ), MMA( * ),
     2        NHB, MTB( * ), MLB( * ), MMB( * ), NHG, MTG( * ),
     3        MLG( * ), MMG( * ), LH, NTHP, NPHP, MMAX
      CHARACTER *(4) TVHI( MAXNVI )
      DOUBLE PRECISION CVI( MAXNVI ), GAUX( NTHP ), GAUW( NTHP )
      DOUBLE PRECISION PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     1                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     2                 FF1( 2*NPHP )
      DOUBLE PRECISION VF1( NPHP, NTHP, 3), VF2( NPHP, NTHP, 3),
     1                 SF( NPHP, NTHP ), SHC( LH*( LH + 2) )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IHA, IHB, IHG, LA, LB, LG, MA, MB, MG, ICSA, ICSB,
     1        ICSG, IQSTA, IQSTB, INDA, INDB, INDG, INDSHC
C     .
      LOGICAL OK
      DOUBLE PRECISION ZCOEF, LOW, SQQ, SSS, STS
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      OK = .TRUE.
C     .
C     . IHA is the 'alpha' harmonic.
C     . This is the velocity.
C     .
      DO IHA = 1, NHA
C       .
C       . Go onto next alpha harmonic if it is not a velocity term
C       .
        IF ( MTA( IHA ).NE.1 .AND. MTA( IHA ).NE.2 ) GOTO 50
C       .
C       . Consider the case of alpha poloidal
C       .
        IF ( MTA( IHA ).EQ.1 ) THEN
          LA = MLA( IHA )
          IF ( MMA( IHA ).LT.0 ) THEN
            MA   = -MMA( IHA )
            ICSA = 2
          ELSE
            MA   = MMA( IHA )
            ICSA = 1
          ENDIF
          INDA   = INDSHC( LA, MA, ICSA )
C         .
C         . Now evaluate q_{alpha} in VF1
C         .
          IQSTA = 1
C         .
          CALL SHVECT ( LA, MA, ICSA, IQSTA, VF1, GAUX, PA, DPA,
     1                  NTHP, NPHP, LH )
C         .
C         . Now loop IHB around harmonics in solution
C         . vector and for temperature
C         .
          DO IHB = 1, NHB
C           .
C           . Go onto next beta harmonic if it is not temperature
C           .
            IF ( MTB( IHB ).NE.3 ) GOTO 51
C           .
C           . Calculate indices of temperature harm. beta
C           .
            LB = MLB( IHB )
            IF ( MMB( IHB ).LT.0 ) THEN
              MB   = -MMB( IHB )
              ICSB = 2
            ELSE
              MB   = MMB( IHB )
              ICSB = 1
            ENDIF
            INDB   = INDSHC( LB, MB, ICSB )
C           .
C           . Evaluate q_{beta} in VF2
C           .
            IQSTB = 1
C           .
            CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX, PA, DPA,
     1                    NTHP, NPHP, LH )
C           .
C           . Calculate scalar product in SF
C           .
            CALL VFDP ( VF1, VF2, SF, NPHP, NTHP )
C           .
C           . Now transform this function into scalar
C           . spherical harmonic coefficients
C           .
            CALL FORSST ( SHC, SF, GAUW, PA, FF1, LH, MMAX,
     1                    NTHP, NPHP, ZCOEF )
C           .
C           . Now loop 'gamma' harmonic around the
C           . output harmonics in our matrix.
C           .
            DO IHG = 1, NHG
C             .
C             . Calculate indices of temperature harm. gamma
C             .
              IF ( MTG( IHG ).NE.3 ) GOTO 52
C             .
              LG = MLG( IHG )
              IF ( MMG( IHG ).LT.0 ) THEN
                MG   = -MMG( IHG )
                ICSG = 2
              ELSE
                MG   = MMG( IHG )
                ICSG = 1
              ENDIF
C             .
              INDG   = INDSHC( LG, MG, ICSG )
C             .
              IF ( INDG.EQ.0 ) THEN
                SQQ = ZCOEF
              ELSE
                SQQ = SHC( INDG )
              ENDIF
C             .
              IF ( DABS( SQQ ).LT.LOW ) GOTO 52
C             .
C             . OK - we have an S^{abg}_{qq} interaction
C             .
              CALL CNTRIC( NVI, MAXNVI, OK )
              IF ( OK ) THEN
                IHNALP( NVI ) = INDA
                IHNBET( NVI ) = INDB
                IHNGAM( NVI ) = INDG
                CVI( NVI )    = SQQ
                TVHI( NVI )   = 'SPQQ'
              ENDIF
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
C         . Now evaluate s_{alpha} in VF1
C         .
          IQSTA = 2
C         .
          CALL SHVECT ( LA, MA, ICSA, IQSTA, VF1, GAUX, PA, DPA,
     1                  NTHP, NPHP, LH )
C         .
C         . Now loop IHB around harmonics in solution
C         . vector and for temperature
C         .
          DO IHB = 1, NHB
C           .
C           . Go onto next beta harmonic if it is not temperature
C           .
            IF ( MTB( IHB ).NE.3 ) GOTO 53
C           .
C           . Also go on to next harmonic if this is the
C           . monopole. This is because we will evaluate
C           . s_{beta} in VF2 which is ofcourse zero for
C           . l = 0.
C           .
            LB = MLB( IHB )
            IF ( LB.EQ.0 ) GOTO 53
C           .
C           . Calculate indices of temperature harm. beta
C           .
            IF ( MMB( IHB ).LT.0 ) THEN
              MB   = -MMB( IHB )
              ICSB = 2
            ELSE
              MB   = MMB( IHB )
              ICSB = 1
            ENDIF
C           .
            INDB   = INDSHC( LB, MB, ICSB )
C           .
C           . Evaluate s_{beta} in VF2
C           .
            IQSTB = 2
C           .
            CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX, PA, DPA,
     1                    NTHP, NPHP, LH )
C           .
C           . Calculate scalar product in SF
C           .
            CALL VFDP ( VF1, VF2, SF, NPHP, NTHP )
C           .
C           . Now transform this function into scalar
C           . spherical harmonic coefficients
C           .
            CALL FORSST ( SHC, SF, GAUW, PA, FF1, LH, MMAX,
     1                    NTHP, NPHP, ZCOEF )
C           .
C           . Now loop 'gamma' harmonic around the
C           . output harmonics in our matrix.
C           .
            DO IHG = 1, NHG
C             .
C             . Calculate indices of temperature harm. gamma
C             .
              IF ( MTG( IHG ).NE.3 ) GOTO 54
C             .
              LG = MLG( IHG )
              IF ( MMG( IHG ).LT.0 ) THEN
                MG   = -MMG( IHG )
                ICSG = 2
              ELSE
                MG   = MMG( IHG )
                ICSG = 1
              ENDIF
C             .
              INDG   = INDSHC( LG, MG, ICSG )
C             .
              IF ( INDG.EQ.0 ) THEN
                SSS = ZCOEF
              ELSE
                SSS = SHC( INDG )
              ENDIF
C             .
              IF ( DABS( SSS ).LT.LOW ) GOTO 54
C             .
C             . OK - we have an S^{abg}_{ss} interaction
C             .
              CALL CNTRIC( NVI, MAXNVI, OK )
              IF ( OK ) THEN
                IHNALP( NVI ) = INDA
                IHNBET( NVI ) = INDB
                IHNGAM( NVI ) = INDG
                CVI( NVI )    = SSS
                TVHI( NVI )   = 'SPSS'
              ENDIF
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
        IF ( MTA( IHA ).EQ.2 ) THEN
          LA = MLA( IHA )
          IF ( MMA( IHA ).LT.0 ) THEN
            MA   = -MMA( IHA )
            ICSA = 2
          ELSE
            MA   = MMA( IHA )
            ICSA = 1
          ENDIF
C         .
          INDA   = INDSHC( LA, MA, ICSA )
C         .
C         . Now evaluate t_{alpha} in VF1
C         .
          IQSTA = 3
C         .
          CALL SHVECT ( LA, MA, ICSA, IQSTA, VF1, GAUX, PA, DPA,
     1                  NTHP, NPHP, LH )
C         .
C         . Now loop IHB around harmonics in solution
C         . vector and for temperature
C         .
          DO IHB = 1, NHB
C           .
C           . Go onto next beta harmonic if it is not temperature
C           .
            IF ( MTB( IHB ).NE.3 ) GOTO 55
C           .
C           . Also go on to next harmonic if this is the
C           . monopole. This is because we will evaluate
C           . s_{beta} in VF2 which is ofcourse zero for
C           . l = 0.
C           .
            LB = MLB( IHB )
            IF ( LB.EQ.0 ) GOTO 55
C           .
C           . Calculate indices of temperature harm. beta
C           .
            IF ( MMB( IHB ).LT.0 ) THEN
              MB   = -MMB( IHB )
              ICSB = 2
            ELSE
              MB   = MMB( IHB )
              ICSB = 1
            ENDIF
C           .
            INDB   = INDSHC( LB, MB, ICSB )
C           .
            IQSTB = 2
C           .
            CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX, PA, DPA,
     1                    NTHP, NPHP, LH )
C           .
C           . Calculate scalar product in SF
C           .
            CALL VFDP ( VF1, VF2, SF, NPHP, NTHP )
C           .
C           . Now transform this function into scalar
C           . spherical harmonic coefficients
C           .
            CALL FORSST ( SHC, SF, GAUW, PA, FF1, LH, MMAX,
     1                    NTHP, NPHP, ZCOEF )
C           .
C           . Now loop 'gamma' harmonic around the
C           . output harmonics in our matrix.
C           .
            DO IHG = 1, NHG
C             .
C             . Calculate indices of temperature harm. gamma
C             .
              IF ( MTG( IHG ).NE.3 ) GOTO 56
C             .
              LG = MLG( IHG )
              IF ( MMG( IHG ).LT.0 ) THEN
                MG   = -MMG( IHG )
                ICSG = 2
              ELSE
                MG   = MMG( IHG )
                ICSG = 1
              ENDIF
C             .
              INDG   = INDSHC( LG, MG, ICSG )
C             .
              IF ( INDG.EQ.0 ) THEN
                STS = ZCOEF
              ELSE
                STS = SHC( INDG )
              ENDIF
C             .
              IF ( DABS( STS ).LT.LOW ) GOTO 56
C             .
C             . OK - we have an S^{abg}_{ts} interaction
C             .
              CALL CNTRIC( NVI, MAXNVI, OK )
              IF ( OK ) THEN
                IHNALP( NVI ) = INDA
                IHNBET( NVI ) = INDB
                IHNGAM( NVI ) = INDG
                CVI( NVI )    = STS
                TVHI( NVI )   = 'SPTS'
              ENDIF
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
      IF ( OK ) RETURN
      PRINT *,' Subroutine VSPCC.'
      PRINT *,' Number of required interactions = ', NVI
      PRINT *,' Maximum interactions allowed    = ', MAXNVI
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************

