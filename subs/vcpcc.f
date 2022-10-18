C*********************************************************************
C subroutine Vector Cross Product Coefficients Calculate *************
C            -      -     -       -            -         *************
C Steve Gibbons Tue Nov 23 07:51:01 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Calculates the coefficients for the non-linear interaction         C
C                                                                    C
C curl( v_{alpha}  x  B_{beta} )                                     C
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
C Interaction] is given as 'CABG' where each of A, B and G is        C
C replaced by Q, S or T depending upon the nature of the vector      C
C harmonics.                                                         C
C                                                                    C
C  The explicit formulae for the interactions 'CABG' is given by     C
C the equations (B.45) through to (B.50) on page 188 of my thesis.   C
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
C  velocity: (3 is temperature but not relevant here),               C
C  4 = poloidal magnetic field and 5 = toroidal magnetic field.      C
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
C                                                                    C
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
      SUBROUTINE VCPCC( NVI, MAXNVI, IHNALP, IHNBET, IHNGAM, NHA, 
     1    MTA, MLA, MMA, NHB, MTB, MLB, MMB, NHG, MTG, MLG, MMG, 
     2    LH, NTHP, NPHP, MMAX, TVHI, CVI, GAUX, GAUW, PA, DPA,
     3    FF1, FF2, FF3, VF1, VF2, VF3, QST )
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
     2                 FF1( 2*NPHP ), FF2( 2*NPHP ), FF3( 2*NPHP )
      DOUBLE PRECISION VF1( NPHP, NTHP, 3), VF2( NPHP, NTHP, 3),
     1                 VF3( NPHP, NTHP, 3 ), QST( LH*( LH + 2), 3 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IHA, IHB, IHG, LA, LB, LG, MA, MB, MG, ICSA, ICSB,
     1        ICSG, IQSTA, IQSTB, INDA, INDB, INDG, INDSHC
C     .
      DOUBLE PRECISION ZCOEF, LOW, 
     1                 CQTS, CQTT, CQSS, CQST, CSTQ, CSQS,
     2                 CSQT, CSSQ, CTTQ, CTQS, CTQT, CTSQ
      PARAMETER ( LOW = 1.0d-9 )
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      OK = .TRUE.
C     .
C     . Check that NVI is valid.
C     .
      IF ( NVI.LT.0 .OR. NVI.GT.MAXNVI ) THEN
        PRINT *,' Subroutine VCPCC.'
        PRINT *,' NVI = ', NVI,' MAXNVI = ', MAXNVI 
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . IHA is the 'alpha' harmonic.
C     . This is the velocity.
C     . vector and so loop around IHA from 1 to NHA
C     .
      DO IHA = 1, NHA
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
C         . Now loop IHB around beta harmonics
C         .
          DO IHB = 1, NHB
C           .
C           . Case of beta toroidal
C           .
            IF ( MTB( IHB ).EQ.5 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLB( IHB )
              IF ( MMB( IHB ).LT.0 ) THEN
                MB   = -MMB( IHB )
                ICSB = 2
              ELSE
                MB   = MMB( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
C             . Evaluate t_{beta} in VF2
C             .
              IQSTB = 3
C             .
              CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                      PA, DPA, NTHP, NPHP, LH )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX )
C             .
C             . QST now contains C_{qts}^{abg} in comp. 2
C             . QST now contains C_{qtt}^{abg} in comp. 3
C             .
C             . Loop around gamma harmonics
C             .
              DO IHG = 1, NHG
                LG = MLG( IHG )
                IF ( MMG( IHG ).LT.0 ) THEN
                  MG   = -MMG( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMG( IHG )
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
                IF ( MTG( IHG ).EQ.5 .AND. DABS( CQTS ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CQTS
                    TVHI( NVI )   = 'CQTS'
                  ENDIF
C                 .
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{qtt}^{abg} interaction
C               .
                IF ( MTG( IHG ).EQ.4 .AND. DABS( CQTT ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CQTT
                    TVHI( NVI )   = 'CQTT'
                  ENDIF
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
            ENDIF
C           .
C           .
C           . Case of beta poloidal
C           .
            IF ( MTB( IHB ).EQ.4 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLB( IHB )
              IF ( MMB( IHB ).LT.0 ) THEN
                MB   = -MMB( IHB )
                ICSB = 2
              ELSE
                MB   = MMB( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
C             . Evaluate s_{beta} in VF2
C             . (No need to evaluate q_{beta}
C             . - there is no Q x Q )
C             .
              IQSTB = 2
C             .
              CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                      PA, DPA, NTHP, NPHP, LH )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX )
C             .
C             . QST now contains C_{qss}^{abg} in comp. 2
C             . QST now contains C_{qst}^{abg} in comp. 3
C             .
C             . Loop around gamma harmonics
C             .
              DO IHG = 1, NHG
                LG = MLG( IHG )
                IF ( MMG( IHG ).LT.0 ) THEN
                  MG   = -MMG( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMG( IHG )
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
                IF ( MTG( IHG ).EQ.5 .AND. DABS( CQSS ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CQSS
                    TVHI( NVI )   = 'CQSS'
                  ENDIF
C                 .
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{qst}^{abg} interaction
C               .
                IF ( MTG( IHG ).EQ.4 .AND. DABS( CQST ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CQST
                    TVHI( NVI )   = 'CQST'
                  ENDIF
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
          CALL SHVECT ( LA, MA, ICSA, IQSTA, VF1, GAUX, PA, DPA,
     1                  NTHP, NPHP, LH )
C         .
C         . Now loop IHB around beta harmonics
C         .
          DO IHB = 1, NHB
C           .
C           . Case of beta toroidal
C           .
            IF ( MTB( IHB ).EQ.5 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLB( IHB )
              IF ( MMB( IHB ).LT.0 ) THEN
                MB   = -MMB( IHB )
                ICSB = 2
              ELSE
                MB   = MMB( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
C             . Evaluate t_{beta} in VF2
C             .
              IQSTB = 3
C             .
              CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                      PA, DPA, NTHP, NPHP, LH )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX )
C             .
C             . QST now contains C_{stq}^{abg} in comp. 1
C             .
C             . Loop around gamma harmonics
C             .
              DO IHG = 1, NHG
                LG = MLG( IHG )
                IF ( MMG( IHG ).LT.0 ) THEN
                  MG   = -MMG( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMG( IHG )
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
                IF ( MTG( IHG ).EQ.5 .AND. DABS( CSTQ ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CSTQ
                    TVHI( NVI )   = 'CSTQ'
                  ENDIF
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
            ENDIF
C           .
C           . Case of beta poloidal
C           .
            IF ( MTB( IHB ).EQ.4 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLB( IHB )
              IF ( MMB( IHB ).LT.0 ) THEN
                MB   = -MMB( IHB )
                ICSB = 2
              ELSE
                MB   = MMB( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
C             . Evaluate q_{beta} in VF2
C             .
              IQSTB = 1
C             .
              CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                      PA, DPA, NTHP, NPHP, LH )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX )
C             .
C             . QST now contains C_{sqs}^{abg} in comp. 2
C             . QST now contains C_{sqt}^{abg} in comp. 3
C             .
C             . Loop around gamma harmonics
C             .
              DO IHG = 1, NHG
                LG = MLG( IHG )
                IF ( MMG( IHG ).LT.0 ) THEN
                  MG   = -MMG( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMG( IHG )
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
                IF ( MTG( IHG ).EQ.5 .AND. DABS( CSQS ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CSQS
                    TVHI( NVI )   = 'CSQS'
                  ENDIF
C                 .
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{sqt}^{abg} interaction
C               .
                IF ( MTG( IHG ).EQ.4 .AND. DABS( CSQT ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CSQT
                    TVHI( NVI )   = 'CSQT'
                  ENDIF
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
              CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                      PA, DPA, NTHP, NPHP, LH )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX )
C             .
C             . QST now contains C_{ssq}^{abg} in comp. 1
C             .
C             . Loop around gamma harmonics
C             .
              DO IHG = 1, NHG
                LG = MLG( IHG )
                IF ( MMG( IHG ).LT.0 ) THEN
                  MG   = -MMG( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMG( IHG )
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
                IF ( MTG( IHG ).EQ.5 .AND. DABS( CSSQ ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CSSQ
                    TVHI( NVI )   = 'CSSQ'
                  ENDIF
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
        IF ( MTA( IHA ).EQ.2 ) THEN
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
C         . Now evaluate t_{alpha} in VF1
C         .
          IQSTA = 3
C         .
          CALL SHVECT ( LA, MA, ICSA, IQSTA, VF1, GAUX, PA, DPA,
     1                  NTHP, NPHP, LH )
C         .
C         . Now loop IHB around beta harmonics
C         .
          DO IHB = 1, NHB
C           .
C           . Case of beta toroidal
C           .
            IF ( MTB( IHB ).EQ.5 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLB( IHB )
              IF ( MMB( IHB ).LT.0 ) THEN
                MB   = -MMB( IHB )
                ICSB = 2
              ELSE
                MB   = MMB( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
C             . Evaluate t_{beta} in VF2
C             .
              IQSTB = 3
C             .
              CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                      PA, DPA, NTHP, NPHP, LH )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX )
C             .
C             . QST now contains C_{ttq}^{abg} in comp. 1
C             .
C             . Loop around gamma harmonics
C             .
              DO IHG = 1, NHG
                LG = MLG( IHG )
                IF ( MMG( IHG ).LT.0 ) THEN
                  MG   = -MMG( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMG( IHG )
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
                IF ( MTG( IHG ).EQ.5 .AND. DABS( CTTQ ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CTTQ
                    TVHI( NVI )   = 'CTTQ'
                  ENDIF
C                 .
                ENDIF
C               .
              ENDDO
C             (end loop ihg = 1, nhg)
C             .
            ENDIF
C           .
C           .
C           . Case of beta poloidal
C           .
            IF ( MTB( IHB ).EQ.4 ) THEN
C             .
C             . Calculate indices of velocity harm. beta
C             .
              LB = MLB( IHB )
              IF ( MMB( IHB ).LT.0 ) THEN
                MB   = -MMB( IHB )
                ICSB = 2
              ELSE
                MB   = MMB( IHB )
                ICSB = 1
              ENDIF
              INDB   = INDSHC( LB, MB, ICSB )
C             .
C             . Evaluate q_{beta} in VF2
C             .
              IQSTB = 1
C             .
              CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                      PA, DPA, NTHP, NPHP, LH )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX )
C             .
C             . QST now contains C_{tqs}^{abg} in comp. 2
C             . QST now contains C_{tqt}^{abg} in comp. 3
C             .
C             . Loop around gamma harmonics
C             .
              DO IHG = 1, NHG
                LG = MLG( IHG )
                IF ( MMG( IHG ).LT.0 ) THEN
                  MG   = -MMG( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMG( IHG )
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
                IF ( MTG( IHG ).EQ.5 .AND. DABS( CTQS ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CTQS
                    TVHI( NVI )   = 'CTQS'
                  ENDIF
C                 .
                ENDIF
C               .
C               . Case of Gamma poloidal
C               . See if we have C_{tqt}^{abg} interaction
C               .
                IF ( MTG( IHG ).EQ.4 .AND. DABS( CTQT ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CTQT
                    TVHI( NVI )   = 'CTQT'
                  ENDIF
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
              CALL SHVECT ( LB, MB, ICSB, IQSTB, VF2, GAUX,
     1                      PA, DPA, NTHP, NPHP, LH )
C             .
C             . Calculate VF3 = VF1 x VF2
C             .
              CALL VFCP ( VF1, VF2, VF3, NPHP, NTHP )
C             .
C             . Transform back into QST coefficients
C             .
              CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FF1, FF2,
     1                    FF3, ZCOEF, LH, NTHP, NPHP, MMAX )
C             .
C             . QST now contains C_{tsq}^{abg} in comp. 1
C             .
C             . Loop around gamma harmonics
C             .
              DO IHG = 1, NHG
                LG = MLG( IHG )
                IF ( MMG( IHG ).LT.0 ) THEN
                  MG   = -MMG( IHG )
                  ICSG = 2
                ELSE
                  MG   = MMG( IHG )
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
                IF ( MTG( IHG ).EQ.5 .AND. DABS( CTSQ ).GT.LOW ) THEN
C                 .
                  CALL CNTRIC( NVI, MAXNVI, OK )
                  IF ( OK ) THEN
                    IHNALP( NVI ) = INDA
                    IHNBET( NVI ) = INDB
                    IHNGAM( NVI ) = INDG
                    CVI( NVI )    = CTSQ
                    TVHI( NVI )   = 'CTSQ'
                  ENDIF
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
      IF ( OK ) RETURN
      PRINT *,' Subroutine VCPCC.'
      PRINT *,' Number of required interactions = ', NVI
      PRINT *,' Maximum interactions allowed    = ', MAXNVI
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************

