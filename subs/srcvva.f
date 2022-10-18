C*********************************************************************
C subroutine Stored derivative Recalled coefficient Curl Vector   ****
C            -                 -                    -    -        ****
C                                           cross curl Vector Add ****
C                                                      -      -   ****
C Steve Gibbons Mon Feb  7 10:12:34 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given a Vector, VA, which is stored in V0A (indexed by INARRA,     C
C MHTA, MHLA, MHMA - 1st and 2nd deriv.s in V1A and V2A resp. ),     C
C a vector, VB, which is stored in V0B (indexed by INARRB, MHTB,     C
C MHLB, MHMB - 1st, 2nd and 3rd deriv.s in V1B, V2B and V3B resp.)   C
C and the arrays                                                     C
C calculated by VCCPCG (i.e. IHNALP, IHNBET, IHNGAM, TVHI, CVI)      C
C to calculate the terms                                             C
C                                                                    C
C  curl (  VA . Grad ) VB = - curl ( VA x ( curl VB )  )             C
C                                                                    C
C and add a multiple FAC of these terms to                           C
C the appropriate elements of VECG (indexed by INARRG, MHTG, MHLG,   C
C MHMG).                                                             C
C                                                                    C
C This can then be used to calculate either the inertial terms       C
C or the Lorentz force in the vorticity equation.                    C
C                                                                    C
C The latter is chosen by selecting CHVMFF to 'MAG' and the former   C
C by selecting CHVMFF as 'VEL'.                                      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NVI      : Number of vector interactions.                      C
C     IHNALP   : Number of alpha harmonics. Dim ( * )                C
C                  (Output by subroutine VCCPCG).                    C
C     IHNBET   : Number of beta harmonics. Dim ( * )                 C
C                  (Output by subroutine VCCPCG).                    C
C     IHNGAM   : Number of gamma harmonics. Dim ( * )                C
C                  (Output by subroutine VCCPCG).                    C
C                                                                    C
C   [Key for MTA, MTB, MTG :- 1 = poloidal velocity, 2 = toroidal    C
C  velocity: 3, 4 and 5 are temp, pol mag. field and tor mag. f.]    C
C                                                                    C
C     INARRA   : Indexing array for gamma vector. Dim ( 3 )          C
C     MTA      : Array length ( * ) - atleast length NHA             C
C                  See above for key. (corresponds to alpha vec.)    C
C     MLA      : Array length ( * ) - atleast length NHA             C
C                  Sph. harm. degree, l.                             C
C     MMA      : Array length ( * ) - atleast length NHA             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C     INARRB   : Indexing array for beta  vector. Dim ( 3 )          C
C     MTB      : Array length ( * ) - atleast length NHB             C
C                  See above for key. (corresponds to beta vector)   C
C     MLB      : Array length ( * ) - atleast length NHB             C
C                  Sph. harm. degree, l.                             C
C     MMB      : Array length ( * ) - atleast length NHB             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C     INARRG   : Indexing array for gamma vector. Dim ( 3 )          C
C     MTG      : Array length ( * ) - atleast length NHG             C
C                  See above for key. (corresponds to gamma vec.)    C
C     MLG      : Array length ( * ) - atleast length NHG             C
C                  Sph. harm. degree, l.                             C
C     MMG      : Array length ( * ) - atleast length NHG             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     ILNP     : Left most node on which operation is to be          C
C                 performed for poloidal gamma.                      C
C     IRNP     : Right most node on which operation is to be         C
C                 performed for poloidal gamma.                      C
C     ILNT     : Left most node on which operation is to be          C
C                 performed for toroidal gamma.                      C
C     IRNT     : Right most node on which operation is to be         C
C                 performed for toroidal gamma.                      C
C                                                                    C
C     NR       : Number of radial grid nodes.                        C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. ( * ).       C
C                 = 'SPQQ', 'SPSS' etc. according to the corresp.    C
C                 vector interaction.                                C
C                  (Output by subroutine VCCPCG).                    C
C                                                                    C
C     CHVMFF    : Velocity/magnetic field select.                    C
C                 Set to either 'VEL' or 'MAG'                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     CVI       : Coefficient of vector interaction. Dim ( * ).      C
C                  (Output by subroutine VCCPCG).                    C
C     V0A       : 0 derivatives of alpha velocity. Dim ( * )         C
C     V1A       : 1 derivatives of alpha velocity. Dim ( * )         C
C     V2A       : 2 derivatives of alpha velocity. Dim ( * )         C
C     V0B       : 0 derivatives of beta velocity. Dim ( * )          C
C     V1B       : 1 derivatives of beta velocity. Dim ( * )          C
C     V2B       : 2 derivatives of beta velocity. Dim ( * )          C
C     V3B       : 3 derivatives of beta velocity. Dim ( * )          C
C     VECG      : Solution vector for gamma function. Dim ( * ).     C
C     XARR      : Dim ( NR ). Radial grid spacings.                  C
C     FAC       : Coefficient of term to be added to VECG.           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SRCVVA( NVI, IHNALP, IHNBET, IHNGAM, INARRA, MTA,
     1          MLA, MMA, INARRB, MTB, MLB, MMB, INARRG, MTG, MLG,
     2          MMG, ILNP, IRNP, ILNT, IRNT, NR, TVHI, CVI, V0A,
     3          V1A, V2A, V0B, V1B, V2B, V3B, VECG, XARR, FAC, CHVMFF)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NVI, IHNALP( * ), IHNBET( * ), IHNGAM( * ), INARRA(3),
     1        MTA( * ), MLA( * ), MMA( * ), INARRB(3),
     2        MTB( * ), MLB( * ), MMB( * ), INARRG(3), NR,
     3        MTG( * ), MLG( * ), MMG( * ), ILNP, IRNP, ILNT, IRNT
      CHARACTER *(4) TVHI( * )
      CHARACTER *(3) CHVMFF
      DOUBLE PRECISION CVI( * ), V0A( * ), V1A( * ), V2A( * ),
     1                 V0B( * ), V1B( * ), V2B( * ), V3B( * ),
     2                 VECG( * ), XARR( NR ), FAC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER LA, LB, IR, IHA, IHB, IHG, INDA, INDB, INDG,
     1        LOCG, INDFUN, ISHCIA, NRRA, NRRB, NRRG, NHA, NHB, NHG,
     2        IPOL, ITOR, LG
      DOUBLE PRECISION RAD, COEF, TERM, SDSVVT, LOW
      PARAMETER ( LOW = 1.0d-9 )
      CHARACTER *(4) REQVIT
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
C     . Early exit for zero value of FAC
C     .
      IF ( DABS( FAC ).LT.LOW ) RETURN
C     .
C     . Let's check validity of CHVMFF
C     .
      IF ( CHVMFF.NE.'VEL' .AND. CHVMFF.NE.'vel' .AND.
     1     CHVMFF.NE.'Vel' .AND. CHVMFF.NE.'MAG' .AND.
     2     CHVMFF.NE.'mag' .AND. CHVMFF.NE.'Mag'  ) THEN
        PRINT *,' Subroutine VCCPCG.'
        PRINT *,' CHVMFF = ', CHVMFF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( CHVMFF.EQ.'VEL' .OR. CHVMFF.EQ.'vel' .OR.
     1     CHVMFF.EQ.'Vel' ) THEN
        IPOL = 1
        ITOR = 2
      ENDIF
C     .
      IF ( CHVMFF.EQ.'MAG' .OR. CHVMFF.EQ.'mag' .OR.
     1     CHVMFF.EQ.'Mag' ) THEN
        IPOL = 4
        ITOR = 5
      ENDIF
C     .
C Just check the number of grid nodes in A, B, and G
C
      NRRA = INARRA( 2 )
      NRRB = INARRB( 2 )
      NRRG = INARRG( 2 )
      IF ( NRRA.NE.NR .OR. NRRB.NE.NR .OR. NRRG.NE.NR ) THEN
        PRINT *,' Subroutine SRCVVA. NR = ', NR
        PRINT *,' NRRA = ',NRRA,' NRRB = ',NRRB,' NRRG = ',NRRG
        PRINT *,' Program aborted,'
        STOP
      ENDIF
C
      NHA = INARRA( 3 )
      NHB = INARRB( 3 )
      NHG = INARRG( 3 )
C
C Begin loop around ALPHA vector harmonics
C
      DO IHA = 1, NHA
C       .
C       . Jump to 700 (end of loop) if A harmonic is not velocity
C       . or magnetic field
C       .
        IF ( MTA( IHA ).NE.IPOL .AND. MTA( IHA ).NE.ITOR ) GOTO 700
C       .
        INDA = ISHCIA( IHA, MLA, MMA )
        LA   = MLA( IHA )
C       .
C       . Now loop around BETA vector harmonics
C       .
        DO IHB = 1, NHB
C         .
C         . Jump to 701 (end of loop) if B harmonic is not vel.
C         . or magnetic field.
C         .
          IF ( MTB( IHB ).NE.IPOL .AND. MTB( IHB ).NE.ITOR ) GOTO 701
C         .
          INDB = ISHCIA( IHB, MLB, MMB )
          LB   = MLB( IHB )
C         .
C         . Now loop around GAMMA harmonics (vorticity.)
C         .
          DO IHG = 1, NHG
C           .
C           . Jump to 702 (end of loop) if G harm. not vort.
C           .
            IF ( MTG( IHG ).NE.1 .AND. MTG( IHG ).NE.2 ) GOTO 702
C           .
            INDG = ISHCIA( IHG, MLG, MMG )
            LG   = MLG( IHG )
C           .
C           . So our alpha harmonic is velocity/magnetic field
C           . Our beta harmonic is velocity/magnetic field
C           . Our gamma harmonic is vorticity
C           .
C           . First case: alpha POLOIDAL
C           .             beta  POLOIDAL
C           .             gamma POLOIDAL
C           .
            IF (     MTA( IHA ).EQ.IPOL       .AND.
     1               MTB( IHB ).EQ.IPOL       .AND.
     2               MTG( IHG ).EQ.1                ) THEN
C              .
C              . Test for C_qtt interaction
C              .
               REQVIT = 'CQTT'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNP, IRNP
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
            ENDIF
C           .
C           . Second case: alpha POLOIDAL
C           .              beta  POLOIDAL
C           .              gamma TOROIDAL
C           .
            IF (     MTA( IHA ).EQ.IPOL       .AND.
     1               MTB( IHB ).EQ.IPOL       .AND.
     2               MTG( IHG ).EQ.2                ) THEN
C              .
C              . Test for C_qts interaction
C              .
               REQVIT = 'CQTS'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNT, IRNT
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
C              . Test for C_stq interaction
C              .
               REQVIT = 'CSTQ'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNT, IRNT
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
            ENDIF
C           .
C           . Third case: alpha POLOIDAL
C           .             beta  TOROIDAL
C           .             gamma POLOIDAL
C           .
            IF (     MTA( IHA ).EQ.IPOL       .AND.
     1               MTB( IHB ).EQ.ITOR       .AND.
     2               MTG( IHG ).EQ.1                ) THEN
C              .
C              . Test for C_qst interaction
C              .
               REQVIT = 'CQST'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNP, IRNP
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
C              . Test for C_sqt interaction
C              .
               REQVIT = 'CSQT'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNP, IRNP
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
            ENDIF
C           .
C           . Fourth case: alpha POLOIDAL
C           .              beta  TOROIDAL
C           .              gamma TOROIDAL
C           .
            IF (     MTA( IHA ).EQ.IPOL       .AND.
     1               MTB( IHB ).EQ.ITOR       .AND.
     2               MTG( IHG ).EQ.2                ) THEN
C              .
C              . Test for C_ssq interaction
C              .
               REQVIT = 'CSSQ'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNT, IRNT
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
C              . Test for C_qss interaction
C              .
               REQVIT = 'CQSS'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNT, IRNT
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
C              . Test for C_sqs interaction
C              .
               REQVIT = 'CSQS'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNT, IRNT
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
            ENDIF
C           .
C           . N.B. there is no case alpha TOR, beta POL and 
C           . gamma POL. These are all identically zero.
C           .
C           . Fifth case: alpha TOROIDAL
C           .             beta  POLOIDAL
C           .             gamma TOROIDAL
C           .
            IF (     MTA( IHA ).EQ.IPOL       .AND.
     1               MTB( IHB ).EQ.ITOR       .AND.
     2               MTG( IHG ).EQ.2                ) THEN
C              .
C              . Test for C_ttq interaction
C              .
               REQVIT = 'CTTQ'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNT, IRNT
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
            ENDIF
C           .
C           . Sixth case: alpha TOROIDAL
C           .             beta  TOROIDAL
C           .             gamma POLOIDAL
C           .
            IF (     MTA( IHA ).EQ.ITOR       .AND.
     1               MTB( IHB ).EQ.ITOR       .AND.
     2               MTG( IHG ).EQ.1                ) THEN
C              .
C              . Test for C_tqt interaction
C              .
               REQVIT = 'CTQT'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNP, IRNP
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
            ENDIF
C           .
C           . Seventh case: alpha TOROIDAL
C           .               beta  TOROIDAL
C           .               gamma TOROIDAL
C           .
            IF (     MTA( IHA ).EQ.ITOR       .AND.
     1               MTB( IHB ).EQ.ITOR       .AND.
     2               MTG( IHG ).EQ.2                ) THEN
C              .
C              . Test for C_tqs interaction
C              .
               REQVIT = 'CTQS'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNT, IRNT
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
C              . Test for C_tsq interaction
C              .
               REQVIT = 'CTSQ'
               CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C              .
               IF ( OK ) THEN
               DO IR = ILNT, IRNT
                 RAD  = XARR( IR )
                 LOCG = INDFUN( IR, IHG, INARRG )
                 TERM = SDSVVT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                          IHB, INARRB, LG, V0A, V1A, V2A, V0B,
     2                          V1B, V2B, V3B )
                 VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
               ENDDO
               ENDIF
C              .
            ENDIF
C           .
 702      CONTINUE
          ENDDO
C         .
 701    CONTINUE
        ENDDO
C       .
 700  CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************
