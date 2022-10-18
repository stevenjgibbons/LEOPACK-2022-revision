C*********************************************************************
C subroutine Stored derivative Recalled coefficient Velocity dot  ****
C            -                 -                    -             ****
C                                                  Grad Theta Add ****
C                                                  -    -     -   ****
C Steve Gibbons Mon Feb  7 09:27:27 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given a velocity which is stored in V0A (indexed by INARRA, MHTA,  C
C MHLA, MHMA - first derivatives in V1A), a temperature which is     C
C stored in V0B (indexed by INARRB, MHTB, MHLB, MHMB - first         C
C derivatives in V1B) and the arrays calculated by                   C
C VSPCC (i.e. IHNALP, IHNBET, IHNGAM, TVHI, CVI) to calculate the    C
C v . Grad ( theta ) terms and add a multiple FAC of these terms to  C
C the appropriate elements of VECG (indexed by INARRG, MHTG, MHLG,   C
C MHMG).                                                             C
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
C                  (Output by subroutine VSPCC).                     C
C     IHNBET   : Number of beta harmonics. Dim ( * )                 C
C                  (Output by subroutine VSPCC).                     C
C     IHNGAM   : Number of gamma harmonics. Dim ( * )                C
C                  (Output by subroutine VSPCC).                     C
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
C                                                                    C
C     INARRB   : Indexing array for beta  vector. Dim ( 3 )          C
C     MTB      : Array length ( * ) - atleast length NHB             C
C                  See above for key. (corresponds to beta vector)   C
C     MLB      : Array length ( * ) - atleast length NHB             C
C                  Sph. harm. degree, l.                             C
C     MMB      : Array length ( * ) - atleast length NHB             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     INARRG   : Indexing array for gamma vector. Dim ( 3 )          C
C     MTG      : Array length ( * ) - atleast length NHG             C
C                  See above for key. (corresponds to gamma vec.)    C
C     MLG      : Array length ( * ) - atleast length NHG             C
C                  Sph. harm. degree, l.                             C
C     MMG      : Array length ( * ) - atleast length NHG             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     ILN      : Left most node on which operation is to be          C
C                 performed.                                         C
C     IRN      : Right most node on which operation is to be         C
C                 performed.                                         C
C                                                                    C
C     NR       : Number of radial grid nodes.                        C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. ( * ).       C
C                 = 'SPQQ', 'SPSS' etc. according to the corresp.    C
C                 vector interaction.                                C
C                  (Output by subroutine VSPCC).                     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     CVI       : Coefficient of vector interaction. Dim ( * ).      C
C                  (Output by subroutine VSPCC).                     C
C     V0A       : Solution vector for alpha function. Dim ( * ).     C
C     V1A       : First derivatives for alpha function. Dim ( * ).   C
C     V0B       : Solution vector for beta  function. Dim ( * ).     C
C     V1B       : First derivatives for beta  function. Dim ( * ).   C
C     VECG      : Solution vector for gamma function. Dim ( * ).     C
C     XARR      : Dim ( NR ). Radial grid spacings.                  C
C     FAC       : Coefficient of term to be added to VECG.           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SRVGTA( NVI, IHNALP, IHNBET, IHNGAM, INARRA, MTA,
     1          MLA, MMA, INARRB, MTB, MLB, MMB, INARRG,
     2          MTG, MLG, MMG, ILN, IRN, NR, TVHI, CVI, V0A, V1A,
     3          V0B, V1B, VECG, XARR, FAC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NVI, IHNALP( * ), IHNBET( * ), IHNGAM( * ), INARRA(3),
     1        MTA( * ), MLA( * ), MMA( * ), INARRB(3),
     2        MTB( * ), MLB( * ), MMB( * ), INARRG(3),
     3        MTG( * ), MLG( * ), MMG( * ), ILN, IRN, NR
      CHARACTER *(4) TVHI( * )
      DOUBLE PRECISION CVI( * ), V0A( * ), V1A( * ), V0B( * ),
     1                 V1B( * ), VECG( * ), XARR( NR ), FAC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER LA, LB, IR, IHA, IHB, IHG, INDA, INDB, INDG,
     1        LOCG, INDFUN, ISHCIA, NRRA, NRRB, NRRG, NHA, NHB, NHG
      DOUBLE PRECISION RAD, COEF, TERM, SDSVGT, LOW
      CHARACTER *(4) REQVIT
      PARAMETER ( LOW = 1.0d-9 )
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
C     . Early exit for zero value of FAC
C     .
      IF ( DABS( FAC ).LT.LOW ) RETURN
C
C Just check the number of grid nodes in A, B, and G
C
      NRRA = INARRA( 2 )
      NRRB = INARRB( 2 )
      NRRG = INARRG( 2 )
      IF ( NRRA.NE.NR .OR. NRRB.NE.NR .OR. NRRG.NE.NR ) THEN
        PRINT *,' Subroutine SRVGTA. NR = ', NR
        PRINT *,' NRRA = ',NRRA,' NRRB = ',NRRB,' NRRG = ',NRRG
        PRINT *,' Program aborted,'
        STOP
      ENDIF
C
      NHA = INARRA( 3 )
      NHB = INARRB( 3 )
      NHG = INARRG( 3 )
C
C Begin loop around ALPHA velocity harmonics
C
      DO IHA = 1, NHA
C       .
C       . Jump to 700 (end of loop) if A harmonic is not velocity
C       .
        IF ( MTA( IHA ).NE.1 .AND. MTA( IHA ).NE.2 ) GOTO 700
C       .
        INDA = ISHCIA( IHA, MLA, MMA )
        LA   = MLA( IHA )
C       .
C       . Now loop around BETA temperature harmonics
C       .
        DO IHB = 1, NHB
C         .
C         . Jump to 701 (end of loop) if B harmonic is not temp.
C         .
          IF ( MTB( IHB ).NE.3 ) GOTO 701
C         .
          INDB = ISHCIA( IHB, MLB, MMB )
          LB   = MLB( IHB )
C         .
C         . Now loop around GAMMA harmonics (temp.)
C         .
          DO IHG = 1, NHG
C           .
C           . Jump to 702 (end of loop) if G harm. not temp.
C           .
            IF ( MTG( IHG ).NE.3 ) GOTO 702
C           .
            INDG = ISHCIA( IHG, MLG, MMG )
C           .
C           . So our alpha harmonic is velocity
C           . Our beta harmonic is temperature
C           . Our gamma harmonic is temperature
C           .
            IF ( MTA( IHA ).EQ.2 ) GOTO 703
C           .
C           . Case alpha is poloidal
C           . First test for S_qq interaction
C           .
            REQVIT = 'SPQQ'
            CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C           .
            IF ( .NOT. OK ) GOTO 704
C           .
            DO IR = ILN, IRN
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSVGT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                 IHB, INARRB, V0A, V1A, V0B, V1B )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
C           .
 704        CONTINUE
C           .
C           . Test for S_ss interaction
C           .
            REQVIT = 'SPSS'
            CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C           .
            IF ( .NOT. OK ) GOTO 703
C           .
            DO IR = ILN, IRN
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSVGT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                 IHB, INARRB, V0A, V1A, V0B, V1B )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
C           .
 703        CONTINUE
C           .
C           . Case alpha is toroidal
C           . Test for S_ts interaction
C           .
            REQVIT = 'SPTS'
            CALL VICEXR( NVI, IHNALP, IHNBET, IHNGAM, INDA,
     1                INDB, INDG, TVHI, REQVIT, CVI, COEF, OK )
C           .
            IF ( .NOT. OK ) GOTO 702
C           .
            DO IR = ILN, IRN
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSVGT( IR, RAD, REQVIT, LA, IHA, INARRA, LB,
     1                 IHB, INARRB, V0A, V1A, V0B, V1B )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
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

