C*********************************************************************
C subroutine Stored derivative Recalled coefficient curl (K Cross  ***
C            -                 -                          - -      ***
C                                                    Velocity) Add ***
C                                                    -         -   ***
C Steve Gibbons Mon Feb  7 10:22:05 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given a velocity which is stored in V0A (indexed by INARRA, MTA,   C
C MLA, MMA - 1st and 2nd deriv.s in V1A and V2A resp.) and the       C
C arrays calculated by CFVICC (i.e. IHNALP,                          C
C IHNGAM, TVHI, CVI) to calculate the curl ( k x v ) terms and add   C
C a multiple FAC of these terms to the appropriate elements of VECG  C
C (indexed by INARRG, MTG, MLG, MMG).                                C
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
C                  (Output by subroutine CFVICC).                    C
C     IHNGAM   : Number of gamma harmonics. Dim ( * )                C
C                  (Output by subroutine CFVICC).                    C
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
C                 = 'KCTS', 'KCTT' etc. according to the corresp.    C
C                 vector interaction.                                C
C                  (Output by subroutine CFVICC).                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     CVI       : Coefficient of vector interaction. Dim ( * ).      C
C                  (Output by subroutine CFVICC).                    C
C     V0A       : 0 derivatives of alpha velocity. Dim ( * )         C
C     V1A       : 1 derivatives of alpha velocity. Dim ( * )         C
C     V2A       : 2 derivatives of alpha velocity. Dim ( * )         C
C     XARR      : Dim ( NR ). Radial grid spacings.                  C
C     FAC       : Coefficient of term to be added to VECG.           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SRKCVA( NVI, IHNALP, IHNGAM, INARRA, MTA, MLA, MMA,
     1                 INARRG, MTG, MLG, MMG, ILNP, IRNP, ILNT, IRNT,
     2                 NR, TVHI, CVI, V0A, V1A, V2A, XARR, FAC, VECG)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NVI, IHNALP( * ), IHNGAM( * ), INARRA(3),
     1        MTA( * ), MLA( * ), MMA( * ), INARRG(3), NR,
     2        MTG( * ), MLG( * ), MMG( * ), ILNP, IRNP, ILNT, IRNT
      CHARACTER *(4) TVHI( * )
      DOUBLE PRECISION CVI( * ), V0A( * ), V1A( * ), V2A( * ), 
     1                 VECG( * ), XARR( NR ), FAC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER LA, IR, IHA, IHG, INDA, INDG,
     1        LOCG, INDFUN, ISHCIA, NRRA, NRRG, NHA, NHG, LG
      DOUBLE PRECISION RAD, COEF, TERM, SDSKCV, LOW
      CHARACTER *(4) REQVIT
      LOGICAL OK
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
C     . Early exit for zero value of FAC
C     .
      IF ( DABS( FAC ).LT.LOW ) RETURN
C
C Just check the number of grid nodes in A and G
C
      NRRA = INARRA( 2 )
      NRRG = INARRG( 2 )
      IF ( NRRA.NE.NR .OR. NRRG.NE.NR ) THEN
        PRINT *,' Subroutine SRKCVA. NR = ', NR
        PRINT *,' NRRA = ',NRRA,' NRRG = ',NRRG
        PRINT *,' Program aborted,'
        STOP
      ENDIF
C
      NHA = INARRA( 3 )
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
C       . Now loop around GAMMA harmonics (vorticity).
C       .
        DO IHG = 1, NHG
C         .
C         . Jump to 702 (end of loop) if G harm. not pol. or tor.
C         .
          IF ( MTG( IHG ).NE.1 .AND. MTG( IHG ).NE.2 ) GOTO 702
C         .
          INDG = ISHCIA( IHG, MLG, MMG )
          LG   = MLG( IHG )
C         .
C         . So our alpha harmonic is velocity
C         . Our gamma harmonic is vorticity
C         .
C         . First case: alpha poloidal, gamma poloidal
C         .
          IF ( MTA( IHA ).EQ.1 .AND. MTG( IHG ).EQ.1 ) THEN
C           .
C           . First test for K_qt interaction.
C           .
            REQVIT = 'KCQT'
            CALL VLICER( NVI, IHNALP, IHNGAM, INDA, INDG, TVHI,
     1                   REQVIT, CVI, COEF, OK )
C           .
            IF ( OK ) THEN
            DO IR = ILNP, IRNP
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSKCV( IR, RAD, REQVIT, LA, IHA, INARRA, LG,
     1                 V0A, V1A, V2A )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
            ENDIF
C           .
C           . Test for K_st interaction.
C           .
            REQVIT = 'KCST'
            CALL VLICER( NVI, IHNALP, IHNGAM, INDA, INDG, TVHI,
     1                   REQVIT, CVI, COEF, OK )
C           .
            IF ( OK ) THEN
            DO IR = ILNP, IRNP
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSKCV( IR, RAD, REQVIT, LA, IHA, INARRA, LG,
     1                 V0A, V1A, V2A )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
            ENDIF
C           .
          ENDIF
C         .
C         . Second case: alpha poloidal, gamma toroidal
C         .
          IF ( MTA( IHA ).EQ.1 .AND. MTG( IHG ).EQ.2 ) THEN
C           .
C           . Test for K_sq interaction.
C           .
            REQVIT = 'KCSQ'
            CALL VLICER( NVI, IHNALP, IHNGAM, INDA, INDG, TVHI,
     1                   REQVIT, CVI, COEF, OK )
C           .
            IF ( OK ) THEN
            DO IR = ILNT, IRNT
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSKCV( IR, RAD, REQVIT, LA, IHA, INARRA, LG,
     1                 V0A, V1A, V2A )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
            ENDIF
C           .
C           . Test for K_ss interaction.
C           .
            REQVIT = 'KCSS'
            CALL VLICER( NVI, IHNALP, IHNGAM, INDA, INDG, TVHI,
     1                   REQVIT, CVI, COEF, OK )
C           .
            IF ( OK ) THEN
            DO IR = ILNT, IRNT
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSKCV( IR, RAD, REQVIT, LA, IHA, INARRA, LG,
     1                 V0A, V1A, V2A )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
            ENDIF
C           .
C           . Test for K_qs interaction.
C           .
            REQVIT = 'KCQS'
            CALL VLICER( NVI, IHNALP, IHNGAM, INDA, INDG, TVHI,
     1                   REQVIT, CVI, COEF, OK )
C           .
            IF ( OK ) THEN
            DO IR = ILNT, IRNT
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSKCV( IR, RAD, REQVIT, LA, IHA, INARRA, LG,
     1                 V0A, V1A, V2A )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
            ENDIF
C           .
          ENDIF
C         .
C         . Third case: alpha toroidal, gamma poloidal
C         .
          IF ( MTA( IHA ).EQ.2 .AND. MTG( IHG ).EQ.1 ) THEN
C           .
C           . First test for K_tt interaction.
C           .
            REQVIT = 'KCTT'
            CALL VLICER( NVI, IHNALP, IHNGAM, INDA, INDG, TVHI,
     1                   REQVIT, CVI, COEF, OK )
C           .
            IF ( OK ) THEN
            DO IR = ILNP, IRNP
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSKCV( IR, RAD, REQVIT, LA, IHA, INARRA, LG,
     1                 V0A, V1A, V2A )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
            ENDIF
C           .
          ENDIF
C         .
C         . Fourth case: alpha toroidal, gamma toroidal
C         .
          IF ( MTA( IHA ).EQ.2 .AND. MTG( IHG ).EQ.2 ) THEN
C           .
C           . First test for K_ts interaction.
C           .
            REQVIT = 'KCTS'
            CALL VLICER( NVI, IHNALP, IHNGAM, INDA, INDG, TVHI,
     1                   REQVIT, CVI, COEF, OK )
C           .
            IF ( OK ) THEN
            DO IR = ILNT, IRNT
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSKCV( IR, RAD, REQVIT, LA, IHA, INARRA, LG,
     1                 V0A, V1A, V2A )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
            ENDIF
C           .
C           . First test for K_tq interaction.
C           .
            REQVIT = 'KCTQ'
            CALL VLICER( NVI, IHNALP, IHNGAM, INDA, INDG, TVHI,
     1                   REQVIT, CVI, COEF, OK )
C           .
            IF ( OK ) THEN
            DO IR = ILNT, IRNT
              RAD  = XARR( IR )
              LOCG = INDFUN( IR, IHG, INARRG )
              TERM = SDSKCV( IR, RAD, REQVIT, LA, IHA, INARRA, LG,
     1                 V0A, V1A, V2A )
              VECG( LOCG ) = VECG( LOCG ) + FAC*TERM*COEF
            ENDDO
            ENDIF
C           .
          ENDIF
C         .
 702    CONTINUE
        ENDDO
C       .
 700  CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************

