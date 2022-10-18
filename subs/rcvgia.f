C*********************************************************************
C subroutine Recalled Coefficient Velocity dot Grad Inhom. temp. Add *
C            -        -           -            -    -            -   *
C Steve Gibbons Tue Feb  1 16:33:26 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given a velocity which is stored in VECA (indexed by INARRA, MHTA, C
C MHLA, MHMA, MHPA), a temperature which is stored in VECB (indexed  C
C by INARRB, MHTB, MHLB, MHMB, MHPB) and the arrays calculated by    C
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
C     MPA      : Array length ( * ) - finite diff. scheme index. (A) C
C                                                                    C
C     INARRB   : Indexing array for beta  vector. Dim ( 3 )          C
C     MTB      : Array length ( * ) - atleast length NHB             C
C                  See above for key. (corresponds to beta vector)   C
C     MLB      : Array length ( * ) - atleast length NHB             C
C                  Sph. harm. degree, l.                             C
C     MMB      : Array length ( * ) - atleast length NHB             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C     MPB      : Array length ( * ) - finite diff. scheme index. (B) C
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
C     NBN       : Number of bounding nodes. See ASVDR.               C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C     NCFM      : Leading dimension of SVFDC. At least (2*NBN+1)     C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     MIB      : Dim ( * ) - atleast atleast length NHB.             C
C                For each temperature harmonic, IHB, MIB(IHB) gives  C
C                the index of the array CAFIT which stores the       C
C                coefficients for that radial function.              C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C                                                                    C
C                If MIB( IH ) = IITH, then CA, CB and CC are         C
C                respectively stored in CAFIT( 1, IITH ),            C
C                CAFIT( 2, IITH ) and CAFIT( 3, IITH ).              C
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
C     VECA      : Solution vector for alpha function. Dim ( * ).     C
C     VECB      : Solution vector for beta  function. Dim ( * ).     C
C     VECG      : Solution vector for gamma function. Dim ( * ).     C
C     XARR      : Dim ( NR ). Radial grid spacings.                  C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C     FAC       : Coefficient of term to be added to VECG.           C
C     CAFIT     : Dimension ( 3, * ). See MIB.                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RCVGIA( NVI, IHNALP, IHNBET, IHNGAM, INARRA, MTA,
     1          MLA, MMA, MPA, INARRB, MTB, MLB, MMB, MPB, INARRG,
     2          MTG, MLG, MMG, ILN, IRN, NR, NBN, NDRVS, NDRVM, NCFM,
     3          NDCS, TVHI, CVI, VECA, VECB, VECG, XARR, SVFDC, FAC,
     4          MIB, CAFIT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NVI, IHNALP( * ), IHNBET( * ), IHNGAM( * ), INARRA(3),
     1        MTA( * ), MLA( * ), MMA( * ), MPA( * ), INARRB(3),
     2        MTB( * ), MLB( * ), MMB( * ), MPB( * ), INARRG(3),
     3        MTG( * ), MLG( * ), MMG( * ), ILN, IRN, NR, NBN
      INTEGER NDRVS, NDRVM, NCFM, NDCS, MIB( * )
      CHARACTER *(4) TVHI( * )
      DOUBLE PRECISION CVI( * ), VECA( * ), VECB( * ), VECG( * ),
     1                 XARR( NR ), SVFDC( NCFM, NR, NDRVM+1, NDCS ),
     2                 FAC, CAFIT( 3, * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER LA, LB, IR, ISA, ISB, IHA, IHB, IHG, INDA, INDB, INDG,
     1        LOCG, INDFUN, ISHCIA, NRRA, NRRB, NRRG, NHA, NHB, NHG,
     2        IITH
      DOUBLE PRECISION RAD, COEF, TERM, SHNVGI, RI, RO, CA, CB, CC,
     1                 LOW
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
C Just check the number of grid nodes in A, B, and G
C
      NRRA = INARRA( 2 )
      NRRB = INARRB( 2 )
      NRRG = INARRG( 2 )
      IF ( NRRA.NE.NR .OR. NRRB.NE.NR .OR. NRRG.NE.NR ) THEN
        PRINT *,' Subroutine RCVGIA. NR = ', NR
        PRINT *,' NRRA = ',NRRA,' NRRB = ',NRRB,' NRRG = ',NRRG
        PRINT *,' Program aborted,'
        STOP
      ENDIF
C
      RI = XARR(  1 )
      RO = XARR( NR )
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
        ISA  = MPA( IHA )
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
          ISB  = MPB( IHB )
          IITH = MIB( IHB )
C         .
          CA     = CAFIT( 1, IITH )
          CB     = CAFIT( 2, IITH )
          CC     = CAFIT( 3, IITH )
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
              TERM = SHNVGI( IR, RAD, NBN, REQVIT, NR, NDRVS, NDRVM,
     1                 NCFM, NDCS, SVFDC, LA, VECA, IHA, ISA, INARRA,
     2                 LB, VECB, IHB, ISB, INARRB,
     3                 RI, RO, CA, CB, CC )
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
              TERM = SHNVGI( IR, RAD, NBN, REQVIT, NR, NDRVS, NDRVM,
     1                 NCFM, NDCS, SVFDC, LA, VECA, IHA, ISA, INARRA,
     2                 LB, VECB, IHB, ISB, INARRB,
     3                 RI, RO, CA, CB, CC )
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
              TERM = SHNVGI( IR, RAD, NBN, REQVIT, NR, NDRVS, NDRVM,
     1                 NCFM, NDCS, SVFDC, LA, VECA, IHA, ISA, INARRA,
     2                 LB, VECB, IHB, ISB, INARRB,
     3                 RI, RO, CA, CB, CC )
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

