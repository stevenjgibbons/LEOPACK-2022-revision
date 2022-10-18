C*********************************************************************
C subroutine Stored derivative Solution Vector CurL ******************
C            -                 -        -      -  - ******************
C Steve Gibbons Tue Oct 26 10:12:07 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Adds a multiple of FAC* the curl of all harmonics IHI in solution  C
C vector VI0 (first and second derivatives in VI1 and VI2 resp.)     C
C which satisfy MHTI( IHI ) = ICMPI to the corresponding             C
C harmonics IHO in vector VO with MHTO( IHO ) = ICMPO.               C
C                                                                    C
C Now MHTI defines what each scalar function in a solution vector    C
C represents.                                                        C
C                                                                    C
C         MHTI( IH ) = 1 --> harmonic is poloidal velocity           C
C         MHTI( IH ) = 2 --> harmonic is toroidal velocity           C
C         MHTI( IH ) = 3 --> harmonic is temperature.                C
C         MHTI( IH ) = 4 --> harmonic is poloidal magnetic field     C
C         MHTI( IH ) = 5 --> harmonic is toroidal magnetic field     C
C                                                                    C
C and similarly MHTO for the resulting vector.                       C
C                                                                    C
C The program aborts if ICMPI = 3 as curl has no meaning.            C
C                                                                    C
C Taking the curl of a toroidal harmonic with radial function        C
C x( r ) results in a poloidal harmonic with r.f. x(r)               C
C                                                                    C
C Hence ICMPI = 2(5) would normally correspond to ICMPO = 1(4)       C
C                                                                    C
C                                                                    C
C Taking the curl of a poloidal harmonic with radial function        C
C x( r ) results in a toroidal harmonic with r.f. -D_l x(r)          C
C                                                                    C
C Hence ICMPI = 1(4) would normally correspond to ICMPO = 2(5)       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARRI    : Int. parameter array corresponding to VI.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARRI( 1 ) = IFORMF                               C
C                 INARRI( 2 ) = NRI     See INDFUN for details       C
C                 INARRI( 3 ) = NHI                                  C
C                                                                    C
C     MHTI      : Array length ( * ) - atleast length NHI            C
C                  See above for key. (corresponds to input vec.)    C
C     MHLI      : Array length ( * ) - atleast length NHI            C
C                  Sph. harm. degree, l.                             C
C     MHMI      : Array length ( * ) - atleast length NHI            C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     ICMPI     : The type of harmonic we wish to take curl of.      C
C                                                 see above.         C
C                                                                    C
C     Note that INARRO, MHTO, MHLO, MHMO and ICMPO correspond to     C
C    exactly the above variables but for the output and not input    C
C    vectors. ICMPI cannot be 3. NRI must equal NRO and must both    C
C    equal NR.                                                       C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VI0       : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NRI*NHI                     C
C     VI1       : First derivs of VI0. Dim ( * ) ( input )           C
C     VI2       : Second derivs of VI0. Dim ( * ) ( input )          C
C     VO        : Solution vector. Dim ( * )  ( output )             C
C                 Length must be atleast NRO*NHO                     C
C     FAC       : Multiplier of curl to be added.                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SSVCL( NR, VI0, VI1, VI2, INARRI, MHTI, MHLI, MHMI,
     1                  ICMPI, VO, INARRO, MHTO, MHLO, MHMO, ICMPO,
     2                  FAC, ILNR, IRNR, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARRI( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        ICMPI, INARRO( * ), MHTO( * ), MHLO( * ),
     2        MHMO( * ), ICMPO, ILNR, IRNR
      DOUBLE PRECISION VI0( * ), VI1( * ), VI2( * ), VO( * ),
     1                 FAC, XARR( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IFORMI, IFORMO, NRI, NRO, NHI, NHO, IHI, IHO, IR,
     1        INDO, INDFUN, L, INDI
      DOUBLE PRECISION LOW, RAD, DL, D0F, D1F, D2F
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     .
C     . Check ICMPI
C     .
C
      IF ( ICMPI.NE.1 .AND. ICMPI.NE.2 .AND. ICMPI.NE.4
     1                .AND. ICMPI.NE.5           )    THEN
         PRINT *,' Subroutine SSVCL.'
         PRINT *,' ICMPI = ', ICMPI
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IFORMI = INARRI( 1 )
      NRI    = INARRI( 2 )
      NHI    = INARRI( 3 )
C
      IFORMO = INARRO( 1 )
      NRO    = INARRO( 2 )
      NHO    = INARRO( 3 )
C
      IF ( IFORMI.EQ.1 .OR. IFORMI.EQ.2 .OR.
     1     IFORMO.EQ.1 .OR. IFORMO.EQ.2 ) THEN
         PRINT *,' Subroutine SSVCL.'
         PRINT *,' IFORMI = ', IFORMI
         PRINT *,' IFORMO = ', IFORMO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( NRI.NE.NR .OR. NRO.NE.NR ) THEN
         PRINT *,' Subroutine SSVCL.'
         PRINT *,' NR   = ', NR
         PRINT *,' NRI  = ', NRI
         PRINT *,' NRO  = ', NRO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C     . early exit?
C
      IF ( DABS(FAC).LT.LOW ) RETURN
C     .
C     . Loop around in harmonics
C     .
      DO IHI = 1, NHI
        IF ( MHTI( IHI ).NE.ICMPI ) GOTO 50
C
C ok - so our harmonic is of the correct type
C so let's look for its corresponding harmonic
C in VO
C
        DO IHO = 1, NHO
          IF (   MHTO( IHO ).EQ.ICMPO          .AND.
     1           MHLO( IHO ).EQ.MHLI( IHI )    .AND.
     2           MHMO( IHO ).EQ.MHMI( IHI )         ) THEN
C
C o.k. we've found the corresponding harmonic
C
            DO IR = ILNR, IRNR
C
C Find locations of source and destination vectors ...
C
             INDO = INDFUN( IR, IHO, INARRO )
C
             IF ( MHTI( IHI ).EQ.1 .OR. MHTI( IHI ).EQ.4 ) THEN
c              .
c              . we are doing poloidal --> toroidal
c              . so need to add -D_l operator to RHS
c              . This requires we take derivatives
c              .
               L  = MHLI( IHI )
               RAD = XARR( IR )
               INDI = INDFUN( IR, IHI, INARRI )
               D0F = VI0( INDI )
               D1F = VI1( INDI )
               D2F = VI2( INDI )
               VO( INDO ) = VO( INDO ) 
     1           - FAC*DL( L, RAD, D0F, D1F, D2F )
c              .
             ELSE
c              .
c              . we are doing toroidal --> poloidal
c              . so we just add on the value of the term
c              .
               INDI = INDFUN( IR, IHI, INARRI )
               D0F = VI0( INDI )
               VO( INDO ) = VO( INDO ) + FAC*D0F
c              .
             ENDIF
C
            ENDDO
C
          ENDIF 
        ENDDO
C
 50   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
