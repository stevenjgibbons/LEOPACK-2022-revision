C*********************************************************************
C subroutine Stored derivative Solution Vector LaPlacian *************
C            -                 -        -      - -       *************
C Steve Gibbons Mon Feb  7 08:40:26 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Adds a multiple of FAC* the laplacian of all harmonics IHI in      C
C solution vector VI0 which satisfy MHTI( IHI ) = ICMP to the        C
C corresponding harmonics IHO in vector VO with MHTO( IHO ) = ICMP.  C
C                                                                    C
C The first and second derivatives (precalculated by CASVDR) are     C
C stored in VI1 and VI2.                                             C
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
C The lapacian of a scalar with radial function f(r)                 C
C is a scalar with radial function D_l f(r) - l is degree of S.H.    C
C                                                                    C
C The lapacian of a toroidal harmonic with radial function f(r)      C
C is toroidal with radial function D_l f(r) - l is degree of S.H.    C
C                                                                    C
C The lapacian of a poloidal harmonic with radial function f(r)      C
C is poloidal with radial function D_l f(r) - l is degree of S.H.    C
C                                                                    C
C  This routine is therefore valid for all forms of harmonics        C
C (itype = 1, 2, 3, 4 and 5) and the result must go in the same      C
C type of harmonic as the original.                                  C
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
C     ICMP     : The type of harmonic we wish to take curl of.       C
C                                                 see above.         C
C                                                                    C
C     Note that INARRO, MHTO, MHLO and MHMO correspond to            C
C    exactly the above variables but for the output and not input    C
C    vectors. NRI must equal NRO and must both equal NR.             C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VI0       : Solution vector. Dim ( * ) ( input ). 0 derivs.    C
C                 Length must be atleast NRI*NHI                     C
C     VI1       : Solution vector. Dim ( * ) ( input ). 1 derivs.    C
C                 Length must be atleast NRI*NHI                     C
C     VI2       : Solution vector. Dim ( * ) ( input ). 2 derivs.    C
C                 Length must be atleast NRI*NHI                     C
C     VO        : Solution vector. Dim ( * )  ( output )             C
C                 Length must be atleast NRO*NHO                     C
C     FAC       : Multiplier of curl to be added.                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SSVLP( NR, VI0, VI1, VI2, INARRI, MHTI, MHLI, MHMI,
     1                  ICMP, VO, INARRO, MHTO, MHLO, MHMO, FAC,
     2                  ILNR, IRNR, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARRI( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        ICMP, INARRO( * ), MHTO( * ), MHLO( * ),
     2        MHMO( * ), ILNR, IRNR
      DOUBLE PRECISION VI0( * ), VI1( * ), VI2( * ), VO( * ), FAC,
     1                 XARR( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IFORMI, IFORMO, NRI, NRO, NHI, NHO, IHI, IHO, IR,
     1        INDO, INDI, INDFUN, L
      DOUBLE PRECISION LOW, RAD, DL, D0F, D1F, D2F
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     .
C     . However do check ICMP
C     .
C
      IF ( ICMP.NE.1 .AND. ICMP.NE.2 .AND. ICMP.NE.3 .AND.
     1        ICMP.NE.4 .AND. ICMP.NE.5           )    THEN
         PRINT *,' Subroutine SSVLP.'
         PRINT *,' ICMP = ', ICMP
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
         PRINT *,' Subroutine SSVLP.'
         PRINT *,' IFORMI = ', IFORMI
         PRINT *,' IFORMO = ', IFORMO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( NRI.NE.NR .OR. NRO.NE.NR ) THEN
         PRINT *,' Subroutine SSVLP.'
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
        IF ( MHTI( IHI ).NE.ICMP ) GOTO 50
C
C ok - so our harmonic is of the correct type
C so let's look for its corresponding harmonic
C in VO
C
        DO IHO = 1, NHO
          IF (   MHTO( IHO ).EQ.ICMP          .AND.
     1           MHLO( IHO ).EQ.MHLI( IHI )   .AND.
     2           MHMO( IHO ).EQ.MHMI( IHI )         ) THEN
C
C o.k. we've found the corresponding harmonic
C
            DO IR = ILNR, IRNR
C
C Find locations of source and destination vectors ...
C
             INDO = INDFUN( IR, IHO, INARRO )
             INDI = INDFUN( IR, IHI, INARRI )
C            .
             L   = MHLI( IHI )
             RAD = XARR( IR )
             D0F = VI0( INDI )
             D1F = VI1( INDI )
             D2F = VI2( INDI )
             VO( INDO ) = VO( INDO )
     1                       + FAC*DL( L, RAD, D0F, D1F, D2F )
C            .
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
