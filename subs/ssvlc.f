C*********************************************************************
C subroutine Stored derivative Solution Vector Laplcian and Curl *****
C            -                 -        -      -            -    *****
C Steve Gibbons Tue Oct 26 12:09:14 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Adds a multiple of FAC* the curl of all Laplacian of all the       C
C velocity harmonics of the vector VI0 to the appropriate elements   C
C of the solution vector, VO.                                        C
C The first, second, third and fourth derivatives (precalculated by  C
C CASVDR) are stored in VI1, VI2, VI3 and VI4.                       C
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
C Taking the curl of a Laplacian of a toroidal harmonic with         C
C radial function x(r) results in a poloidal harmonic with scalar    C
C function D_l x(r).                                                 C
C                                                                    C
C Taking the curl of a Laplacian of a poloidal harmonic with         C
C radial function x(r) results in a toroidal harmonic with scalar    C
C function -D_l^2 x(r).                                              C
C                                                                    C
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
C     Note that INARRO, MHTO, MHLO and MHMO correspond to            C
C    exactly the above variables but for the output and not input    C
C    vectors. NRI must equal NRO and must both equal NR.             C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VI0       : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NRI*NHI                     C
C     VI1       : First  derivatives of VI0 vector. Dim ( * ).       C
C     VI2       : Second derivatives of VI0 vector. Dim ( * ).       C
C     VI3       : Third  derivatives of VI0 vector. Dim ( * ).       C
C     VI4       : Fourth derivatives of VI0 vector. Dim ( * ).       C
C     VO        : Solution vector. Dim ( * )  ( output )             C
C                 Length must be atleast NRO*NHO                     C
C     FAC       : Multiplier of curl to be added.                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SSVLC( NR, VI0, VI1, VI2, VI3, VI4, INARRI, MHTI,
     1                  MHLI, MHMI, VO, INARRO, MHTO, MHLO, MHMO,
     2                  FAC, ILNR, IRNR, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARRI( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        INARRO( * ), MHTO( * ), MHLO( * ), MHMO( * ),
     2        ILNR, IRNR
      DOUBLE PRECISION VI0( * ), VI1( * ), VI2( * ), VI3( * ),
     1                 VI4( * ), VO( * ), FAC, XARR( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IFORMI, IFORMO, NRI, NRO, NHI, NHO, IHI, IHO, IR,
     1        INDO, INDFUN, L, ICMPO, INDI
      DOUBLE PRECISION LOW, RAD, DL, D0F, D1F, D2F, D3F, D4F, DLDL
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
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
         PRINT *,' Subroutine SSVLC.'
         PRINT *,' IFORMI = ', IFORMI
         PRINT *,' IFORMO = ', IFORMO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( NRI.NE.NR .OR. NRO.NE.NR ) THEN
         PRINT *,' Subroutine SSVLC.'
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
        IF ( MHTI( IHI ).NE.1 .AND. MHTI( IHI ).NE.2 ) GOTO 50
        IF ( MHTI( IHI ).EQ.1 ) ICMPO = 2
        IF ( MHTI( IHI ).EQ.2 ) ICMPO = 1
C
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
             INDI = INDFUN( IR, IHI, INARRI )
C
             IF ( MHTI( IHI ).EQ.1 ) THEN
c              .
c              . we are doing poloidal --> toroidal
c              . so need to add -D_l^2 operator to RHS
c              . This requires we take up to 4^{th} derivatives
c              .
               L   = MHLI( IHI )
               RAD = XARR( IR )
               D0F = VI0( INDI )
               D2F = VI2( INDI )
               D3F = VI3( INDI )
               D4F = VI4( INDI )
               VO( INDO ) = VO( INDO ) 
     1           - FAC*DLDL( L, RAD, D0F, D2F, D3F, D4F )
c              .
             ELSE
c              .
c              . we are doing toroidal --> poloidal
c              . so we add D_l operator to RHS.
c              . need up to 2^{nd} derivatives
c              .
               L   = MHLI( IHI )
               RAD = XARR( IR )
               D0F = VI0( INDI )
               D1F = VI1( INDI )
               D2F = VI2( INDI )
               VO( INDO ) = VO( INDO )
     1           + FAC*DL( L, RAD, D0F, D1F, D2F )
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
