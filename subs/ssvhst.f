C*********************************************************************
C subroutine Stored derivative Solution Vector Heat Source Term ******
C            -                 -        -      -    -      -    ******
C Steve Gibbons Thu Nov 18 16:52:11 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C If \nabla^2 T_0( r ) = C (C is a constant), then T_0 has the       C
C general solution                                                   C
C                                                                    C
C            CB1 r^2                                                 C
C T_0( r ) = ------- - CB2 r^{-1}      (cb1 and cb2 are constants)   C
C               2                                                    C
C                                                                    C
C                                                                    C
C and so                                                             C
C                                                                    C
C              [                                     ]               C
C \nabla T_0 = [  CB1 r  +  CB2 r^{-2}     0     0   ]               C
C              [                       ,      ,      ]               C
C                                                                    C
C                                                                    C
C                                          [         CB2  ]          C
C and so v . \nabla T_0 = l(l+1)P(r) Y_l^m [ CB1 +  ----- ]          C
C                                          [         r^3  ]          C
C                                                                    C
C SSVHST adds this term to the appropriate temperature radial funcs. C
C                                                                    C
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
C     INARRI    : Int. parameter array corresponding to VI0.         C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VI0       : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NRI*NHI                     C
C     VO        : Solution vector. Dim ( * )  ( output )             C
C                 Length must be atleast NRO*NHO                     C
C     FAC       : Multiplier of term to be added.                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C     CB1       : Constant - see above.                              C
C     CB2       : Constant - see above.                              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SSVHST( NR, VI0, INARRI, MHTI, MHLI, MHMI, VO,
     1                   INARRO, MHTO, MHLO, MHMO, FAC, ILNR, IRNR,
     2                   XARR, CB1, CB2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARRI( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        INARRO( * ), MHTO( * ), MHLO( * ),
     2        MHMO( * ), ILNR, IRNR
      DOUBLE PRECISION VI0( * ), VO( * ), FAC, XARR( NR ), CB1, CB2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NRI, NRO, NHI, NHO, IR, L, IHI, IHO,
     1        INDFUN, INDO, INDI
      DOUBLE PRECISION LOW, RAD, D0F, F1, R3
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NRI    = INARRI( 2 )
      NHI    = INARRI( 3 )
C
      NRO    = INARRO( 2 )
      NHO    = INARRO( 3 )
C
      IF ( NRI.NE.NR .OR. NRO.NE.NR ) THEN
         PRINT *,' Subroutine SSVHST.'
         PRINT *,' NR   = ', NR
         PRINT *,' NRI  = ', NRI
         PRINT *,' NRO  = ', NRO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C     . early exit?
C
      IF ( ABS(FAC).LT.LOW ) RETURN
C     .
C     . Loop around in harmonics
C     .
      DO IHI = 1, NHI
        IF ( MHTI( IHI ).NE.1 ) GOTO 50
C
C ok - so our harmonic is poloidal velocity
C so let's look for the corresponding temperature
C harmonic in VO
C
        DO IHO = 1, NHO
          IF (   MHTO( IHO ).EQ.3              .AND.
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
              L  = MHLI( IHI )
              RAD = XARR( IR )
              IF ( ABS( RAD ).LT.LOW ) THEN
                PRINT *,' Subroutine SSVHST.'
                PRINT *,' Rad at node ',IR,' is ',RAD
                PRINT *,' Division by zero imminent.'
                PRINT *,' Program aborted.'
                STOP
              ENDIF
              R3  = RAD*RAD*RAD
              D0F = VI0( INDI )
              F1  = D0F*DBLE( L )*DBLE( L + 1 )
              VO( INDO ) = VO( INDO ) + FAC*( CB1 + CB2/R3 )*F1
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
