C*********************************************************************
C subroutine Non-uniform Grid Solution Vector LaPlacian **************
C            -           -    -        -      - -       **************
C Steve Gibbons Mon Sep 27 13:09:50 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Adds a multiple of FAC* the laplacian of all harmonics IHI in      C
C solution vector VI which satisfy MHTI( IHI ) = ICMP to the         C
C corresponding harmonics IHO in vector VO with MHTO( IHO ) = ICMP.  C
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
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NDRVS     : Highest derivative stored in FDCM.                 C
C                 Must be atleast 2.                                 C
C                                                                    C
C     NFDCM     : Leading dim of FDCM. See FDCMBD.                   C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C     ILNC      : First radial node which may be used for f.d. form. C
C     IRNC      : Last radial node which may be used for f.d. form.  C
C                                                                    C
C  ILNC and IRNC correspond to NLMC and NRMC in the call to FDCMBD.  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VI        : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NRI*NHI                     C
C     VO        : Solution vector. Dim ( * )  ( output )             C
C                 Length must be atleast NRO*NHO                     C
C     FAC       : Multiplier of curl to be added.                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     FDCM      : Dim ( NFDCM, NR, NDRVS ). Finite diff. coeff.s     C
C                (Not needed if we are doing tor --> pol).           C
C                Must be formed in advance by a call to FDCMBD.      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NGSVLP( NR, VI, INARRI, MHTI, MHLI, MHMI, ICMP,
     1                   VO, INARRO, MHTO, MHLO, MHMO, FAC,
     2                   NBN, NDRVS, ILNR, IRNR, ILNC, IRNC, FDCM,
     3                   XARR, NFDCM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARRI( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        ICMP, INARRO( * ), MHTO( * ), MHLO( * ), MHMO( * ),
     2        NBN, NDRVS, ILNR, IRNR, ILNC, IRNC, NFDCM
      DOUBLE PRECISION VI( * ), VO( * ), FAC, XARR( NR ),
     1                 FDCM( NFDCM, NR, NDRVS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IFORMI, IFORMO, NRI, NRO, NHI, NHO, IHI, IHO, IR,
     1        INDI, INDO, IHD, INDFUN, L
      DOUBLE PRECISION LOW, RAD, DL, D0F, D1F, D2F, DERV( 2 )
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Don't bother checking ILNR, ILNC, IRNR, IRNC, NDRVS or NBN
C as these errors will all be trapped by NGSVDR
C
C     .
C     . However do check ICMP
C     .
C
      IF ( ICMP.NE.1 .AND. ICMP.NE.2 .AND. ICMP.NE.3 .AND.
     1        ICMP.NE.4 .AND. ICMP.NE.5           )    THEN
         PRINT *,' Subroutine NGSVLP.'
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
         PRINT *,' Subroutine NGSVLP.'
         PRINT *,' IFORMI = ', IFORMI
         PRINT *,' IFORMO = ', IFORMO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( NRI.NE.NR .OR. NRO.NE.NR ) THEN
         PRINT *,' Subroutine NGSVLP.'
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
      IHD = 2
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
             INDI = INDFUN( IR, IHI, INARRI )
             INDO = INDFUN( IR, IHO, INARRO )
             D0F = VI( INDI )
C            .
             L = MHLI( IHI )
             RAD = XARR( IR )
             CALL NGSVDR( VI, IR, IHI, NBN, IHD, NFDCM, NR,
     1                    NDRVS, DERV, ILNC, IRNC, INARRI, FDCM )
             D1F = DERV( 1 )
             D2F = DERV( 2 )
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
