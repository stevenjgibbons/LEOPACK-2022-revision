C*********************************************************************
C subroutine Solution Vector 2 Radial QST array **********************
C            -        -        -      ---       **********************
C Steve Gibbons Mon Sep 27 10:42:18 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C If the vector V contains radial functions to NH spherical harm.s   C
C with type defined by MHT, MHL and MHM, then SVRQST fills the       C
C array RQST with the scaloidal, spheroidal and toroidal coeff.s     C
C of the equivalent representation of the vector.                    C
C The monopole term is not included as it cannot exist in the        C
C poloidal/toroidal decomposition.                                   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NR        : Number of radial grid nodes.                       C
C     LH        : Maximum spherical harmonic degree, l.              C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C     ILNC      : First radial node which may be used for f.d. form. C
C     IRNC      : Last radial node which may be used for f.d. form.  C
C                                                                    C
C  ILNC and IRNC correspond to NLMC and NRMC in the call to FDCMBD.  C
C                                                                    C
C     INARR     : Int. parameter array corresponding to V.           C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NR     See INDFUN for details         C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     NFDCM     : Leading dim of FDCM. See FDCMBD.                   C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     NDRVS     : Highest derivative stored in FDCM.                 C
C                 Must be atleast 1.                                 C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHVMFF    : Velocity/magnetic field select.                    C
C                 Set to either 'VEL' or 'MAG'                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RQST      : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector.                          C
C                  Has dimensions (  LH*(LH+2) ,3, NR ).             C
C              RQST (l*l+2m,1,I) = q_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,1,I) = q_l^mc(r_i)                       C
C              RQST (l*l+2m,2,I) = s_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,2,I) = s_l^mc(r_i)                       C
C              RQST (l*l+2m,3,I) = t_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,3,I) = t_l^mc(r_i)                       C
C                                                                    C
C     V         : Solution vector. Legnth must be atleast NR*NH      C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     FDCM      : Dim ( NFDCM, NR, NDRVS ). Finite diff. coeff.s     C
C                Must be formed in advance by a call to FDCMBD.      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVRQST( NR, LH, V, ILNR, IRNR, ILNC, IRNC, INARR,
     1                   NFDCM, NDRVS, FDCM, RQST, NBN, XARR, 
     2                   MHT, MHL, MHM, CHVMFF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, LH, ILNR, IRNR, ILNC, IRNC, INARR( * ), NFDCM, 
     1        NDRVS, NBN, MHT( * ), MHL( * ), MHM( * )
      DOUBLE PRECISION V( * ), XARR( NR ),
     1                 FDCM( NFDCM, NR, NDRVS ),
     2                 RQST( LH*(LH+2), 3, NR )
      CHARACTER *(3) CHVMFF
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NRR, NH, IFORMF, IND, INDFUN, ICOMP, IR, IHARM,
     1        IPOL, ITOR, IH, L, M, ICS, INDSHC, IHD
      DOUBLE PRECISION ZERO, DERV( 1 ), D0F, D1F, RAD, SQRLL1
      PARAMETER ( ZERO = 0.0d0 )
C
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameters ...
C
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C
C Don't bother checking whether ILNR, IRNR are
C compatible with ILNC, IRNC as SVRQST will call 
C NGSVDR which will trap any errors here.
C
      IF ( NR.NE.NRR ) THEN
        PRINT *,' Subroutine SVRQST'
        PRINT *,' NR         = ', NR
        PRINT *,' INARR( 2 ) = ', NRR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 ) THEN
        PRINT *,' Subroutine SVRQST.'
        PRINT *,' IFORMF = ', IFORMF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Let's check validity of CHVMFF
C     .
      IF ( CHVMFF.NE.'VEL' .AND. CHVMFF.NE.'vel' .AND.
     1     CHVMFF.NE.'Vel' .AND. CHVMFF.NE.'MAG' .AND.
     2     CHVMFF.NE.'mag' .AND. CHVMFF.NE.'Mag'  ) THEN
        PRINT *,' Subroutine SVRQST.'
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
C     . Let's zero the array RQST in all
C     . elements we've requested.
C     .
      DO IR = ILNR, IRNR
        DO ICOMP = 1, 3
          DO IHARM = 1, LH*(LH+2)
            RQST( IHARM, ICOMP, IR ) = ZERO
          ENDDO
        ENDDO
      ENDDO
C     .
C     . Loop around the harmonics
C     .
      DO IH = 1, NH
        IF (          MHT( IH ).NE.IPOL     .AND.
     1                MHT( IH ).NE.ITOR     )    GOTO 50
        L = MHL( IH )
        IF ( L.GT.LH ) THEN
          PRINT *,' Subroutine SVRQST.'
          PRINT *,' Harmonic ',IH,' has L = ',L
          PRINT *,' LH = ', LH
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        IF (          MHM( IH ).GE.0      ) THEN
          M =  MHM( IH )
          ICS = 1
        ELSE
          M = -MHM( IH )
          ICS = 2
        ENDIF
        IHARM = INDSHC( L, M, ICS )
c       .
c       . poloidal harmonics
c       .
        IF ( MHT( IH ).EQ.IPOL ) THEN
C          .
           IHD = 1
           DO IR = ILNR, IRNR
C            .
             RAD = XARR( IR )
C            .
C            . Take derivatives
C            .
             IND = INDFUN( IR, IH, INARR )
             D0F = V( IND )
             CALL NGSVDR( V, IR, IH, NBN, IHD, NFDCM, NR,
     1                    NDRVS, DERV, ILNC, IRNC, INARR, FDCM )
             D1F = DERV( 1 )
C            .
C            . Scaloidal part
C            .
             RQST( IHARM, 1, IR ) = DBLE( L*L + L )*D0F/RAD
C            .
C            . Spheroidal part
C            .
             RQST( IHARM, 2, IR ) = 
     1            SQRLL1( L )*(D0F/RAD + D1F)
C            .
           ENDDO
C          .
        ENDIF
c       .
c       . toroidal harmonics
c       .
        IF ( MHT( IH ).EQ.ITOR ) THEN
C          .
           DO IR = ILNR, IRNR
             IND = INDFUN( IR, IH, INARR )
             D0F = V( IND )
C            .
C            . Toroidal part
C            .
             RQST( IHARM, 3, IR ) = (-1.0d0)*SQRLL1( L )*D0F
C            .
           ENDDO
C          .
        ENDIF
c       .
 50   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

