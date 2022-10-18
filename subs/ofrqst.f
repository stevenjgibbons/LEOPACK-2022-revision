C*********************************************************************
C subroutine Outer core Field 2 Radial QST array *********************
C            -          -       -      ---       *********************
C Steve Gibbons Mon Apr  3 08:00:49 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C If the vector V0 contains radial functions to NH spherical harm.s  C
C with type defined by MHT, MHL and MHM, and V1 the first derivs.,   C
C then OFRQST fills the array RQST with the scaloidal, spheroidal    C
C and toroidal coeff.s of the equivalent representation of the vec.  C
C The monopole term is not included as it cannot exist in the        C
C poloidal/toroidal decomposition.                                   C
C                                                                    C
C V0 and V1 are pre-calculated by CASVDR.                            C
C                                                                    C
C Important note!!!                                                  C
C                                                                    C
C  NR is the number of grid nodes in the velocity and temperature    C
C expansions. NRMF is the number of grid nodes in the magnetic       C
C field expansion. Now if NRIC nodes are used to represent the       C
C field in the inner core, then                                      C
C                                                                    C
C  NR + NRIC = NRMF                                                  C
C                                                                    C
C  XARRM( NRIC + 1 ) = RI (radius of inner core)                     C
C  XARRM( NRIC + i ) = XARR( i ) for all 1 .le. i .le. NR.           C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NR        : Number of radial grid nodes for outer core.        C
C     NRMF      : Number of radial grid nodes for mag. field.        C
C     LH        : Maximum spherical harmonic degree, l.              C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C     (note that ILNR and IRNR refer to the outer core - so ILNR = 2 C
C     means we start at magnetic field grid node NRIC + 2 )          C
C                                                                    C
C     INARRM    : Int. parameter array corresponding to V.           C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRMF   See INDFUN for details         C
C                 INARR( 3 ) = NH                                    C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RQST      : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector in outer core.            C
C                  Has dimensions (  LH*(LH+2) ,3, NR ).             C
C              RQST (l*l+2m,1,I) = q_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,1,I) = q_l^mc(r_i)                       C
C              RQST (l*l+2m,2,I) = s_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,2,I) = s_l^mc(r_i)                       C
C              RQST (l*l+2m,3,I) = t_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,3,I) = t_l^mc(r_i)                       C
C                                                                    C
C     V0        : Solution vector. Legnth must be atleast NRMF*NH    C
C     V1        : Deriv.s of s. v. Legnth must be atleast NRMF*NH    C
C                                                                    C
C     XARRM     : Dim ( NRMF ). i^{th} element gives x_i.            C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE OFRQST( NR, NRMF, LH, V0, V1, ILNR, IRNR, INARRM,
     1                   RQST, XARRM, MHT, MHL, MHM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NRMF, LH, ILNR, IRNR, INARRM( * ),
     1        MHT( * ), MHL( * ), MHM( * )
      DOUBLE PRECISION V0( * ), V1( * ), XARRM( NRMF ),
     1                 RQST( LH*(LH+2), 3, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NRR, NH, IFORMF, ICOMP, IR, IHARM, INDFUN, IND,
     1        IPOL, ITOR, IH, L, M, ICS, INDSHC, NRIC, IRMF
      DOUBLE PRECISION ZERO, D0F, D1F, RAD, SQRLL1, TOL
      PARAMETER ( ZERO = 0.0d0, TOL = 1.0d-9 )
C
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameters ...
C
      IFORMF = INARRM( 1 )
      NRR    = INARRM( 2 )
      NH     = INARRM( 3 )
C
      IF ( NRMF.NE.NRR .OR. NRMF.LT.NR ) THEN
        PRINT *,' Subroutine OFRQST'
        PRINT *,' NR          = ', NR
        PRINT *,' NRMF        = ', NRMF
        PRINT *,' INARRM( 2 ) = ', NRR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 ) THEN
        PRINT *,' Subroutine OFRQST.'
        PRINT *,' IFORMF = ', IFORMF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      NRIC = NRMF - NR
C     .
      IPOL = 4
      ITOR = 5
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
        L  = MHL( IH )
        IF ( L.GT.LH ) THEN
          PRINT *,' Subroutine OFRQST.'
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
           DO IR = ILNR, IRNR
C            .
             IRMF = IR + NRIC
             RAD  = XARRM( IRMF )
             IF ( RAD.LT.TOL ) THEN
               PRINT *,' Subroutine OFRQST.'
               PRINT *,' Radius = ', RAD,' at node ',IR
               PRINT *,' Division by zero imminent!'
               PRINT *,' Program aborted.'
               STOP
             ENDIF
C            .
C            . Take derivatives
C            .
             IND = INDFUN( IRMF, IH, INARRM )
             D0F = V0( IND )
             D1F = V1( IND )
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
C            .
             IRMF = IR + NRIC
             IND  = INDFUN( IRMF, IH, INARRM )
             D0F  = V0( IND )
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

