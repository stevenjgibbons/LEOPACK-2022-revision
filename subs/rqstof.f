C*********************************************************************
C subroutine Radial QST array 2 Outer core magnetic Field ************
C            -      ---         -                   -     ************
C Steve Gibbons Mon Apr  3 08:59:40 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Adds the poloidal and toroidal components of the RQST array        C
C to the appropriate elements of the solution vector.                C
C                                                                    C
C There are NR grid nodes in the outer core and NRMF total grid      C
C nodes for the magnetic field. Then NRIC, the number of nodes for   C
C the inner core, is NRMF - NR.                                      C
C                                                                    C
C Node IR in RQST then corresponds to node IR + NRIC in the field    C
C vector.                                                            C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NR        : Number of radial grid nodes outer core.            C
C     NRMF      : Number of radial nodes for magnetic field vector.  C
C     LH        : Maximum spherical harmonic degree, l.              C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C     INARRM    : Int. parameter array corresponding to V.           C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARRM( 1 ) = IFORMF                               C
C                 INARRM( 2 ) = NRMF   See INDFUN for details        C
C                 INARRM( 3 ) = NH                                   C
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
C                  decomposition of vector.                          C
C                  Has dimensions (  LH*(LH+2) ,3, NR ).             C
C              RQST (l*l+2m,1,I) = q_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,1,I) = q_l^mc(r_i)                       C
C              RQST (l*l+2m,2,I) = s_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,2,I) = s_l^mc(r_i)                       C
C              RQST (l*l+2m,3,I) = t_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,3,I) = t_l^mc(r_i)                       C
C                                                                    C
C     V         : Solution vector. Legnth must be atleast NRMF*NH    C
C                                                                    C
C     XARRM     : Dim ( NRMF ). i^{th} element gives x_i.            C
C     FAC       : Constant scalar to multiply before addition to V.  C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RQSTOF( NR, NRMF, LH, ILNR, IRNR, INARRM, MHT,
     1                   MHL, MHM, RQST, V, XARRM, FAC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NRMF, LH, ILNR, IRNR, INARRM( * ),
     1        MHT( * ), MHL( * ), MHM( * )
      DOUBLE PRECISION V( * ), XARRM( NR ),
     2                 RQST( LH*(LH+2), 3, NR ), FAC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NRR, NH, IFORMF, IND, INDFUN, L, M, ICS, IH, IHARM,
     1        IPOL, ITOR, IR, INDSHC, NRIC, IRMF
      DOUBLE PRECISION D0F, RAD, SQRLL1
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
        PRINT *,' Subroutine RQSTOF'
        PRINT *,' NR         = ', NR
        PRINT *,' NRMF       = ', NRMF
        PRINT *,' INARR( 2 ) = ', NRR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 ) THEN
        PRINT *,' Subroutine RQSTOF.'
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
C     . Loop around the harmonics
C     .
      DO IH = 1, NH
        IF (          MHT( IH ).NE.IPOL     .AND.
     1                MHT( IH ).NE.ITOR     )    GOTO 50
        L = MHL( IH )
        IF ( L.GT.LH ) THEN
          PRINT *,' Subroutine RQSTOF.'
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
             IRMF = IR + NRIC
             RAD  = XARRM( IRMF )
             IND  = INDFUN( IRMF, IH, INARRM )
             D0F  = RQST( IHARM, 1, IR )
             V( IND ) = V( IND ) + FAC*D0F*RAD/DBLE( L*L + L )
           ENDDO
C          .
        ENDIF
c       .
c       . toroidal harmonics
c       .
        IF ( MHT( IH ).EQ.ITOR ) THEN
C          .
           DO IR = ILNR, IRNR
             IRMF = IR + NRIC
             IND  = INDFUN( IRMF, IH, INARRM )
             D0F  = RQST( IHARM, 3, IR )
             V( IND ) = V( IND ) - FAC*D0F/SQRLL1( L )
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
