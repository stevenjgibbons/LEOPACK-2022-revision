C*********************************************************************
C subroutine Radial QST array 2 Solution vector: Optimised ***********
C            -      ---         -                -         ***********
C Steve Gibbons Fri Jul 14 09:39:05 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Adds the poloidal and toroidal components of the RQSTO array       C
C to the appropriate elements of the solution vector, V.             C
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
C     INARR     : Int. parameter array corresponding to V.           C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NR     See INDFUN for details         C
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
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHVMFF    : Velocity/magnetic field select.                    C
C                 Set to either 'VEL' or 'MAG'                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RQSTO     : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector.                          C
C                                                                    C
C              RQSTO (I,l*l+2m,1) = q_l^ms(r_i)                      C
C   RQSTO (I,l*l+2m-1+delta_m0,1) = q_l^mc(r_i)                      C
C              RQSTO (I,l*l+2m,2) = s_l^ms(r_i)                      C
C   RQSTO (I,l*l+2m-1+delta_m0,2) = s_l^mc(r_i)                      C
C              RQSTO (I,l*l+2m,3) = t_l^ms(r_i)                      C
C   RQSTO (I,l*l+2m-1+delta_m0,3) = t_l^mc(r_i)                      C
C                                                                    C
C     V         : Solution vector. Legnth must be atleast NR*NH      C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     FAC       : Constant scalar to multiply before addition to V.  C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RQSTSO( NR, LH, ILNR, IRNR, INARR, MHT, MHL, MHM, 
     1                   CHVMFF, RQSTO, V, XARR, FAC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, LH, ILNR, IRNR, INARR( * ),
     1        MHT( * ), MHL( * ), MHM( * )
      DOUBLE PRECISION V( * ), XARR( NR ),
     2                 RQSTO( NR, LH*(LH+2), 3 ), FAC
      CHARACTER *(3) CHVMFF
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NRR, NH, IFORMF, IND, INDFUN, L, M, ICS, IH, IHARM,
     1        IPOL, ITOR, IR, INDSHC
      DOUBLE PRECISION D0F, RAD, SQRLL1
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
      IF ( NR.NE.NRR ) THEN
        PRINT *,' Subroutine RQSTSO'
        PRINT *,' NR         = ', NR
        PRINT *,' INARR( 2 ) = ', NRR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 ) THEN
        PRINT *,' Subroutine RQSTSO.'
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
        PRINT *,' Subroutine RQSTSO.'
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
C     . Loop around the harmonics
C     .
      DO IH = 1, NH
        IF (          MHT( IH ).NE.IPOL     .AND.
     1                MHT( IH ).NE.ITOR     )    GOTO 50
        L = MHL( IH )
        IF ( L.GT.LH ) THEN
          PRINT *,' Subroutine RQSTSO.'
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
             RAD = XARR( IR )
             IND = INDFUN( IR, IH, INARR )
             D0F = RQSTO( IR, IHARM, 1 )
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
             IND = INDFUN( IR, IH, INARR )
             D0F = RQSTO( IR, IHARM, 3 )
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
