C*********************************************************************
C subroutine Solution Vector PRiNT ***********************************
C            -        -      -- -- ***********************************
C Steve Gibbons Fri Sep 24 13:43:53 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Outputs an eye-readable breakdown of a solution vector with the    C
C radial function of each component being headed by its              C
C description.                                                       C
C                                                                    C
C Writes out to file with logical unit number, LU.                   C
C                                                                    C
C If LU = 6, it goes to standard output, otherwise it goes to        C
C an opened file. SVPRNT does NOT open files.                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARR     : Int. parameter array corresponding to V.           C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR     See INDFUN for details        C
C                 INARR( 3 ) = NH      nrr must equal nr             C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C     MHT( ih ) = 1 if 'ih' is poloidal velocity harmonic.           C
C     MHT( ih ) = 2 if 'ih' is toroidal velocity harmonic.           C
C     MHT( ih ) = 3 if 'ih' is temperature harmonic.                 C
C     MHT( ih ) = 4 if 'ih' is poloidal magnetic field harmonic.     C
C     MHT( ih ) = 5 if 'ih' is toroidal magnetic field harmonic.     C
C                                                                    C
C     MHL       : Spherical harmonic degree, l, of harmonic 'ih'.    C
C                                                                    C
C     MHM       : Spherical harmonic degree, m, of harmonic 'ih'     C
C                if harmonic has (cos m phi dependence) - otherwise  C
C                MHM( ih ) = -m                                      C
C                                                                    C
C     ILNR      : Lowest radial node to output.                      C
C     IRNR      : Highest radial node to output.                     C
C                                                                    C
C     LU        : Logical unit number of file for output.            C
C                                                                    C
C     IZF       : = 1 to display ALL vectors                         C
C                 = 2 to display only non-zero harmonics             C
C                                                                    C
C     IMODE     : Integer flag which sets the format of output       C
C                 Currently, option is limited to IMODE = 1.         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V         : DP vector of dimension ( * ) containing solution   C
C                 Length must be atleast NH*NR                       C
C                                                                    C
C     DPARR     : Array for use by RADVLF                            C
C                                                                    C
C     XARR      : Array ( * ). Only referenced if IFORMF = 3 or 4.   C
C                 In this case dim( xarr ) must be atleast NR and    C
C                 xarr( IR ) will contain the x value at node IR.    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVPRNT ( V, NR, INARR, DPARR, XARR, MHT, MHL, MHM,
     1                    ILNR, IRNR, LU, IZF, IMODE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MHT( * ), MHL( * ), MHM( * ), ILNR,
     1        IRNR, LU, IZF, IMODE
      DOUBLE PRECISION V( * ), DPARR( * ), XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RAD, TOL, HTOT, H
      INTEGER IND, IH, NH, NRR, IR, M, INDFUN
      PARAMETER ( TOL = 1.0d-6 )
      CHARACTER *(23) CHPV, CHTV, CHTH, CHPM, CHTM, CH23
      CHARACTER *(3) CHCOS, CHSIN, CH3
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IMODE.NE.1 ) THEN
        PRINT *,' Subroutine SVPRNT.'
        PRINT *,' IMODE = ', IMODE
        PRINT *,' IMODE = 1 only current option.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NRR = INARR( 2 )
      NH  = INARR( 3 )
      IF ( NRR.NE.NR ) THEN
        PRINT *,' Subroutine SVPRNT.'
        PRINT *,' NRR = ', NRR,' NR = ', NR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( IZF.EQ.2 ) THEN
        WRITE ( LU, 82 ) '# Non-zero radial functions only,'
      ENDIF
C
      CHCOS = 'Cos'
      CHSIN = 'Sin'
C
      CHPV = 'Poloidal velocity      '
      CHTV = 'Toroidal velocity      '
      CHTH = 'Temperature            '
      CHPM = 'Poloidal magnetic field'
      CHTM = 'Toroidal magnetic field'
C             00000000011111111112222
C             12345678901234567890123
C     .
      DO IH = 1, NH
        IF ( IZF.EQ.2 ) THEN
           HTOT = 0.0d0
           DO IR = 1, NR
             IND = INDFUN( IR, IH, INARR )
             HTOT = HTOT + DABS( V( IND ) )
           ENDDO
           IF ( HTOT.LT.TOL ) GOTO 100
        ENDIF
C       .
C       . Decide nature of our harmonic
C       .
        M = IABS( MHM( IH ) )
C       .
        IF ( MHT( IH ).EQ.1 ) CH23 = CHPV
        IF ( MHT( IH ).EQ.2 ) CH23 = CHTV
        IF ( MHT( IH ).EQ.3 ) CH23 = CHTH
        IF ( MHT( IH ).EQ.4 ) CH23 = CHPM
        IF ( MHT( IH ).EQ.5 ) CH23 = CHTM
C       .
        IF ( MHM( IH ).GE.0 ) CH3  = CHCOS
        IF ( MHM( IH ).LT.0 ) CH3  = CHSIN
C       .
        IF ( IMODE.EQ.1 ) WRITE ( LU, 88 )
     1       CH23, MHL( IH ), M, CH3
        IF ( IMODE.EQ.1 ) WRITE ( LU, 89 )
C       .
        DO IR = ILNR, IRNR
          CALL RADVLF( RAD, IR, INARR, DPARR, XARR, H )
          IND = INDFUN( IR, IH, INARR )
          IF ( IMODE.EQ.1 ) WRITE ( LU, 90 ) RAD, V( IND )
        ENDDO
C       .
 100    CONTINUE
      ENDDO
 82   FORMAT(A)
 88   FORMAT('#  ',A23,' L= ',I3,' M= ',I3,' ',A3)
 89   FORMAT('#  ----------------------------------------')
 90   FORMAT(f20.9,f20.9)
C     .
      RETURN
      END
C*********************************************************************
