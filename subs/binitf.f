C*********************************************************************
C Initial Magnetic Field Vector Form *********************************
C*********************************************************************
      SUBROUTINE BINITF( BINIT, INARR, RI, RO, MHT, MHL, MHM, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * )
      DOUBLE PRECISION BINIT( * ), RI, RO, XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER L, M, ICS, INTPA( 1 ), NR, NH, IMAX, IH, IR,
     1        ITYPE, ICOEF, IND, INDFUN
      PARAMETER ( IMAX = 10 )
      DOUBLE PRECISION PMFIRF, PVSFRF, DPRPA( 1 ), COEFS( IMAX ),
     1                 FAC, RAD
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR = INARR( 2 )
      NH = INARR( 3 )
C
C Loop around harmonics
C
      DO IH = 1, NH
        ITYPE = MHT( IH )
        L     = MHL( IH )
        IF ( MHM( IH ).LT.0 ) THEN
          M   = MHM( IH )*(-1)
          ICS = 2
        ELSE
          M   = MHM( IH )
          ICS = 1
        ENDIF
        INTPA( 1 ) = L
C       .
        IF ( ITYPE.NE.4 .AND. ITYPE.NE.5 ) THEN
          PRINT *,' Subroutine BINITF '
          PRINT *,' ITYPE = ', ITYPE
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
C       . Now we do different things
C       . depending upon whether the harmonic
C       . is poloidal or toroidal - first poloidal
C       .
        IF ( ITYPE.EQ.4 ) THEN
c          .
c          . Here must enter COEF( i ) as a
c          . function of l, m etc. for poloidal radial func.
c          .
           DO ICOEF = 1, IMAX
             FAC = 1.0d0/(2.0d0*ICOEF*ICOEF + 1.0d0)
             FAC = (-1.0d0)**(ICOEF+1)*FAC/DBLE( 2.0d0*L )
             COEFS( ICOEF ) = FAC
           ENDDO
c          .
        ELSE
c          .
c          . Here must enter COEF( i ) as a
c          . function of l, m etc. for toroidal radial func.
c          .
           DO ICOEF = 1, IMAX
             FAC = 1.0d0/(2.0d0*ICOEF*ICOEF + 1.0d0)
             FAC = FAC/DBLE( 2.0d0*L )
             COEFS( ICOEF ) = FAC
           ENDDO
c          .
        ENDIF
C       .
        DO IR = 1, NR
          RAD = XARR( IR )
          IND = INDFUN( IR, IH, INARR )
          IF ( ITYPE.EQ.4 ) FAC =
     1       PMFIRF( RAD, RI, RO, IMAX, COEFS, INTPA, DPRPA )
          IF ( ITYPE.EQ.5 ) FAC =
     1       PVSFRF( RAD, RI, RO, IMAX, COEFS, INTPA, DPRPA )
          BINIT( IND ) = FAC
        ENDDO
C       .
      ENDDO
C
      RETURN
      END
C*********************************************************************
