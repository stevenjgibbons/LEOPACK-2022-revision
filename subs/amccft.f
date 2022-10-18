C*********************************************************************
C subroutine Adapted Matrix Curl of Coriolis Force Terms *************
C            -       -      -       -        -     -     *************
C Steve Gibbons Sun Nov 14 17:50:22 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Supplies coefficients to the routine AMLICA for the terms          C
C involving curl ( k x v ) in the vorticity equation.                C
C                                                                    C
C I have derived the forms of the terms in the LATEX user guide      C
C although as an alternative, all these terms are found in eqn.s     C
C (B.55), (B.56), (B.58) and (B.59) of my Ph.D. thesis               C
C                                                                    C
C 'Dynamo Models of the Earth's Magnetic Field'                      C
C     (University of Leeds 1998)                                     C
C                                                                    C
C Set IPARS( 1 ) to  IQSTA  i.e. 1, 2 or 3 depending on whether      C
C the harmonic in the velocity is scaloidal, spheroidal or toroidal. C
C                                                                    C
C Set IPARS( 2 ) to  LA - the spherical harmonic degree of the       C
C input velocity harmonic.                                           C
C                                                                    C
C Set IPARS( 3 ) to  IQSTG  i.e. 1, 2 or 3 depending on whether      C
C the harmonic in the output is scaloidal, spheroidal or toroidal.   C
C                                                                    C
C Set IPARS( 4 ) to  LG - the spherical harmonic degree of the       C
C output function harmonic.                                          C
C                                                                    C
C Set DPARS( 1 ) = K^{\alpha \beta}_{xy} as defined in eqn.s         C
C (B.42), (B.43) and (B.44) of my thesis.                            C
C                                                                    C
C  x is 'q' when IQSTA = 1                                           C
C  x is 's' when IQSTA = 2                                           C
C  x is 't' when IQSTA = 3                                           C
C                                                                    C
C  y is 'q' when IQSTG = 1                                           C
C  y is 's' when IQSTG = 2                                           C
C  y is 't' when IQSTG = 3                                           C
C                                                                    C
C                                                                    C
C                                                                    C
C  All other parameters are looked after by AMLICA.                  C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AMCCFT( CVEC, RAD, IPARS, DPARS, IHD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IPARS( * ), IHD
      DOUBLE PRECISION CVEC( * ), RAD, DPARS( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IHDMIN( 3, 3 ), ND, IQSTA, LA, IQSTG, LG
      DOUBLE PRECISION TOL, FAC, FLA, FLG, SQRLL1
      PARAMETER ( TOL = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      IQSTA = IPARS( 1 )
      LA    = IPARS( 2 )
      IQSTG = IPARS( 3 )
      LG    = IPARS( 4 )
C     .
      FAC   = DPARS( 1 )
C     .
      FLA = SQRLL1( LA )
      FLG = SQRLL1( LG )
C     .
C Check on value of RAD
C     .
      IF ( ABS( RAD ).LT.TOL ) THEN
         PRINT *,' Subroutine AMCCFT.'
         PRINT *,' RAD = ', RAD
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . check value of IHD
C     .
      IHDMIN( 1, 1 ) = 100
C (this routine should not be called with 1 and 1 )
C     .
      IHDMIN( 1, 3 ) = 0
      IHDMIN( 2, 3 ) = 1
      IHDMIN( 2, 1 ) = 1
      IHDMIN( 2, 2 ) = 2
      IHDMIN( 1, 2 ) = 1
      IHDMIN( 3, 3 ) = 0
      IHDMIN( 3, 2 ) = 1
      IHDMIN( 3, 1 ) = 0
C     .
      IF ( IHD.LT.IHDMIN( IQSTA, IQSTG) ) THEN
         PRINT *,' Subroutine AMCCFT.'
         PRINT *,' IHD       = ', IHD
         PRINT *,' IQSTA     = ', IQSTA
         PRINT *,' IQSTG     = ', IQSTG
         PRINT *,' IHDMIN    = ', IHDMIN( IQSTA, IQSTG)
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Zero all coefficients up to IHD
C     .
      DO ND = 1, IHD + 1
        CVEC( ND ) = 0.0d0
      ENDDO
C     .
C     . First K_{qt} term (eq. B.55 in thesis)
C     . 
      IF ( IQSTA.EQ.1 .AND. IQSTG.EQ.3 ) THEN
        CVEC( 1 ) = (-1.0d0)*FLA*FLA*FAC/(RAD*FLG)
        RETURN
      ENDIF
C     .
C     . Now K_{st} term (eq. B.55 in thesis)
C     .
      IF ( IQSTA.EQ.2 .AND. IQSTG.EQ.3 ) THEN
        CVEC( 1 ) = (-1.0d0)*FLA*FAC/(FLG*RAD)
        CVEC( 2 ) = (-1.0d0)*FLA*FAC/FLG
        RETURN
      ENDIF
C     .
C     . Now K_{sq} term (eq. B.56 in thesis)
C     .
      IF ( IQSTA.EQ.2 .AND. IQSTG.EQ.1 ) THEN
        CVEC( 1 ) = FLA*FAC/(RAD*RAD)
        CVEC( 2 ) = FLA*FAC/RAD
        RETURN
      ENDIF
C     .
C     . Now K_{ss} term (eq. B.56 in thesis)
C     .
      IF ( IQSTA.EQ.2 .AND. IQSTG.EQ.2 ) THEN
        CVEC( 2 ) = (-2.0d0)*FLA*FAC/(RAD*FLG)
        CVEC( 3 ) = (-1.0d0)*FLA*FAC/FLG
        RETURN
      ENDIF
C     .
C     . Now K_{qs} term (eq. B.56 in thesis)
C     .
      IF ( IQSTA.EQ.1 .AND. IQSTG.EQ.2 ) THEN
        CVEC( 2 ) = (-1.0d0)*FLA*FLA*FAC/(RAD*FLG)
        RETURN
      ENDIF
C     .
C     . Now K_{tt} term (eq. B.58 in thesis)
C     .
      IF ( IQSTA.EQ.3 .AND. IQSTG.EQ.3 ) THEN
        CVEC( 1 ) = FLA*FAC/FLG
        RETURN
      ENDIF
C     .
C     . Now K_{ts} term (eq. B.59 in thesis)
C     .
      IF ( IQSTA.EQ.3 .AND. IQSTG.EQ.2 ) THEN
        CVEC( 1 ) = FLA*FAC/(FLG*RAD)
        CVEC( 2 ) = FLA*FAC/FLG
        RETURN
      ENDIF
C     .
C     . Now K_{tq} term (eq. B.59 in thesis)
C     .
      IF ( IQSTA.EQ.3 .AND. IQSTG.EQ.1 ) THEN
        CVEC( 1 ) = (-1.0d0)*FLA*FAC/RAD
        RETURN
      ENDIF
C     .
      PRINT *,' Subroutine AMCCFT.'
      PRINT *,' Your options are invalid.'
      PRINT *,' IQSTA = ', IQSTA,' LA = ', LA
      PRINT *,' IQSTG = ', IQSTG,' LG = ', LG
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
