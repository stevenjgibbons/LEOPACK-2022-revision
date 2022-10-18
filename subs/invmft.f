C*********************************************************************
C subroutine Identical Node Velocity cross Magnetic Field Terms ******
C            -         -    -              -        -     -     ******
C Steve Gibbons Thu Jan 27 11:00:31 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Provides, via INNLCA, the matrix terms for the                     C
C curl ( v x B )                                                     C
C terms in the momentum equation.                                    C
C                                                                    C
C All variables other than those listed are dealt with by            C
C INNLCA.                                                            C
C                                                                    C
C In each of the following,                                          C
C                                                                    C
C curl ( V_{A} x B_{B} ) = sum_G B_{G}                               C
C                                                                    C
C  q_A x q_B = 0                                                     C
C  q_A x s_B = sum_G ( C_{qss}^{abg} s_G + C_{qst}^{abg} t_G )       C
C  q_A x t_B = sum_G ( C_{qts}^{abg} s_G + C_{qtt}^{abg} t_G )       C
C  s_A x s_B = sum_G ( C_{ssq}^{abg} q_G )                           C
C  s_A x t_B = sum_G ( C_{stq}^{abg} q_G )                           C
C  t_A x t_B = sum_G ( C_{ttq}^{abg} q_G )                           C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IPARS      : Array containing details of source harmonics.     C
C                                                                    C
C                  IPARS( 1 ) = IQSTA                                C
C                  IPARS( 2 ) = IQSTB                                C
C                  IPARS( 2 ) = IQSTG                                C
C                                                                    C
C This selects the interaction C_{xyz}^{abg} (c.f. equations         C
C (B.45), to (B.50) - and (B.73) to (B.83) in my thesis.             C
C                                                                    C
C                  IPARS( 4 ) = IOLDF                                C
C                                                                    C
C IOLDF = 1 means we are adding the derivative of a velocity A       C
C         function to the matrix terms for a mag. field B term.      C
C                                                                    C
C (IOLDF = 1 is only current option. IOLDF = 2 may be added at       C
C  a later stage if deemed necessary).                               C
C                                                                    C
C                  IPARS( 5 ) = LA (Spherical harm. degree of A )    C
C                  IPARS( 6 ) = LB (Spherical harm. degree of B )    C
C                  IPARS( 7 ) = LG (Spherical harm. degree of G )    C
C                                                                    C
C     IHD        : Highest derivative required of 'new' harmonic     C
C     IHD0       : Highest derivative required of 'old' harmonic     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DPARS      : DPARS( 1 ) contains coefficient C_{xyz}^{abg}     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE INVMFT( CVEC, RAD, IPARS, DPARS, IHD, IRAD, IH0,
     1                   IHD0, NR, NDCS0, IS0, NBN0, INARR0, NDRVS0,
     2                   NDRVM0, NFDCM0, SVFDC0, VEC0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IPARS( * ), IHD, IRAD, IH0, IHD0, NR, NDCS0, IS0, NBN0,
     1        INARR0( * ), NDRVS0, NDRVM0, NFDCM0
      DOUBLE PRECISION CVEC( * ), RAD, DPARS( * ), VEC0( * ),
     1                 SVFDC0( NFDCM0, NR, NDRVM0+1, NDCS0 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IQSTA, IQSTB, IQSTG, IOLDF, LA, LB, LG, I
      DOUBLE PRECISION FAC, DERV( 4 ), RAD2, RAD3, FUNL,
     1                 LOW, FLA, FLB, FLG, SQRLL1, TEMP, TEMP2
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IQSTA  = IPARS( 1 )
      IQSTB  = IPARS( 2 )
      IQSTG  = IPARS( 3 )
      IOLDF  = IPARS( 4 )
      LA     = IPARS( 5 )
      LB     = IPARS( 6 )
      LG     = IPARS( 7 )
C
      FLA    = SQRLL1( LA )
      FLB    = SQRLL1( LB )
      FLG    = SQRLL1( LG )
C
      FAC    = DPARS( 1 )
C
      DO I = 1, IHD+1
        CVEC( I ) = 0.0d0
      ENDDO
C
C Check for zero value of RAD
C
      IF ( RAD.LT.LOW ) THEN
        PRINT *,' Subroutine INVMFT.'
        PRINT *,' RAD = ', RAD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RAD2 = RAD*RAD
      RAD3 = RAD2*RAD
C
C Take derivative of harmonic IH0 at IRAD
C
      CALL ASVDR( VEC0, IRAD, IS0, IH0, NBN0, IHD0, NFDCM0, NR,
     1            NDRVS0, NDRVM0, DERV, INARR0, SVFDC0, NDCS0 )
C
C DERV( 1 ) now contains zero^{th} derivative
C DERV( 2 ) now contains 1^{st} derivative (if requested)
C DERV( 3 ) now contains 2^{rd} derivative (if requested)
C
C Do C_{qst}^{abg} term for old pol. A and new pol. B
C G is poloidal.
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.0 ) THEN
        FUNL  = (-1.0d0)*FLA*FLA*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*DERV( 1 )/RAD2
        CVEC( 2 ) = FUNL*DERV( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{sqt}^{abg} term for old pol. A and new pol. B
C G is poloidal.
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.0 .AND. IHD0.EQ.1 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FLB*FAC/FLG
        TEMP  = DERV( 1 )/RAD + DERV( 2 )
        CVEC( 1 ) = FUNL*TEMP/RAD
        RETURN
      ENDIF
C
C Do C_{ssq}^{abg} term for old pol. A and new pol. B
C G is toroidal.
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.1 ) THEN
        FUNL  = FLA*FLB*FAC
        TEMP  = DERV( 1 )/RAD3 + DERV( 2 )/RAD2
        TEMP2 = DERV( 1 )/RAD2 + DERV( 2 )/RAD
        CVEC( 1 ) = FUNL*TEMP
        CVEC( 2 ) = FUNL*TEMP2
        RETURN
      ENDIF
C
C Do C_{qss}^{abg} term for old pol. A and new pol. B
C G is toroidal.
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.2 .AND. IHD0.EQ.1 ) THEN
        FUNL  = (-1.0d0)*FLA*FLA*FLB*FAC/FLG
        TEMP  = DERV( 2 )/RAD2 - DERV( 1 )/RAD3
        TEMP2 = DERV( 2 )/RAD + DERV( 1 )/RAD2
        CVEC( 1 ) = FUNL*TEMP
        CVEC( 2 ) = FUNL*TEMP2
        CVEC( 3 ) = FUNL*DERV( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{sqs}^{abg} term for old pol. A and new pol. B
C G is toroidal.
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.2 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FLB*FAC/FLG
        TEMP  = DERV( 3 )/RAD + DERV( 2 )/RAD2 - DERV( 1 )/RAD3
        TEMP2 = DERV( 2 )/RAD + DERV( 1 )/RAD2
        CVEC( 1 ) = FUNL*TEMP
        CVEC( 2 ) = FUNL*TEMP2
        RETURN
      ENDIF
C
C Do C_{qtt}^{abg} term for old pol. A and new tor. B
C G is poloidal.
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.0 .AND. IHD0.EQ.0 ) THEN
        FUNL  = FLA*FLA*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*DERV( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{qts}^{abg} term for old pol. A and new tor. B
C G is toroidal.
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.1 ) THEN
        FUNL  = FLA*FLA*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*DERV( 2 )/RAD
        CVEC( 2 ) = FUNL*DERV( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{stq}^{abg} term for old pol. A and new tor. B
C G is toroidal.
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.0 .AND. IHD0.EQ.1 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FAC
        TEMP  = DERV( 1 )/RAD + DERV( 2 )
        CVEC( 1 ) = FUNL*TEMP/RAD
        RETURN
      ENDIF
C
C Do C_{tqt}^{abg} term for old tor. A and new pol. B
C G is poloidal.
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.0 .AND. IHD0.EQ.0 ) THEN
        FUNL  = FLA*FLB*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*DERV( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{tqs}^{abg} term for old tor. A and new pol. B
C G is toroidal.
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.1 ) THEN
        FUNL  = FLA*FLB*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*DERV( 2 )/RAD
        CVEC( 2 ) = FUNL*DERV( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{tsq}^{abg} term for old tor. A and new pol. B
C G is toroidal.
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.0 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FAC
        CVEC( 1 ) = FUNL*DERV( 1 )/RAD2
        CVEC( 2 ) = FUNL*DERV( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{ttq}^{abg} term for old tor. A and new tor. B
C G is toroidal.
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.0 .AND. IHD0.EQ.0 ) THEN
        FUNL  = FLA*FLB*FAC
        CVEC( 1 ) = FUNL*DERV( 1 )/RAD
        RETURN
      ENDIF
C
      PRINT *,' Subroutine INVMFT.'
      PRINT *,' Your parameter choice is illegal.'
      PRINT *,' IQSTA = ', IQSTA,' LA = ', LA
      PRINT *,' IQSTB = ', IQSTB,' LB = ', LB
      PRINT *,' IQSTG = ', IQSTG,' LG = ', LG
      PRINT *,' IOLDF = ', IOLDF,' IHD = ', IHD
      PRINT *,' IHD0  = ', IHD0
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************

