C*********************************************************************
C subroutine Arbitrary Node Velocity cross Curl of Velocity Terms ****
C            -         -    -              -       -        -     ****
C Steve Gibbons Mon Nov 22 11:17:27 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Provides, via ANNLCA, the matrix terms for the                     C
C -  curl ( v x ( curl v )                                           C
C terms in the momentum equation.                                    C
C                                                                    C
C All variables other than those listed are dealt with by            C
C ANNLCA.                                                            C
C                                                                    C
C In each of the following,                                          C
C                                                                    C
C - curl ( V_{A} x ( curl V_{B} ) = sum_G v_{G}                      C
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
C (B.45), to (B.50) - and (B.62) to (B.72) in my thesis.             C
C                                                                    C
C                  IPARS( 4 ) = IOLDF                                C
C                                                                    C
C IOLDF = 1 means we are adding the derivative of a velocity A       C
C         function to the matrix terms for a curl of vel. B term.    C
C                                                                    C
C IOLDF = 2 means we are adding the derivative of a curl of vel. B   C
C         function to the matrix terms for a velocity A term.        C
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
      SUBROUTINE ANVCVT( CVEC, RAD, IPARS, DPARS, IHD, IH0, IHD0,
     1                   NNDS, XARR0, INARR0, VEC0, WORK1, WORK2,
     2                   IWORK, WORKM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IPARS( * ), IHD, IH0, IHD0, NNDS, INARR0( * ),
     1        IWORK( NNDS )
      DOUBLE PRECISION CVEC( * ), RAD, DPARS( * ), VEC0( * ),
     1                 XARR0( * ), WORK1( NNDS ), WORK2( NNDS ),
     2                 WORKM( NNDS, NNDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IQSTA, IQSTB, IQSTG, IOLDF, LA, LB, LG, I
      DOUBLE PRECISION FAC, RAD2, RAD3, RAD4, FUNL,
     1                 LOW, FLA, FLB, FLG, SQRLL1, TEMP, TEMP2, DL
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
      IF ( NNDS.LT.(IHD0+2) ) THEN
        PRINT *,' Subroutine ANVCVT.'
        PRINT *,' NNDS = ', NNDS
        PRINT *,' IHD0 = ', IHD0
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check for zero value of RAD
C
      IF ( RAD.LT.LOW ) THEN
        PRINT *,' Subroutine ANVCVT.'
        PRINT *,' RAD = ', RAD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RAD2 = RAD*RAD
      RAD3 = RAD*RAD*RAD
      RAD4 = RAD*RAD*RAD*RAD
C
C Take derivative of harmonic IH0 at IRAD
C
      CALL SVRINT( RAD, VEC0, XARR0, INARR0, IH0, NNDS, WORK1,
     1             IWORK, WORK2, WORKM )
C
C WORK1( 1 ) now contains zero^{th} derivative
C WORK1( 2 ) now contains 1^{st} derivative (if requested)
C WORK1( 3 ) now contains 2^{rd} derivative (if requested)
C WORK1( 4 ) now contains 3^{rd} derivative (if requested)
C
C Do C_{qtt}^{abg} term for old pol. A and new pol. B
C G is poloidal. (Equation (B.63) in my thesis)
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.2 .AND. IHD0.EQ.0 ) THEN
        FUNL  = FLA*FLA*FLB*FAC/FLG
        TEMP  = FUNL*WORK1( 1 )
        CVEC( 1 ) = TEMP*(-1.0d0)*FLB*FLB/RAD3
        CVEC( 2 ) = TEMP*2.0d0/RAD2
        CVEC( 3 ) = TEMP/RAD
        RETURN
      ENDIF
C
C Do C_{qtt}^{abg} term for new pol. A and old pol. B
C G is poloidal. (Equation (B.63) in my thesis)
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.0 .AND. IHD0.EQ.2 ) THEN
        FUNL  = FLA*FLA*FLB*FAC/FLG
        TEMP2 = DL( LB, RAD, WORK1( 1 ), WORK1( 2 ), WORK1( 3 ) )
        TEMP  = FUNL*TEMP2
        CVEC( 1 ) = TEMP/RAD
        RETURN
      ENDIF
C
C Do C_{qts}^{abg} term for old pol. A and new pol. B
C G is toroidal.  (Equation (B.64) in my thesis)
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.3 .AND. IHD0.EQ.1 ) THEN
        FUNL  = FLA*FLA*FLB*FAC/FLG
        TEMP  = 2.0d0*WORK1( 1 )/RAD4 - WORK1( 2 )/RAD3
        CVEC( 1 ) = FUNL*TEMP*FLB*FLB
        TEMP  = 2.0d0*WORK1( 2 )/RAD2 - 2.0d0*WORK1( 1 )/RAD3 -
     1                  FLB*FLB*WORK1( 1 )/RAD3
        CVEC( 2 ) = FUNL*TEMP
        TEMP  = 2.0d0*WORK1( 1 )/RAD2 + WORK1( 2 )/RAD
        CVEC( 3 ) = FUNL*TEMP
        CVEC( 4 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{qts}^{abg} term for new pol. A and old pol. B
C G is toroidal.  (Equation (B.64) in my thesis)
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.1 .AND. IHD0.EQ.3 ) THEN
        FUNL  = FLA*FLA*FLB*FAC/FLG
        TEMP  = 2.0d0*FLB*FLB*WORK1( 1 )/RAD4 -
     1          FLB*FLB*WORK1( 2 )/RAD3 - 2.0d0*WORK1( 1 )/RAD3 +
     2          2.0d0*WORK1( 3 )/RAD2 + WORK1( 4 )/RAD
        CVEC( 1 ) = FUNL*TEMP
        TEMP2 = DL( LB, RAD, WORK1( 1 ), WORK1( 2 ), WORK1( 3 ) )
        CVEC( 2 ) = FUNL*TEMP2/RAD
        RETURN
      ENDIF
C
C Do C_{stq}^{abg} term for old pol. A and new pol. B
C G is toroidal.   (Equation (B.64) in my thesis)
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.2 .AND. IHD0.EQ.1 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FAC
        TEMP  = (-1.0d0)*FLB*FLB*( WORK1( 1 )/RAD + WORK1( 2 ) )
        CVEC( 1 ) = FUNL*TEMP/RAD3
        TEMP  = 2.0d0*( WORK1( 1 )/RAD3 + WORK1( 2 )/RAD2 )
        CVEC( 2 ) = FUNL*TEMP
        TEMP  = WORK1( 1 )/RAD2 + WORK1( 2 )/RAD
        CVEC( 3 ) = FUNL*TEMP
        RETURN
      ENDIF
C
C Do C_{stq}^{abg} term for new pol. A and old pol. B
C G is toroidal.   (Equation (B.64) in my thesis)
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.1 .AND. IHD0.EQ.2 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FAC
        TEMP  = DL( LB, RAD, WORK1( 1 ), WORK1( 2 ), WORK1( 3 ) )
        CVEC( 1 ) = FUNL*TEMP/RAD2
        CVEC( 2 ) = FUNL*TEMP/RAD
        RETURN
      ENDIF
C
C Do C_{qst}^{abg} term for old pol. A and new tor. B
C G is poloidal. (Equation (B.66) in my thesis)
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.0 ) THEN
        FUNL  = FLA*FLA*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*WORK1( 1 )/RAD2
        CVEC( 2 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{qst}^{abg} term for new pol. A and old tor. B
C G is poloidal. (Equation (B.66) in my thesis)
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.0 .AND. IHD0.EQ.1 ) THEN
        FUNL  = FLA*FLA*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*( WORK1( 1 )/RAD2 + WORK1( 2 )/RAD )
        RETURN
      ENDIF
C
C Do C_{sqt}^{abg} term for old pol. A and new tor. B
C G is poloidal. (Equation (B.66) in my thesis)
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.0 .AND. IHD0.EQ.1 ) THEN
        FUNL  = FLB*FLB*FLA*FAC/FLG
        CVEC( 1 ) = FUNL*( WORK1( 1 )/RAD2 + WORK1( 2 )/RAD )
        RETURN
      ENDIF
C
C Do C_{sqt}^{abg} term for new pol. A and old tor. B
C G is poloidal. (Equation (B.66) in my thesis)
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.1 .AND. IHD0.EQ.0 ) THEN
        FUNL  = FLB*FLB*FLA*FAC/FLG
        CVEC( 1 ) = FUNL*WORK1( 1 )/RAD2
        CVEC( 2 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{ssq}^{abg} term for old pol. A and new tor. B
C G is toroidal. (Equation (B.67) in my thesis)
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.1 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FAC
        CVEC( 1 ) = FUNL*( WORK1( 1 )/RAD3 + WORK1( 2 )/RAD2 )
        CVEC( 2 ) = FUNL*( WORK1( 1 )/RAD2 + WORK1( 2 )/RAD )
        RETURN
      ENDIF
C
C Do C_{ssq}^{abg} term for new pol. A and old tor. B
C G is toroidal. (Equation (B.67) in my thesis)
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.1 .AND. IHD0.EQ.1 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FAC
        CVEC( 1 ) = FUNL*( WORK1( 1 )/RAD3 + WORK1( 2 )/RAD2 )
        CVEC( 2 ) = FUNL*( WORK1( 1 )/RAD2 + WORK1( 2 )/RAD )
        RETURN
      ENDIF
C
C Do C_{qss}^{abg} term for old pol. A and new tor. B
C G is toroidal. (Equation (B.67) in my thesis)
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.2 .AND. IHD0.EQ.1 ) THEN
        FUNL  = FLA*FLA*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*( WORK1( 2 )/RAD2 - WORK1( 1 )/RAD3 )
        CVEC( 2 ) = FUNL*( WORK1( 2 )/RAD + WORK1( 1 )/RAD2 )
        CVEC( 3 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{qss}^{abg} term for new pol. A and old tor. B
C G is toroidal. (Equation (B.67) in my thesis)
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.1 .AND. IHD0.EQ.2 ) THEN
        FUNL  = FLA*FLA*FLB*FAC/FLG
        TEMP  = WORK1( 3 )/RAD + WORK1( 2 )/RAD2 - WORK1( 1 )/RAD3
        CVEC( 1 ) = FUNL*TEMP
        CVEC( 2 ) = FUNL*( WORK1( 2 )/RAD + WORK1( 1 )/RAD2 )
        RETURN
      ENDIF
C
C Do C_{sqs}^{abg} term for old pol. A and new tor. B
C G is toroidal. (Equation (B.67) in my thesis)
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.2 ) THEN
        FUNL  = FLB*FLB*FLA*FAC/FLG
        TEMP  = WORK1( 3 )/RAD + WORK1( 2 )/RAD2 - WORK1( 1 )/RAD3
        CVEC( 1 ) = FUNL*TEMP
        CVEC( 2 ) = FUNL*( WORK1( 2 )/RAD + WORK1( 1 )/RAD2 )
        RETURN
      ENDIF
C
C Do C_{sqs}^{abg} term for new pol. A and old tor. B
C G is toroidal. (Equation (B.67) in my thesis)
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.2 .AND. IHD0.EQ.1 ) THEN
        FUNL  = FLB*FLB*FLA*FAC/FLG
        CVEC( 1 ) = FUNL*( WORK1( 2 )/RAD2 - WORK1( 1 )/RAD3 )
        CVEC( 2 ) = FUNL*( WORK1( 2 )/RAD + WORK1( 1 )/RAD2 )
        CVEC( 3 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{ttq}^{abg} term for old tor. A and new pol. B
C G is toroidal. (Equation (B.69) in my thesis)
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.2 .AND. IHD0.EQ.0 ) THEN
        FUNL  = FLA*FLB*FAC
        CVEC( 1 ) = FUNL*(-1.0d0)*FLB*FLB*WORK1( 1 )/RAD3
        CVEC( 2 ) = FUNL*(2.0d0)*WORK1( 1 )/RAD2
        CVEC( 3 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{ttq}^{abg} term for new tor. A and old pol. B
C G is toroidal. (Equation (B.69) in my thesis)
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.3 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.0 .AND. IHD0.EQ.2 ) THEN
        FUNL  = FLA*FLB*FAC
        TEMP  = DL( LB, RAD, WORK1( 1 ), WORK1( 2 ), WORK1( 3 ) )
        CVEC( 1 ) = FUNL*TEMP/RAD
        RETURN
      ENDIF
C
C Do C_{tqt}^{abg} term for old tor. A and new tor. B
C G is poloidal. (Equation (B.71) in my thesis)
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.0 .AND. IHD0.EQ.0 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{tqt}^{abg} term for new tor. A and old tor. B
C G is poloidal. (Equation (B.71) in my thesis)
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.3 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.0 .AND. IHD0.EQ.0 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{tqs}^{abg} term for old tor. A and new tor. B
C G is toroidal.  (Equation (B.72) in my thesis)
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.1 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*WORK1( 2 )/RAD
        CVEC( 2 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{tqs}^{abg} term for new tor. A and old tor. B
C G is toroidal.  (Equation (B.72) in my thesis)
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.1 .AND. IQSTG.EQ.2 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.1 .AND. IHD0.EQ.1 ) THEN
        FUNL  = (-1.0d0)*FLA*FLB*FLB*FAC/FLG
        CVEC( 1 ) = FUNL*WORK1( 2 )/RAD
        CVEC( 2 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{tsq}^{abg} term for old tor. A and new tor. B
C G is toroidal.  (Equation (B.72) in my thesis)
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.1 .AND. IHD.EQ.1 .AND. IHD0.EQ.0 ) THEN
        FUNL  = FLA*FLB*FAC
        CVEC( 1 ) = FUNL*WORK1( 1 )/RAD2
        CVEC( 2 ) = FUNL*WORK1( 1 )/RAD
        RETURN
      ENDIF
C
C Do C_{tsq}^{abg} term for new tor. A and old tor. B
C G is toroidal.  (Equation (B.72) in my thesis)
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.2 .AND. IQSTG.EQ.1 .AND.
     1     IOLDF.EQ.2 .AND. IHD.EQ.0 .AND. IHD0.EQ.1 ) THEN
        FUNL  = FLA*FLB*FAC
        CVEC( 1 ) = FUNL*( WORK1( 1 )/RAD2 + WORK1( 2 )/RAD )
        RETURN
      ENDIF
C
      PRINT *,' Subroutine ANVCVT.'
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

