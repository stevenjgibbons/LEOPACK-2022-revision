C*********************************************************************
C subroutine Identical Node Velocity dot Gradient of Inhomogeneous  **
C            -         -    -            -           -              **
C                                                 Temperature terms **
C                                                 -                 **
C Steve Gibbons Mon Nov 22 08:00:19 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Provides, via INNLCA, the matrix terms for the v.Grad( T )         C
C terms in the heat equation.                                        C
C                                                                    C
C T is the temperature given by theta + f( r ) where                 C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C (see ITFA)                                                         C
C                                                                    C
C All variables other than those listed are dealt with by            C
C INNLCA.                                                            C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IPARS      : Array containing details of source harmonics.     C
C                                                                    C
C                  IPARS( 1 ) = IQSTA                                C
C                  IPARS( 2 ) = IQSTB                                C
C                                                                    C
C 'A' is the velocity harmonic and 'B' is the temperature harmonic.  C
C                                                                    C
C  IQSTA = 1 for x = 'q', = 2 for x = 's' and 3 for x = 't'.         C
C  IQSTB = 1 for y = 'q', and 2 for y = 's'. (IQSTB is never 3).     C
C                                                                    C
C This selects the interaction S_{xy}^{abg} (c.f. equations          C
C (B.39), (B.40) and (B.41) - and (B.52) to (B.53) in my thesis.     C
C                                                                    C
C                  IPARS( 3 ) = IOLDF                                C
C                                                                    C
C IOLDF = 1 means we are adding the derivative of a velocity         C
C         function to the matrix element for a temperature term.     C
C                                                                    C
C IOLDF = 2 means we are adding the derivative of a temperature      C
C         function to the matrix element for a velocity term.        C
C                                                                    C
C                  IPARS( 4 ) = LA (Spherical harm. degree of A )    C
C                  IPARS( 5 ) = LB (Spherical harm. degree of B )    C
C                                                                    C
C     IHD        : Highest derivative required of 'new' harmonic     C
C     IHD0       : Highest derivative required of 'old' harmonic     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DPARS      : DPARS( 1 ) contains coefficient S_{xy}^{abg}      C
C                  DPARS( 2 ) contains RI (inner boundary radius)    C
C                  DPARS( 3 ) contains RO (inner boundary radius)    C
C                  DPARS( 4 ) contains CA (see above)                C
C                  DPARS( 5 ) contains CB (see above)                C
C                  DPARS( 6 ) contains CC (see above)                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE INVGIT( CVEC, RAD, IPARS, DPARS, IHD, IRAD, IH0,
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
      INTEGER IQSTA, IQSTB, IOLDF, LA, LB, I
      DOUBLE PRECISION FAC, DERV( 2 ), RAD2, LOW, FLA, FLB,
     1                 SQRLL1, RI, RO, CA, CB, CC
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IQSTA  = IPARS( 1 )
      IQSTB  = IPARS( 2 )
      IOLDF  = IPARS( 3 )
      LA     = IPARS( 4 )
      LB     = IPARS( 5 )
C
      FLA    = SQRLL1( LA )
      FLB    = SQRLL1( LB )
C
      FAC    = DPARS( 1 )
      RI     = DPARS( 2 )
      RO     = DPARS( 3 )
      CA     = DPARS( 4 )
      CB     = DPARS( 5 )
      CC     = DPARS( 6 )
C
      DO I = 1, IHD+1
        CVEC( I ) = 0.0d0
      ENDDO
C
C Check for zero value of RAD
C
      IF ( RAD.LT.LOW ) THEN
        PRINT *,' Subroutine INVGIT.'
        PRINT *,' RAD = ', RAD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RAD2 = RAD*RAD
C
C Take derivative of harmonic IH0 at IRAD
C
      CALL ASVDR( VEC0, IRAD, IS0, IH0, NBN0, IHD0, NFDCM0, NR,
     1            NDRVS0, NDRVM0, DERV, INARR0, SVFDC0, NDCS0 )
C
C DERV( 1 ) now contains zero^{th} derivative of theta
C DERV( 2 ) now contains 1^{st} derivative (if requested) of theta
C Now add on terms from f(r) if so required
C
      IF ( IOLDF.EQ.2 ) THEN
        CALL ITFA( RAD, RI, RO, CA, CB, CC, DERV, IHD0 )
      ENDIF
C
C Do S_{qq}^{abg} term for new theta_b and old p_a
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.1 .AND. IOLDF.EQ.1 .AND.
     1     IHD.EQ.1 .AND. IHD0.EQ.0 ) THEN
        CVEC( 2 ) = FLA*FLA*FAC*DERV( 1 )/RAD
        RETURN
      ENDIF
C
C Do S_{qq}^{abg} term for old theta_b and new p_a
C
      IF ( IQSTA.EQ.1 .AND. IQSTB.EQ.1 .AND. IOLDF.EQ.2 .AND.
     1     IHD.EQ.0 .AND. IHD0.EQ.1 ) THEN
        CVEC( 1 ) = FLA*FLA*FAC*DERV( 2 )/RAD
        RETURN
      ENDIF
C
C Do S_{ss}^{abg} term for new theta_b and old p_a
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.2 .AND. IOLDF.EQ.1 .AND.
     1     IHD.EQ.0 .AND. IHD0.EQ.1 ) THEN
        CVEC( 1 ) = FLA*FLB*FAC*( DERV(2)/RAD + DERV(1)/RAD2 )
        RETURN
      ENDIF
C
C Do S_{ss}^{abg} term for old theta_b and new p_a
C
      IF ( IQSTA.EQ.2 .AND. IQSTB.EQ.2 .AND. IOLDF.EQ.2 .AND.
     1     IHD.EQ.1 .AND. IHD0.EQ.0 ) THEN
        CVEC( 1 ) = FLA*FLB*FAC*DERV(1)/RAD2
        CVEC( 2 ) = FLA*FLB*FAC*DERV(1)/RAD
        RETURN
      ENDIF
C
C Do S_{ts}^{abg} term for new theta_b and old p_a
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.2 .AND. IOLDF.EQ.1 .AND.
     1     IHD.EQ.0 .AND. IHD0.EQ.0 ) THEN
        CVEC( 1 ) = (-1.0d0)*FLA*FLB*FAC*DERV(1)/RAD
        RETURN
      ENDIF
C
C Do S_{ts}^{abg} term for old theta_b and new p_a
C
      IF ( IQSTA.EQ.3 .AND. IQSTB.EQ.2 .AND. IOLDF.EQ.2 .AND.
     1     IHD.EQ.0 .AND. IHD0.EQ.0 ) THEN
        CVEC( 1 ) = (-1.0d0)*FLA*FLB*FAC*DERV(1)/RAD
        RETURN
      ENDIF
C
      PRINT *,' Subroutine INVGIT.'
      PRINT *,' Your parameter choice is illegal.'
      PRINT *,' IQSTA = ', IQSTA,' LA = ', LA
      PRINT *,' IQSTB = ', IQSTB,' LB = ', LB
      PRINT *,' IOLDF = ', IOLDF,' IHD = ', IHD
      PRINT *,' IHD0  = ', IHD0
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************

