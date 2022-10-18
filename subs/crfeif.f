C*********************************************************************
C double precision function      Chebyshev Radial Function Energy ****
C                                -         -      -        -      ****
C                                   Integrand Function            ****
C                                   -         -                   ****
C Steve Gibbons Fri May  5 18:20:44 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C The vector CCV is a Chebyshev coefficient vector which is governed C
C by the array INARR. For the radial function IH, CRFEIF will return C
C a value of the integrand for the energy.                           C
C                                                                    C
C If IPARR( 1 ) = 1, we are evaluating energy in a poloidal harm.    C
C If IPARR( 1 ) = 2, we are evaluating energy in a toroidal harm.    C
C                                                                    C
C IPARR( 2 ) = L ... the spherical harmonic degree.                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Format of CCV. See CHINDF.                         C
C     IPARR     : See above.                                         C
C                                                                    C
C     IH        : Number of harmonic to be evaluated.                C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     R         : Radius.                                            C
C                                                                    C
C     CCV       : Chebyshev coefficient vector. Dim ( * ) ( input )  C
C                 Length must be atleast NR*NH                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION CRFEIF( INARR, IPARR, IH, R, CCV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IH, IPARR( * ), INARR( * )
      DOUBLE PRECISION CRFEIF, R, CCV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFLAG, L
      DOUBLE PRECISION Y, DY, D2Y, DKA, SQRLL1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFLAG = IPARR( 1 )
      L     = IPARR( 2 )
C
C Check IFLAG
C
      IF ( IFLAG.NE.1 .AND. IFLAG.NE.2 ) THEN
        PRINT *,' Function CRFEIF.'
        PRINT *,' IFLAG = ', IFLAG
        PRINT *,' Option not (yet) implemented.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Set constants
C
      DKA   = SQRLL1( L )
C
C First, we evaluate the zero^{th}, first and second
C derivatives of the radial function (Y, DY, D2Y respectively)
C
      CALL CCVSPE( INARR, IH, R, CCV, Y, DY, D2Y )
C
C First the case of poloidal vector
C
      IF ( IFLAG.EQ.1 ) THEN
        CRFEIF = DKA*DKA*Y*Y + (Y + R*DY)*(Y + R*DY)
      ENDIF
C
C Now the case of toroidal vector
C
      IF ( IFLAG.EQ.2 ) THEN
        CRFEIF = R*R*Y*Y
      ENDIF
C
      RETURN
      END
C*********************************************************************
