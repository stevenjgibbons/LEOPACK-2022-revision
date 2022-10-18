C*********************************************************************
C subroutine Scalar Function 2 Vector Gradient xtra special vector C *
C            -      -        - -      -                            - *
C Steve Gibbons Tue Feb 13 09:28:18 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Prepares the 2 arrays                                              C
C                                                                    C
C        VGFA( 3, NH3, NTHP )                                        C
C        IVGFA( 2, NH3 )                                             C
C                                                                    C
C  for use by the gradient transform routine SF2VGD.                 C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NH3       : Number of scalar functions                         C
C     ML3       : Array dim ( NH3 ). Sph. harm degree, L.            C
C     MM3       : Array dim ( NH3 ). Sph. harm order, M, or -M.      C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     NTHP      : Number of theta points.                            C
C                                                                    C
C     IVGFA     : Dim ( 2, NH3 ). Locations for SF2VGD.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     VGFA      : Dim ( 3, NH3, NTHP ). Coeff.s for SF2VGD           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SF2VGC( LH, NH3, ML3, MM3, M0, NTHP,
     1                   GAUX, PA, DPA, VGFA, IVGFA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          LH, NH3, ML3( NH3 ), MM3( NH3 ), M0, NTHP,
     1                 IVGFA( 2, NH3 )
      DOUBLE PRECISION GAUX( NTHP ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2, NTHP ),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2, NTHP )
      DOUBLE PRECISION VGFA( 3, NH3, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ITHETA, IP, ICS, L, M, MPS, IH3,
     1                 INDCOS, INDSIN
      DOUBLE PRECISION X, SINE, DALF, DALFD, TERM1
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C ............. start to loop around theta points ...................
      DO ITHETA = 1, NTHP
        X      = GAUX( ITHETA )
        SINE   = DSQRT( 1.0d0 - X*X )
C
C ............. start to loop around Harmonics ......................
C
        DO IH3 = 1, NH3 
C
           L      = ML3( IH3 )
C
C First do monopole term
C (We do not fill in any coefficients for this term)
C
           IF ( L.EQ.0 ) GOTO 50
C
           M      = MM3( IH3 )
           ICS    = 1
C          .
C          . Modify M and ICS if the harmonic is sin (m phi) dep.
C          .
           IF ( M.LT.0   ) ICS = 2
           IF ( ICS.EQ.2 ) M   = -M
C          .
C          . Store associated Legendre Functions
C          .
           IP     = L*(L+1)/2+M+1
           DALF   = PA( IP, ITHETA )
           DALFD  = DPA( IP, ITHETA )
C          .
C          . Make sure that harmonic is
C          . compatible with the symmetry imposition
C          . MPS (pseudo M) = M/M0
C          .
           MPS = M/M0
           IF ( MPS*M0.NE.M ) THEN
             PRINT *,' Subroutine SF2VGC.'
             PRINT *,' MM3(',IH3,') = ', MM3( IH3 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           IF ( ICS.EQ.1 ) THEN
C            .
C            . We are dealing with a cos ( m phi ) term
C            .
             IF ( M.EQ.0 ) THEN
               INDCOS         = 1
               INDSIN         = INDCOS + 1
               VGFA( 1, IH3, ITHETA ) = DALF
               VGFA( 2, IH3, ITHETA ) = DALFD
               VGFA( 3, IH3, ITHETA ) = 0.0d0
               IVGFA( 1, IH3 ) = INDCOS
               IVGFA( 2, IH3 ) = INDSIN
             ELSE
               INDCOS         = 1 + 2*M/M0
               INDSIN         = INDCOS + 1
               TERM1          = (-1.0d0)*DBLE( M )*DALF/SINE
               VGFA( 1, IH3, ITHETA ) = DALF
               VGFA( 2, IH3, ITHETA ) = DALFD
               VGFA( 3, IH3, ITHETA ) = TERM1
               IVGFA( 1, IH3 ) = INDCOS
               IVGFA( 2, IH3 ) = INDSIN
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS         = 1 + 2*M/M0
             INDSIN         = INDCOS + 1
             TERM1          = DBLE( M )*DALF/SINE
             VGFA( 1, IH3, ITHETA ) = DALF
             VGFA( 2, IH3, ITHETA ) = DALFD
             VGFA( 3, IH3, ITHETA ) = TERM1
             IVGFA( 1, IH3 ) = INDSIN
             IVGFA( 2, IH3 ) = INDCOS
C            .
           ENDIF
C          .
 50     CONTINUE
        ENDDO
C ............. ended looping around Harmonics ......................
      ENDDO
C ............. ended looping around theta points ...................
C
      RETURN
      END
C*********************************************************************
