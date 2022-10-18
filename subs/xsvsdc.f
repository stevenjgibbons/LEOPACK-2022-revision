C*********************************************************************
C subroutine Xtra Special Vector function 2 Spectrally Decomp. vec A *
C            -    -       -                 -          -           - *
C Steve Gibbons Tue Feb 13 14:02:50 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Creates the four arrays                                            C
C                                                                    C
C       PFA( 3, NPH, NTHP ), JPFA( 2, NPH ),                         C
C       TFA( 2, NTH, NTHP ), JTFA( 2, NTH )                          C
C                                                                    C
C required by the transform routine XSVSDD.                          C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NTHP      : Number of theta points.                            C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NPH       : Number of poloidal harmonics [number of radial     C
C                  functions Q(r) and S(r)]                          C
C     MLP       : Array dim ( NPH ). Sph. harm degree, L.            C
C     MMP       : Array dim ( NPH ). Sph. harm order, M, or -M.      C
C     NTH       : Number of toroidal harmonics [number of radial     C
C                  functions T(r)]                                   C
C     MLT       : Array dim ( NTH ). Sph. harm degree, L.            C
C     MMT       : Array dim ( NTH ). Sph. harm order, M, or -M.      C
C                                                                    C
C     JPFA      : Dim(2,NPH). Locations for use by XSVSDD            C
C     JTFA      : Dim(2,NTH). Locations for use by XSVSDD            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUX      : Dim (NTHP). Gauss weights from the GAUWTS routine. C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                  Both PA and DPA are formed by the routine SCHNLA. C
C                                                                    C
C     PFA       : Dim ( 3, NPH, NTHP ). Coeffs for XSVSDD.           C
C     TFA       : Dim ( 2, NTH, NTHP ). Coeffs for XSVSDD.           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XSVSDC( NTHP, M0, LH, NPH, MLP, MMP, NTH, MLT, MMT,
     1                   GAUX, GAUW, PA, DPA, PFA, TFA, JPFA, JTFA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHP, M0, LH, NPH, MLP( NPH ), MMP( NPH ), NTH,
     1        MLT( NTH ), MMT( NTH ), JPFA( 2, NPH ), JTFA( 2, NTH )
      DOUBLE PRECISION GAUX( NTHP ), GAUW( NTHP ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP ),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP )
      DOUBLE PRECISION PFA( 3, NPH, NTHP ), TFA( 2, NTH, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHETA, MPS, ICS, IHP, IHT, IP, L, M, INDCOS, INDSIN
      DOUBLE PRECISION X, SINE, TERM, WEIGHT, W1, W2,
     1                 DALF, DALFD, DLFAC, DLSQR
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C ....................................................................
C ...... now start to loop around theta points .......................
      DO ITHETA = 1, NTHP
C
         X = GAUX( ITHETA )
         SINE   = DSQRT( 1.0d0 - X*X )
         WEIGHT = GAUW( ITHETA )
C        .
C ...................................................................
C .                  Now let's loop around the Harmonics..          .
C ...................................................................
C ............. First do poloidal harmonics (Q and S radial func.s) .
         DO IHP = 1, NPH
C          .
C          . Calculate factors including L.
C          . Protect against division by zero.
C          .
           L      = MLP( IHP )
           IF ( L.LT.1 ) THEN
             PRINT *,' Subroutine XSVSDC.'
             PRINT *,' MLP(',IHP,') = ', L
             PRINT *,' Program aborted.'
             STOP
           ENDIF
           DLFAC  = DBLE( L )
           DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
           W1     = 0.25d0*WEIGHT*(DLFAC+DLFAC+1.0d0)
           W2     = W1/DLSQR
C          .
C          . Find wavenumber, m.
C          .
           M      = MMP( IHP )
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
             PRINT *,' Subroutine XSVSDC.'
             PRINT *,' MMP(',IHP,') = ', MMP( IHP )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           IF ( ICS.EQ.1 ) THEN
C            .
C            . We are dealing with a cos ( m phi ) term
C            .
             IF ( M.EQ.0 ) THEN
               INDCOS  = 1
               INDSIN  = INDCOS + 1
               PFA( 1, IHP, ITHETA ) = W1*DALF
               PFA( 2, IHP, ITHETA ) = W2*DALFD
               PFA( 3, IHP, ITHETA ) = 0.0d0
               JPFA( 1, IHP ) = INDCOS
               JPFA( 2, IHP ) = INDSIN
             ELSE
               INDCOS  = 1 + 2*M/M0
               INDSIN  = INDCOS + 1
               TERM    = (-1.0d0)*W2*DBLE( M )*DALF/SINE
               PFA( 1, IHP, ITHETA ) = W1*DALF
               PFA( 2, IHP, ITHETA ) = W2*DALFD
               PFA( 3, IHP, ITHETA ) = TERM
               JPFA( 1, IHP ) = INDCOS
               JPFA( 2, IHP ) = INDSIN
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS  = 1 + 2*M/M0
             INDSIN  = INDCOS + 1
             TERM    = W2*DBLE( M )*DALF/SINE
             PFA( 1, IHP, ITHETA ) = W1*DALF
             PFA( 2, IHP, ITHETA ) = W2*DALFD
             PFA( 3, IHP, ITHETA ) = TERM
             JPFA( 1, IHP ) = INDSIN
             JPFA( 2, IHP ) = INDCOS
C            .
           ENDIF
C          .
         ENDDO
C        .
C        . End loop ihp = 1, nhp
C        .
C............. Now do toroidal harmonics (T radial func.s) .
         DO IHT = 1, NTH
C          .
C          . Calculate factors including L.
C          . Protect against division by zero.
C          .
           L      = MLT( IHT )
           IF ( L.LT.1 ) THEN
             PRINT *,' Subroutine XSVSDC.'
             PRINT *,' MLT(',IHT,') = ', L
             PRINT *,' Program aborted.'
             STOP
           ENDIF
           DLFAC  = DBLE( L )
           DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
           W1     = 0.25d0*WEIGHT*(DLFAC+DLFAC+1.0d0)
           W2     = W1/DLSQR
C          .
C          . Find wavenumber, m.
C          .
           M      = MMT( IHT )
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
             PRINT *,' Subroutine XSVSDC.'
             PRINT *,' MMP(',IHP,') = ', MMP( IHP )
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
               TFA( 1, IHT, ITHETA ) = 0.0d0
               TFA( 2, IHT, ITHETA ) = W2*DALFD
               JTFA( 1, IHT ) = INDSIN
               JTFA( 2, IHT ) = INDCOS
             ELSE
               INDCOS  = 1 + 2*M/M0
               INDSIN  = INDCOS + 1
               TERM    = W2*DBLE( M )*DALF/SINE
               TFA( 1, IHT, ITHETA ) = TERM
               TFA( 2, IHT, ITHETA ) = W2*DALFD
               JTFA( 1, IHT ) = INDSIN
               JTFA( 2, IHT ) = INDCOS
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS  = 1 + 2*M/M0
             INDSIN  = INDCOS + 1
             TERM    = (-1.0d0)*W2*DBLE( M )*DALF/SINE
             TFA( 1, IHT, ITHETA ) = TERM
             TFA( 2, IHT, ITHETA ) = W2*DALFD
             JTFA( 1, IHT ) = INDCOS
             JTFA( 2, IHT ) = INDSIN
C            .
           ENDIF
C          .
         ENDDO
C        .
C        . End loop iht = 1, nht
C        .
C ...................................................................
C .                  Ended looping around the harmonics.            .
C ...................................................................
      ENDDO
C ...... ended looping around theta points ...........................
C
      RETURN
      END
C*********************************************************************
