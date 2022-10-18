C*********************************************************************
C subroutine Radial Spectrally Decomposed Vector 2 xtra special A ****
C            -      -          -          -      -              - ****
C Steve Gibbons Mon Feb 12 10:14:54 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Prepares the 5 arrays                                              C
C                                                                    C
C         FTF1P( NPH, NTHP )                                         C
C         FTF2P( NPH, NTHP )                                         C
C         FTF3P( NPH, NTHP )                                         C
C         FTF2T( NTH, NTHP )                                         C
C         FTF3T( NTH, NTHP )                                         C
C                                                                    C
C for use by the vector transform routine RSDV2B.                    C
C                                                                    C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                  Both PA and DPA are formed by the routine SCHNLA. C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     FTF1P     : Dim (NPH,NTHP). Coeff.s required by RSDV2B         C
C     FTF2P     : Dim (NPH,NTHP). Coeff.s required by RSDV2B         C
C     FTF3P     : Dim (NPH,NTHP). Coeff.s required by RSDV2B         C
C     FTF2T     : Dim (NTH,NTHP). Coeff.s required by RSDV2B         C
C     FTF3T     : Dim (NTH,NTHP). Coeff.s required by RSDV2B         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RSDV2A( NTHP, M0, LH, NPH, MLP, MMP, NTH, MLT, MMT,
     1             GAUX, PA, DPA, FTF1P, FTF2P, FTF3P, FTF2T, FTF3T )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NTHP, M0, LH, NPH, MLP( NPH ),
     1                 MMP( NPH ), NTH, MLT( NTH ), MMT( NTH )
      DOUBLE PRECISION FTF1P( NPH, NTHP ), GAUX( NTHP ),
     1                 FTF2P( NPH, NTHP ), FTF3P( NPH, NTHP ),
     2                 FTF2T( NTH, NTHP ), FTF3T( NTH, NTHP ),
     3                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     4                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER          ITHETA, IP, ICS, L, M, MPS, IHP, IHT
      DOUBLE PRECISION ZERO, X, SINE, TERM1, DALF, DALFD,
     1                 DLFAC, DLSQR
C
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C   
C ............. start to loop around theta points ...................
      DO ITHETA = 1, NTHP
C
       X = GAUX( ITHETA )
C        SINE = DSIN( ACOS( X ) )
C replace line above with following line:
C (hopefully! equivalent and faster
C S.J.G. Wed Nov 15 10:29:32 WET 2000
C
       SINE = DSQRT( 1.0d0 - X*X )
C
C ............. start to loop around Harmonics ......................
C ............. First do poloidal harmonics (Q and S radial func.s) .
       DO IHP = 1, NPH
C        .
C        . Calculate factors including L.
C        . Protect against division by zero.
C        .
         L      = MLP( IHP )
         IF ( L.LT.1 ) THEN
           PRINT *,' Subroutine RSDV2A.'
           PRINT *,' MLP(',IHP,') = ', L
           PRINT *,' Program aborted.'
           STOP
         ENDIF
         DLFAC  = DBLE( L )
         DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
C        .
C        . Find wavenumber, m.
C        .
         M      = MMP( IHP )
         ICS    = 1
C        .
C        . Modify M and ICS if the harmonic is sin (m phi) dep.
C        .
         IF ( M.LT.0   ) ICS = 2
         IF ( ICS.EQ.2 ) M   = -M
C        .
C        . Store associated Legendre Functions
C        .
         IP     = L*(L+1)/2+M+1
         DALF   = PA( IP, ITHETA )
         DALFD  = DPA( IP, ITHETA )
C        .
C        . Make sure that harmonic is
C        . compatible with the symmetry imposition
C        . MPS (pseudo M) = M/M0
C        .
         MPS = M/M0
         IF ( MPS*M0.NE.M ) THEN
           PRINT *,' Subroutine RSDV2A.'
           PRINT *,' MMP(',IHP,') = ', MMP( IHP )
           PRINT *,' Program aborted.'
           STOP
         ENDIF
C        .
         IF ( ICS.EQ.1 ) THEN
C          .
C          . We are dealing with a cos ( m phi ) term
C          .
           IF ( M.EQ.0 ) THEN
             FTF1P( IHP, ITHETA ) = DALF
             FTF2P( IHP, ITHETA ) = DALFD/DLSQR
             FTF3P( IHP, ITHETA ) = ZERO
           ELSE
             TERM1          = (-1.0d0)*DBLE( M )*DALF/SINE
             FTF1P( IHP, ITHETA ) = DALF
             FTF2P( IHP, ITHETA ) = DALFD/DLSQR
             FTF3P( IHP, ITHETA ) = TERM1/DLSQR
           ENDIF
C          .
         ELSE
C          .
C          . We are dealing with a sin ( m phi ) term
C          . m should _never_ be zero here!!
C          .
           TERM1          = DBLE( M )*DALF/SINE
           FTF1P( IHP, ITHETA ) = DALF
           FTF2P( IHP, ITHETA ) = DALFD/DLSQR
           FTF3P( IHP, ITHETA ) = TERM1/DLSQR
C          .
         ENDIF
C        .
 70      CONTINUE
       ENDDO
C      .
C      . End loop ihp = 1, nhp
C
C............. Now do toroidal harmonics (T radial func.s) .
       DO IHT = 1, NTH
C        .
C        . Calculate factors including L.
C        . Protect against division by zero.
C        .
         L      = MLT( IHT )
         IF ( L.LT.1 ) THEN
           PRINT *,' Subroutine RSDV2A.'
           PRINT *,' MLT(',IHT,') = ', L
           PRINT *,' Program aborted.'
           STOP
         ENDIF
         DLFAC  = DBLE( L )
         DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
C        .
C        . Find wavenumber, m.
C        .
         M      = MMT( IHT )
         ICS    = 1
C        .
C        . Modify M and ICS if the harmonic is sin (m phi) dep.
C        .
         IF ( M.LT.0   ) ICS = 2
         IF ( ICS.EQ.2 ) M   = -M
C        .
C        . Store associated Legendre Functions
C        .
         IP     = L*(L+1)/2+M+1
         DALF   = PA( IP, ITHETA )
         DALFD  = DPA( IP, ITHETA )
C        .
C        . Make sure that harmonic is
C        . compatible with the symmetry imposition
C        . MPS (pseudo M) = M/M0
C        .
         MPS = M/M0
         IF ( MPS*M0.NE.M ) THEN
           PRINT *,' Subroutine RSDV2A.'
           PRINT *,' MMT(',IHT,') = ', MMT( IHT )
           PRINT *,' Program aborted.'
           STOP
         ENDIF
C        .
         IF ( ICS.EQ.1 ) THEN
C          .
C          . We are dealing with a cos ( m phi ) term
C          .
           IF ( M.EQ.0 ) THEN
             FTF2T( IHT, ITHETA ) = ZERO
             FTF3T( IHT, ITHETA ) = DALFD/DLSQR
           ELSE
             TERM1          = DBLE( M )*DALF/SINE
             FTF2T( IHT, ITHETA ) = TERM1/DLSQR
             FTF3T( IHT, ITHETA ) = DALFD/DLSQR
           ENDIF
C          .
         ELSE
C          .
C          . We are dealing with a sin ( m phi ) term
C          . m should _never_ be zero here!!
C          .
           TERM1          = (-1.0d0)*DBLE( M )*DALF/SINE
           FTF2T( IHT, ITHETA ) = TERM1/DLSQR
           FTF3T( IHT, ITHETA ) = DALFD/DLSQR
C          .
         ENDIF
C        .
 71      CONTINUE
       ENDDO
C        .
C        . End loop iht = 1, nht
C        .
C ............. ended looping around Harmonics ......................
C
      ENDDO
C ............. ended looping around theta points ...................
C
      RETURN
      END
C*********************************************************************

