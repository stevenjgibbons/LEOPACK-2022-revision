C*********************************************************************
C subroutine Vector Function 2 Spectrally Decomposed Vector **********
C            -      -          -          -          -      **********
C Steve Gibbons Wed Nov 15 14:13:44 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C The double precision array DVF contains a vector function of       C
C dimensions ( 3, NPHP, NTHP, NR ). VF2SDV transforms this back into C
C a spectral space in the following format.                          C
C                                                                    C
C The double precision arrays QRF, SRF and TRF define sets of        C
C radial functions which constitute a vector with the following      C
C formalism:                                                         C
C                                                                    C
C  {\bm v} =  \sum_{ihp = 1, NPH}   Q_{ihp}(r) {\bm q}_{ihp} +       C
C             \sum_{ihp = 1, NPH}   S_{ihp}(r) {\bm s}_{ihp} +       C
C             \sum_{iht = 1, NTH}   T_{iht}(r) {\bm t}_{iht}         C
C                                                                    C
C Now,                                                               C
C                                                                    C
C   {\bm q}_{ihp} = [ Y_{ihp} ,  0  ,  0 ],                          C
C                                                                    C
C   {\bm s}_{ihp} = FAC.[ 0, \partial Y_{ihp}/\partial \theta,       C
C                (\sin \theta)^{-1} \partial Y_{ihp}/\partial \phi ] C
C                                                                    C
C  where FAC = 1/sqrt( L ( L + 1 ) )                                 C
C      with L = MLP( ihp ) and                                       C
C  Y_{ihp} = P_L^M( \cos \theta ) cos ( M \phi ) for MMP( ihp ) = M  C
C     or                                                             C
C  Y_{ihp} = P_L^M( \cos \theta ) sin ( M \phi ) for MMP( ihp ) = -M C
C    where the P_L^M are Schmidt quasi-normalised associated         C
C    Legendre functions (see routine SCHNLA ).                       C
C                                                                    C
C      and                                                           C
C                                                                    C
C   {\bm t}_{ihp} = FAC.[   0 ,                                      C
C              - (\sin \theta)^{-1} \partial Y_{iht}/\partial \phi,  C
C                        \partial Y_{iht}/\partial \theta ]          C
C                                                                    C
C  where FAC = 1/sqrt( L ( L + 1 ) )                                 C
C      with L = MLT( iht ) and                                       C
C  Y_{iht} = P_L^M( \cos \theta ) cos ( M \phi ) for MMT( iht ) = M  C
C     or                                                             C
C  Y_{iht} = P_L^M( \cos \theta ) sin ( M \phi ) for MMT( iht ) = -M C
C    where the P_L^M are Schmidt quasi-normalised associated         C
C    Legendre functions (see routine SCHNLA ).                       C
C                                                                    C
C The radial functions Q_{ihp}(r) and S_{ihp}(r) are stored          C
C respectively in QRF and SRF. The radius at grid node IR is given   C
C by XARR( IR ). The value of Q_{ihp}(ir) is stored in the element   C
C             QRF( ( ihp - 1 )*NR + IR )                             C
C and similarly, S_{ihp}(ir) is stored in the element                C
C             SRF( ( ihp - 1 )*NR + IR )                             C
C                                                                    C
C Likewise, the function T_{iht}(r) is stored in TRF with            C
C  T_{iht}(ir) in element  TRF( ( iht - 1 )*NR + IR )                C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NTHP      : Number of theta points.                            C
C     NPHP      : Number of phi points for perform Fourier transform C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NR        : Number of radial grid nodes.                       C
C     ILNR      : Lowest radial node to be acted upon.               C
C     IRNR      : Highest radial node to be acted upon.              C
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
C     QRF       : Dim ( NR*NPH ) Radial functions Q_{iph}(r)         C
C     SRF       : Dim ( NR*NPH ) Radial functions S_{iph}(r)         C
C     TRF       : Dim ( NR*NTH ) Radial functions T_{ith}(r)         C
C                                                                    C
C  (qrf, srf and trf are all set to zero on entry)                   C
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
C     DVF       : Vector Function. An array of dimensions            C
C                ( 3, NPHP, NTHP, NR) which contain the R, THETA     C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. DVF( 2, IPHI, ITHETA, NR ) is the Theta  C
C                  compontent of the vector at (ir, iphi, itheta).   C
C                                                                    C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHP )                         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VF2SDV( NTHP, NPHP, M0, LH, NR, ILNR, IRNR, NPH, MLP,
     1                   MMP, NTH, MLT, MMT, QRF, SRF, TRF, GAUX,
     2                   GAUW, PA, DPA, DVF, FTF1, FTF2, FTF3 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHP, NPHP, M0, LH, NR, ILNR, IRNR, NPH, MLP( NPH ),
     1        MMP( NPH ), NTH, MLT( NTH ), MMT( NTH )
      DOUBLE PRECISION QRF( * ), SRF( * ), TRF( * ),
     1                 DVF( 3, NPHP, NTHP, NR), GAUX( NTHP )
      DOUBLE PRECISION FTF1( 2*NPHP ), GAUW( NTHP ),
     1                 FTF2( 2*NPHP ), FTF3( 2*NPHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ISIGN, ITHETA, IR, MPS, ICS, IHP, IHT, INDLOC,
     1        INDCOS, INDSIN, IP, ILEN, L, M, IPHI
      DOUBLE PRECISION ZERO, X, SINE, TERM, WEIGHT, W1, W2,
     1                 DALF, DALFD, QRFVAL, SRFVAL, TRFVAL,
     2                 DLFAC, DLSQR
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C ............ set arrays QRF, SRF og TRF to zero ....................
      ILEN = NPH*NR
      CALL DVECZ( QRF, ILEN )
      CALL DVECZ( SRF, ILEN )
      ILEN = NTH*NR
      CALL DVECZ( TRF, ILEN )
C     .
C     . Loop around radial grid points
C     .
      DO IR = ILNR, IRNR
C ....................................................................
C ...... now start to loop around theta points .......................
       DO ITHETA = 1, NTHP
C
         X = GAUX( ITHETA )
         SINE   = DSQRT( 1.0d0 - X*X )
         WEIGHT = GAUW( ITHETA )
C        .
C        . Zero the arrays ftf1, ftf2 and ftf3
C        .
         ILEN = 2*NPHP
         CALL DVECZ( FTF1, ILEN )
         CALL DVECZ( FTF2, ILEN )
         CALL DVECZ( FTF3, ILEN )
C
C .................. firstly enter the RVF information into FTF1, FTF2,
C .................. FTF3 for r, theta, phi components respectively of
C .................. VF.
C
         DO IPHI = 1, NPHP
           INDLOC = 2*IPHI - 1
           FTF1( INDLOC ) = DVF( 1, IPHI, ITHETA, IR )
           FTF2( INDLOC ) = DVF( 2, IPHI, ITHETA, IR )
           FTF3( INDLOC ) = DVF( 3, IPHI, ITHETA, IR )
         ENDDO
C
C .................. now perform Forward Discreet Fourier Transforms
C .................. on FTF1, FTF2, FTF3
C
         ISIGN = 1
         CALL FFTRLV ( FTF1, NPHP, ISIGN )
         CALL FFTRLV ( FTF2, NPHP, ISIGN )
         CALL FFTRLV ( FTF3, NPHP, ISIGN )
C ...................................................................
C .                  Now let's loop around the Harmonics..          .
C ...................................................................
C ............. First do poloidal harmonics (Q and S radial func.s) .
         DO IHP = 1, NPH
           QRFVAL = ZERO
           SRFVAL = ZERO
C          .
C          . Calculate factors including L.
C          . Protect against division by zero.
C          .
           L      = MLP( IHP )
           IF ( L.LT.1 ) THEN
             PRINT *,' Subroutine VF2SDV.'
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
             PRINT *,' Subroutine VF2SDV.'
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
               QRFVAL  = W1*DALF*FTF1( INDCOS )
               SRFVAL  = W2*DALFD*FTF2( INDCOS )
             ELSE
               INDCOS  = 2*MPS + 1
               INDSIN  = 2*MPS + 2
               TERM    = W2*DBLE( M )*DALF/SINE
               QRFVAL  = W1*DALF*FTF1( INDCOS )
               SRFVAL  = W2*DALFD*FTF2( INDCOS ) - TERM*FTF3( INDSIN )
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS  = 2*MPS + 1
             INDSIN  = 2*MPS + 2
             TERM    = W2*DBLE( M )*DALF/SINE
             QRFVAL  = W1*DALF*FTF1( INDSIN )
             SRFVAL  = W2*DALFD*FTF2( INDSIN ) + TERM*FTF3( INDCOS )
C            .
           ENDIF
C          .
           INDLOC        = IHP*NR - NR + IR
           QRF( INDLOC ) = QRF( INDLOC ) + QRFVAL
           SRF( INDLOC ) = SRF( INDLOC ) + SRFVAL
         ENDDO
C        .
C        . End loop ihp = 1, nhp
C        .
C............. Now do toroidal harmonics (T radial func.s) .
         DO IHT = 1, NTH
           TRFVAL = ZERO
C          .
C          . Calculate factors including L.
C          . Protect against division by zero.
C          .
           L      = MLT( IHT )
           IF ( L.LT.1 ) THEN
             PRINT *,' Subroutine VF2SDV.'
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
             PRINT *,' Subroutine VF2SDV.'
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
               TRFVAL  = W2*DALFD*FTF3( INDCOS )
             ELSE
               INDCOS  = 2*MPS + 1
               INDSIN  = 2*MPS + 2
               TERM    = W2*DBLE( M )*DALF/SINE
               TRFVAL  = W2*DALFD*FTF3( INDCOS ) + TERM*FTF2( INDSIN )
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS  = 2*MPS + 1
             INDSIN  = 2*MPS + 2
             TERM    = W2*DBLE( M )*DALF/SINE
             TRFVAL  = W2*DALFD*FTF3( INDSIN ) - TERM*FTF2( INDCOS )
C            .
           ENDIF
C          .
           INDLOC        = IHT*NR - NR + IR
           TRF( INDLOC ) = TRF( INDLOC ) + TRFVAL
         ENDDO
C        .
C        . End loop iht = 1, nht
C        .
C ...................................................................
C .                  Ended looping around the harmonics.            .
C ...................................................................
       ENDDO
C ...... ended looping around theta points ...........................
C ....................................................................
      ENDDO
C
C Ended looping around radial grid nodes ...
C
      RETURN
      END
C*********************************************************************
