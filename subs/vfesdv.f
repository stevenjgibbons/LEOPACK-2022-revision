C*********************************************************************
C subroutine Vector Function 2 Equal Spectrally Decomposed Vector ****
C            -      -          -     -          -          -      ****
C Steve Gibbons Fri Nov 17 14:08:22 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Performs largely the same task as the subroutine VF2SDV,           C
C except that whereas VF2SDV can have different harmonic sets        C
C for the poloidal and toroidal parts of the vector, VFESDV          C
C assumes that they are identical - this is more efficient if        C
C both equatorial symmetries are included in the solution.           C
C                                                                    C
C The double precision array DVF contains a vector function of       C
C dimensions ( 3, NPHP, NTHP, NR ). VFESDV transforms this back into C
C a spectral space in the following format.                          C
C                                                                    C
C The double precision arrays QRF, SRF and TRF define sets of        C
C radial functions which constitute a vector with the following      C
C formalism:                                                         C
C                                                                    C
C  {\bm v} =  \sum_{ih = 1, NH}   Q_{ih}(r) {\bm q}_{ih} +           C
C             \sum_{ih = 1, NH}   S_{ih}(r) {\bm s}_{ih} +           C
C             \sum_{ih = 1, NH}   T_{ih}(r) {\bm t}_{ih}             C
C                                                                    C
C Now,                                                               C
C                                                                    C
C   {\bm q}_{ih} = [ Y_{ih} ,  0  ,  0 ],                            C
C                                                                    C
C   {\bm s}_{ih} = FAC.[ 0, \partial Y_{ih}/\partial \theta,         C
C                (\sin \theta)^{-1} \partial Y_{ih}/\partial \phi ]  C
C                                                                    C
C      and                                                           C
C                                                                    C
C   {\bm t}_{ih} = FAC.[   0 ,                                       C
C              - (\sin \theta)^{-1} \partial Y_{ih}/\partial \phi,   C
C                        \partial Y_{ih}/\partial \theta ]           C
C                                                                    C
C  where FAC = 1/sqrt( L ( L + 1 ) )                                 C
C      with L = ML( ih ) and                                         C
C  Y_{ih} = P_L^M( \cos \theta ) cos ( M \phi ) for MM( ih ) = M     C
C     or                                                             C
C  Y_{ih} = P_L^M( \cos \theta ) sin ( M \phi ) for MM( ih ) = -M    C
C    where the P_L^M are Schmidt quasi-normalised associated         C
C    Legendre functions (see routine SCHNLA ).                       C
C                                                                    C
C The radial functions Q_{ih}(r), S_{ih}(r) and T_{ih}(r) are stored C
C respectively in QRF, SRF and TRF. The radius at grid node IR is    C
C given by XARR( IR ). The value of Q_{ih}(ir), S_{ih}(ir) and       C
C T_{ih}(ir) are respectively stored in elements                     C
C             QRF( ( ih - 1 )*NR + IR ),                             C
C             SRF( ( ih - 1 )*NR + IR ) and                          C
C             TRF( ( ih - 1 )*NR + IR ).                             C
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
C     NH        : Number of harmonics [number of radial              C
C                  functions Q(r), S(r) and T(r)]                    C
C     ML        : Array dim ( NH ). Sph. harm degree, L.             C
C     MM        : Array dim ( NH ). Sph. harm order, M, or -M.       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     QRF       : Dim ( NR*NH ) Radial functions Q_{ih}(r)           C
C     SRF       : Dim ( NR*NH ) Radial functions S_{ih}(r)           C
C     TRF       : Dim ( NR*NH ) Radial functions T_{ih}(r)           C
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
      SUBROUTINE VFESDV( NTHP, NPHP, M0, LH, NR, ILNR, IRNR, NH, ML,
     1                   MM, QRF, SRF, TRF, GAUX, GAUW, PA, DPA, DVF,
     2                   FTF1, FTF2, FTF3 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHP, NPHP, M0, LH, NR, ILNR, IRNR, NH, ML( NH ),
     1        MM( NH )
      DOUBLE PRECISION QRF( * ), SRF( * ), TRF( * ),
     1                 DVF( 3, NPHP, NTHP, NR), GAUX( NTHP )
      DOUBLE PRECISION FTF1( 2*NPHP ), GAUW( NTHP ),
     1                 FTF2( 2*NPHP ), FTF3( 2*NPHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ISIGN, ITHETA, IR, MPS, ICS, IH, INDLOC,
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
      ILEN = NH*NR
      CALL DVECZ( QRF, ILEN )
      CALL DVECZ( SRF, ILEN )
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
C 
         DO IH = 1, NH
           QRFVAL = ZERO
           SRFVAL = ZERO
           TRFVAL = ZERO
C          .
C          . Calculate factors including L.
C          . Protect against division by zero.
C          .
           L      = ML( IH )
           IF ( L.LT.1 ) THEN
             PRINT *,' Subroutine VFESDV.'
             PRINT *,' ML(',IH,') = ', L
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
           M      = MM( IH )
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
             PRINT *,' Subroutine VFESDV.'
             PRINT *,' MM(',IH,') = ', MM( IH )
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
               TRFVAL  = W2*DALFD*FTF3( INDCOS )
             ELSE
               INDCOS  = 2*MPS + 1
               INDSIN  = 2*MPS + 2
               TERM    = W2*DBLE( M )*DALF/SINE
               QRFVAL  = W1*DALF*FTF1( INDCOS )
               SRFVAL  = W2*DALFD*FTF2( INDCOS ) - TERM*FTF3( INDSIN )
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
             QRFVAL  = W1*DALF*FTF1( INDSIN )
             SRFVAL  = W2*DALFD*FTF2( INDSIN ) + TERM*FTF3( INDCOS )
             TRFVAL  = W2*DALFD*FTF3( INDSIN ) - TERM*FTF2( INDCOS )
C            .
           ENDIF
C          .
           INDLOC        = IH*NR - NR + IR
           QRF( INDLOC ) = QRF( INDLOC ) + QRFVAL
           SRF( INDLOC ) = SRF( INDLOC ) + SRFVAL
           TRF( INDLOC ) = TRF( INDLOC ) + TRFVAL
         ENDDO
C        .
C        . End loop ih = 1, nh
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
