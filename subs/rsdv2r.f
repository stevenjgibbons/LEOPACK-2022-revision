C*********************************************************************
C subroutine Radial Spectrally Decomposed Vector 2 Real space ********
C            -      -          -          -      - -          ********
C Steve Gibbons Sat Jul 22 16:18:28 BST 2000                         C
C____________________________________________________________________C
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
C   {\bm t}_{iht} = FAC.[   0 ,                                      C
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
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C                                                                    C
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
C     DVF       : Vector Function. An array of dimensions            C
C                ( 3, NPHP, NTHP, NR) which contain the R, THETA     C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. DVF( 2, IPHI, ITHETA, NR ) is the Theta  C
C                  compontent of the vector at (ir, iphi, itheta).   C
C                                                                    C
C                 DVF  is completely zeroed on entry to RSDV2R.      C
C____________________________________________________________________C
C                                                                    C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHP )                         C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RSDV2R( NTHP, NPHP, M0, LH, NR, ILNR, IRNR, NPH, MLP,
     1                   MMP, NTH, MLT, MMT, QRF, SRF, TRF, GAUX, PA,
     2                   DPA, DVF, FTF1, FTF2, FTF3 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHP, NPHP, M0, LH, NR, ILNR, IRNR, NPH, MLP( NPH ),
     1        MMP( NPH ), NTH, MLT( NTH ), MMT( NTH )
      DOUBLE PRECISION QRF( * ), SRF( * ), TRF( * ),
     1                 DVF( 3, NPHP, NTHP, NR), GAUX( NTHP )
      DOUBLE PRECISION FTF1( 2*NPHP ),
     1                 FTF2( 2*NPHP ), FTF3( 2*NPHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER          IOP, ISIGN, ITHREE, ITHETA, LENV, IPHI, IR, 
     1                 IP, ICS, L, M, INDLOC, INDCOS, INDSIN, MPS,
     2                 IHP, IHT
      DOUBLE PRECISION ZERO, X, SINE, TERM1, DALF, DALFD,
     1                 DLFAC, DLSQR, QRFVAL, SRFVAL, SRFV2,
     2                 TRFVAL, TRFV2
C
      PARAMETER ( ZERO = 0.0d0, ITHREE = 3 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C   
      IOP = 0
C ............. set DVF array to zero
      CALL QUADOP( DVF, ZERO, ITHREE, NPHP, NTHP, NR, IOP )
C
C     .................. start loop around radial grid nodes
      DO IR = ILNR, IRNR
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
C ............. set all of FTF1, FTF2, FTF3 to zero
C
         LENV = 2*NPHP
C
         CALL VECOP ( FTF1, ZERO, LENV, IOP )
         CALL VECOP ( FTF2, ZERO, LENV, IOP )
         CALL VECOP ( FTF3, ZERO, LENV, IOP )
C
C ............. start to loop around Harmonics ......................
C ............. First do poloidal harmonics (Q and S radial func.s) .
         DO IHP = 1, NPH
           INDLOC = IHP*NR - NR + IR
           QRFVAL = QRF( INDLOC )
           SRFVAL = SRF( INDLOC )
           IF ( QRFVAL.EQ.ZERO .AND. SRFVAL.EQ.ZERO ) GOTO 70
C          .
C          . Calculate factors including L.
C          . Protect against division by zero.
C          .
           L      = MLP( IHP )
           IF ( L.LT.1 ) THEN
             PRINT *,' Subroutine RSDV2R.'
             PRINT *,' MLP(',IHP,') = ', L
             PRINT *,' Program aborted.'
             STOP
           ENDIF
           DLFAC  = DBLE( L )
           DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
           SRFV2  = SRFVAL/DLSQR
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
             PRINT *,' Subroutine RSDV2R.'
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
               FTF1( INDCOS ) = FTF1( INDCOS ) + DALF*QRFVAL
               FTF2( INDCOS ) = FTF2( INDCOS ) + DALFD*SRFV2
             ELSE
               INDCOS         = 2*MPS + 1
               INDSIN         = 2*MPS + 2
               TERM1          = DBLE( M )*DALF/SINE
               FTF1( INDCOS ) = FTF1( INDCOS ) + DALF*QRFVAL
               FTF2( INDCOS ) = FTF2( INDCOS ) + DALFD*SRFV2
               FTF3( INDSIN ) = FTF3( INDSIN ) - TERM1*SRFV2
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS         = 2*MPS + 1
             INDSIN         = 2*MPS + 2
             TERM1          = DBLE( M )*DALF/SINE
             FTF1( INDSIN ) = FTF1( INDSIN ) + DALF*QRFVAL
             FTF2( INDSIN ) = FTF2( INDSIN ) + DALFD*SRFV2
             FTF3( INDCOS ) = FTF3( INDCOS ) + TERM1*SRFV2
C            .
           ENDIF
C          .
 70        CONTINUE
         ENDDO
C        .
C        . End loop ihp = 1, nhp
C
C............. Now do toroidal harmonics (T radial func.s) .
         DO IHT = 1, NTH
           INDLOC = IHT*NR - NR + IR
           TRFVAL = TRF( INDLOC )
           IF ( TRFVAL.EQ.ZERO ) GOTO 71
C          .
C          . Calculate factors including L.
C          . Protect against division by zero.
C          .
           L      = MLT( IHT )
           IF ( L.LT.1 ) THEN
             PRINT *,' Subroutine RSDV2R.'
             PRINT *,' MLT(',IHT,') = ', L
             PRINT *,' Program aborted.'
             STOP
           ENDIF
           DLFAC  = DBLE( L )
           DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
           TRFV2  = TRFVAL/DLSQR
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
             PRINT *,' Subroutine RSDV2R.'
             PRINT *,' MMT(',IHT,') = ', MMT( IHT )
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
               FTF3( INDCOS ) = FTF3( INDCOS ) + DALFD*TRFV2
             ELSE
               INDCOS         = 2*MPS + 1
               INDSIN         = 2*MPS + 2
               TERM1          = DBLE( M )*DALF/SINE
               FTF2( INDSIN ) = FTF2( INDSIN ) + TERM1*TRFV2
               FTF3( INDCOS ) = FTF3( INDCOS ) + DALFD*TRFV2
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS         = 2*MPS + 1
             INDSIN         = 2*MPS + 2
             TERM1          = DBLE( M )*DALF/SINE
             FTF2( INDCOS ) = FTF2( INDCOS ) - TERM1*TRFV2
             FTF3( INDSIN ) = FTF3( INDSIN ) + DALFD*TRFV2
C            .
           ENDIF
C          .
 71        CONTINUE
         ENDDO
C        .
C        . End loop iht = 1, nht
C        .
C ............. ended looping around Harmonics ......................
C
C ............... now perform Fourier Transforms on FTF1, FTF2, FTF3
C
         ISIGN = -1
         CALL FFTRLV( FTF1, NPHP, ISIGN )
         CALL FFTRLV( FTF2, NPHP, ISIGN )
         CALL FFTRLV( FTF3, NPHP, ISIGN )
C
C ...................................................................
         DO IPHI = 1, NPHP
           INDLOC = 2*IPHI - 1
           DVF( 1, IPHI, ITHETA, IR ) = FTF1( INDLOC )
           DVF( 2, IPHI, ITHETA, IR ) = FTF2( INDLOC )
           DVF( 3, IPHI, ITHETA, IR ) = FTF3( INDLOC )
         ENDDO
C
        ENDDO
C ............. ended looping around theta points ...................
      ENDDO
C     . Ended loop around radial grid nodes
C
      RETURN
      END
C*********************************************************************

