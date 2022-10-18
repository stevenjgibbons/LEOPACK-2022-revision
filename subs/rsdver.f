C*********************************************************************
C subroutine Radial Spectrally Decomposed Vector 2 Equal Real space **
C            -      -          -          -        -     -          **
C Steve Gibbons Fri Nov 17 13:16:29 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Performs largely the same task as the subroutine RSDV2R,           C
C except that whereas RSDV2R can have different harmonic sets        C
C for the poloidal and toroidal parts of the vector, RSDVER          C
C assumes that they are identical - this is more efficient if        C
C both equatorial symmetries are included in the solution.           C
C                                                                    C
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
C     DVF       : Vector Function. An array of dimensions            C
C                ( 3, NPHP, NTHP, NR) which contain the R, THETA     C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. DVF( 2, IPHI, ITHETA, NR ) is the Theta  C
C                  compontent of the vector at (ir, iphi, itheta).   C
C                                                                    C
C                 DVF  is completely zeroed on entry to RSDVER.      C
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
      SUBROUTINE RSDVER( NTHP, NPHP, M0, LH, NR, ILNR, IRNR, NH, ML,
     1                   MM, QRF, SRF, TRF, GAUX, PA, DPA, DVF, FTF1,
     2                   FTF2, FTF3 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHP, NPHP, M0, LH, NR, ILNR, IRNR, NH, ML( NH ),
     1        MM( NH )
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
     1                 IP, ICS, L, M, INDLOC, INDCOS, INDSIN, MPS, IH
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
C
         DO IH = 1, NH
           INDLOC = IH*NR - NR + IR
           QRFVAL = QRF( INDLOC )
           SRFVAL = SRF( INDLOC )
           TRFVAL = TRF( INDLOC )
           IF ( QRFVAL.EQ.ZERO .AND. SRFVAL.EQ.ZERO .AND.
     1                 TRFVAL.EQ.ZERO    ) GOTO 70
C          .
C          . Calculate factors including L.
C          . Protect against division by zero.
C          .
           L      = ML( IH )
           IF ( L.LT.1 ) THEN
             PRINT *,' Subroutine RSDVER.'
             PRINT *,' ML(',IH,') = ', L
             PRINT *,' Program aborted.'
             STOP
           ENDIF
           DLFAC  = DBLE( L )
           DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
           SRFV2  = SRFVAL/DLSQR
           TRFV2  = TRFVAL/DLSQR
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
             PRINT *,' Subroutine RSDVER.'
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
               INDCOS         = 1
               FTF1( INDCOS ) = FTF1( INDCOS ) + DALF*QRFVAL
               FTF2( INDCOS ) = FTF2( INDCOS ) + DALFD*SRFV2
               FTF3( INDCOS ) = FTF3( INDCOS ) + DALFD*TRFV2
             ELSE
               INDCOS         = 2*MPS + 1
               INDSIN         = 2*MPS + 2
               TERM1          = DBLE( M )*DALF/SINE
               FTF1( INDCOS ) = FTF1( INDCOS ) + DALF*QRFVAL
               FTF2( INDCOS ) = FTF2( INDCOS ) + DALFD*SRFV2
               FTF2( INDSIN ) = FTF2( INDSIN ) + TERM1*TRFV2
               FTF3( INDSIN ) = FTF3( INDSIN ) - TERM1*SRFV2
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
             FTF1( INDSIN ) = FTF1( INDSIN ) + DALF*QRFVAL
             FTF2( INDSIN ) = FTF2( INDSIN ) + DALFD*SRFV2
             FTF2( INDCOS ) = FTF2( INDCOS ) - TERM1*TRFV2
             FTF3( INDSIN ) = FTF3( INDSIN ) + DALFD*TRFV2
             FTF3( INDCOS ) = FTF3( INDCOS ) + TERM1*SRFV2
C            .
           ENDIF
C          .
 70        CONTINUE
         ENDDO
C        .
C        . End loop ih = 1, nh
C
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

