C*********************************************************************
C subroutine Scalar Function 2 Vector Gradient Xtra special vector ***
C            -      -        - -      -        -                   ***
C Steve Gibbons Thu Nov 30 13:53:45 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C A scalar function psi( r, \theta, \phi ) is decomposed             C
C into spherical harmonics such that                                 C
C                                                                    C
C  psi = \sum_{ih3 = 1, nh3} psi_{ih3}(r) Y_{ih3}                    C
C                                                                    C
C  Y_{ih3} = P_L^M( \cos \theta ) cos ( M \phi ) for MM3( ih3 ) = M  C
C     or                                                             C
C  Y_{ih3} = P_L^M( \cos \theta ) sin ( M \phi ) for MM3( ih3 ) = -M C
C                                                                    C
C    where L = ML3( ih3 ) and                                        C
C    the P_L^M are Schmidt quasi-normalised associated               C
C    Legendre functions (see routine SCHNLA ).                       C
C                                                                    C
C The radial functions psi_{ih3}(r) are stored in the vector SV3     C
C with psi_{ih3}( ir ) being stored in the element                   C
C                                                                    C
C  SV3( ( ih3 - 1 )*NR + IR )                                        C
C                                                                    C
C The radial derivative, d psi_{ih3}/dr at grid node ir is stored    C
C in the element    DSV3( ( ih3 - 1 )*NR + IR )                      C
C                                                                    C
C SF2VGX evaluates at grid nodes IR = ILNR to IRNR, the gradient of  C
C psi                                                                C
C                                                                    C
C with (\nabla \psi)_r        = d \psi / d r                         C
C with (\nabla \psi)_{\theta} = r^{-1} d \psi / d {\theta}           C
C with (\nabla \psi)_{\phi}   = {r sin theta}^{-1} d\psi / d{\phi}   C
C                                                                    C
C In the output array, XVF, the r, theta and phi components          C
C of (\nabla \psi) are stored in                                     C
C                                                                    C
C    XSV( ICMR, IPHP, ITHP, IR ),                                    C
C    XSV( ICMT, IPHP, ITHP, IR ) and                                 C
C    XSV( ICMP, IPHP, ITHP, IR ) respectively.                       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ILNR      : First radial grid node to be acted upon.           C
C     IRNR      : Final radial grid node to be acted upon.           C
C     NR        : Total number of radial grid nodes.                 C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NH3       : Number of scalar functions                         C
C     ML3       : Array dim ( NH3 ). Sph. harm degree, L.            C
C     MM3       : Array dim ( NH3 ). Sph. harm order, M, or -M.      C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     NTHP      : Number of theta points.                            C
C     NPHP      : Number of phi points.                              C
C     NCMX      : Maximum number of components stored in XSV         C
C     ICMR      : Index for radial component (see XSV)               C
C     ICMT      : Index for theta component (see XSV)                C
C     ICMP      : Index for phi component (see XSV)                  C
C                                                                    C
C ICMR, ICMT and ICMP must ofcourse be distinct and between 1 and    C
C NCMX.                                                              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV3       : Dim ( NR*NH3 ) Radial functions psi_{ip3}(r)       C
C     DSV3      : Dim ( NR*NH3 ) Radial derivatives d psi_{ip3}/dr   C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C                                                                    C
C     XARR      : Array length NR of radial node values.             C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     XSV       : eXtra Special Vector. An array of dimensions       C
C                ( NCMX, NPHP, NTHP, NR)                             C
C                 At the radial grid node, IR, and theta point       C
C                 ITHE and phi point IPHI, the radial component      C
C                 of the vector is stored in                         C
C                   XSV( ICMR, IPHI, ITHE, IR )                      C
C                 The theta component is stored in                   C
C                   XSV( ICMT, IPHI, ITHE, IR )                      C
C                 The phi component is stored in                     C
C                   XSV( ICMP, IPHI, ITHE, IR )                      C
C                                                                    C
C Note that XSV is not initialised on entry ...                      C
C                                                                    C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Have dimensions ( 2*NPHP )                        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SF2VGX( ILNR, IRNR, NR, LH, NH3, ML3, MM3, M0, NTHP,
     1                   NPHP, NCMX, ICMR, ICMT, ICMP, SV3, DSV3,
     2                   GAUX, XARR, PA, DPA, XSV, FTF1, FTF2, FTF3 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          ILNR, IRNR, NR, LH, NH3, ML3( NH3 ), ICMP,
     1                 MM3( NH3 ), M0, NTHP, NPHP, NCMX, ICMR, ICMT
      DOUBLE PRECISION SV3( * ), DSV3( * ), GAUX( NTHP ), XARR( NR ),
     1                 XSV( NCMX, NPHP, NTHP, NR )
      DOUBLE PRECISION FTF1( 2*NPHP ),
     1                 FTF2( 2*NPHP ), FTF3( 2*NPHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2, NTHP ),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IOP, ISIGN, ITHETA, LENV, IPHI,
     1                 IP, ICS, L, M, INDSIN, INDCOS, INDLOC, MPS,
     2                 IH3, IR
      DOUBLE PRECISION ZERO, X, SINE, TOL, DZCOEF, COEF, DCOEF,
     1                 DALF, DALFD, TERM1, RAD
      PARAMETER ( ZERO = 0.0d0, TOL = 1.0d-8 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C   Check the dimensions of XSV
C
      IF ( ICMR.EQ.ICMT .OR. ICMR.EQ.ICMP .OR. ICMT.EQ.ICMP .OR.
     1     ICMR.LT.1 .OR. ICMT.LT.1 .OR. ICMP.LT.1 .OR.
     2     ICMR.GT.NCMX .OR. ICMT.GT.NCMX .OR. ICMP.GT.NCMX  ) THEN
        PRINT *,' Subroutine SF2VGX'
        PRINT *,' NCMX = ', NCMX
        PRINT *,' ICMR = ', ICMR
        PRINT *,' ICMT = ', ICMT
        PRINT *,' ICMP = ', ICMP
        PRINT *,' Program aborted,'
        STOP
      ENDIF
C
C Pre-empt a division by zero with RAD
C
      IOP   = 0
      ISIGN = -1
C
C ............. start to loop around theta points ...................
      DO IR = ILNR, IRNR
       RAD = XARR( IR )
C
C Uncomment following lines for this safety check:
c      IF ( DABS( RAD ).LT.TOL ) THEN
c        PRINT *,' Subroutine SF2VGX.'
c        PRINT *,' RAD = ', RAD,' Division by zero imminent.'
c        PRINT *,' Program stopped.'
c        STOP
c      ENDIF
C ............. start to loop around theta points ...................
       DO ITHETA = 1, NTHP
        X      = GAUX( ITHETA )
        SINE   = DSQRT( 1.0d0 - X*X )
        DZCOEF = ZERO
C
C ............. set FTF1, FTF2, FTF3 all to zero
        LENV = 2*NPHP
        CALL VECOP ( FTF1, ZERO, LENV, IOP )
        CALL VECOP ( FTF2, ZERO, LENV, IOP )
        CALL VECOP ( FTF3, ZERO, LENV, IOP )
C
C ............. start to loop around Harmonics ......................
C
        DO IH3 = 1, NH3 
           INDLOC = IH3*NR - NR + IR
C
           L      = ML3( IH3 )
C
C First do monopole term
C
           IF ( L.EQ.0 ) THEN
             DZCOEF = DSV3( INDLOC )
             GOTO 50
           ENDIF
C
           COEF   = SV3( INDLOC )
           DCOEF  = DSV3( INDLOC )
C
           IF ( COEF.EQ.ZERO .AND. DCOEF.EQ.ZERO ) GOTO 50
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
C Uncomment following lines to include this safety check:
c          IF ( MPS*M0.NE.M ) THEN
c            PRINT *,' Subroutine SF2VGX.'
c            PRINT *,' MM3(',IH3,') = ', MM3( IH3 )
c            PRINT *,' Program aborted.'
c            STOP
c          ENDIF
C          .
           IF ( ICS.EQ.1 ) THEN
C            .
C            . We are dealing with a cos ( m phi ) term
C            .
             IF ( M.EQ.0 ) THEN
               INDCOS         = 1
               FTF1( INDCOS ) = FTF1( INDCOS ) + DALF*DCOEF
               FTF2( INDCOS ) = FTF2( INDCOS ) + DALFD*COEF/RAD
             ELSE
               INDCOS         = 2*MPS + 1
               INDSIN         = INDCOS + 1
c              INDSIN         = 2*MPS + 2
               TERM1          = DBLE( M )*DALF/(RAD*SINE)
               FTF1( INDCOS ) = FTF1( INDCOS ) + DALF*DCOEF
               FTF2( INDCOS ) = FTF2( INDCOS ) + COEF*DALFD/RAD
               FTF3( INDSIN ) = FTF3( INDSIN ) - COEF*TERM1
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS         = 2*MPS + 1
             INDSIN         = INDCOS + 1
c            INDSIN         = 2*MPS + 2
             TERM1          = DBLE( M )*DALF/(RAD*SINE)
             FTF1( INDSIN ) = FTF1( INDSIN ) + DALF*DCOEF
             FTF2( INDSIN ) = FTF2( INDSIN ) + COEF*DALFD/RAD
             FTF3( INDCOS ) = FTF3( INDCOS ) + COEF*TERM1
C            .
           ENDIF
C          .
 50     CONTINUE
        ENDDO
C ............. ended looping around Harmonics ......................
C 
C ............... now perform Fourier Transforms on FTF1, FTF2, FTF3
        CALL FFTRLV ( FTF1, NPHP, ISIGN )
        CALL FFTRLV ( FTF2, NPHP, ISIGN )
        CALL FFTRLV ( FTF3, NPHP, ISIGN )
C ...................................................................
        DO IPHI = 1, NPHP
          INDLOC = 2*IPHI - 1
          XSV( ICMR, IPHI, ITHETA, IR ) = FTF1( INDLOC ) + DZCOEF
          XSV( ICMT, IPHI, ITHETA, IR ) = FTF2( INDLOC )
          XSV( ICMP, IPHI, ITHETA, IR ) = FTF3( INDLOC )
        ENDDO
C
       ENDDO
C ............. ended looping around theta points ...................
      ENDDO
C ............. ended looping around radial grid nodes ..............
C
      RETURN
      END
C*********************************************************************
