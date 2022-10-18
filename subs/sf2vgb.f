C*********************************************************************
C subroutine Scalar Function 2 Vector Gradient Xtra special vector B *
C            -      -        - -      -        -                   - *
C Steve Gibbons Tue Feb 13 10:05:33 WET 2001                         C
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
C SF2VGB evaluates at grid nodes IR = ILNR to IRNR, the gradient of  C
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
C Requires the 3 arrays VGFA1( NH3, NTHP ), VGFA2( NH3, NTHP ) and   C
C VGFA3( NH3, NTHP ) which have been prepared by the routine SF2VGA. C
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
C     XARR      : Array length NR of radial node values.             C
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
C     VGFA1     : Dim ( NH3, NTHP ). Coeff.s from SF2VGA             C
C     VGFA2     : Dim ( NH3, NTHP ). Coeff.s from SF2VGA             C
C     VGFA3     : Dim ( NH3, NTHP ). Coeff.s from SF2VGA             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SF2VGB( ILNR, IRNR, NR, NH3, ML3, MM3, M0, NTHP,
     1                   NPHP, NCMX, ICMR, ICMT, ICMP, SV3, DSV3,
     2                   XARR, XSV, FTF1, FTF2, FTF3,
     3                   VGFA1, VGFA2, VGFA3 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          ILNR, IRNR, NR, NH3, ML3( NH3 ), ICMP,
     1                 MM3( NH3 ), M0, NTHP, NPHP, NCMX, ICMR, ICMT
      DOUBLE PRECISION SV3( * ), DSV3( * ), XARR( NR ),
     1                 XSV( NCMX, NPHP, NTHP, NR )
      DOUBLE PRECISION FTF1( 2*NPHP ), FTF2( 2*NPHP ), FTF3( 2*NPHP )
      DOUBLE PRECISION VGFA1( NH3, NTHP ),
     1                 VGFA2( NH3, NTHP ),
     2                 VGFA3( NH3, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ISIGN, ITHETA, LENV, IPHI, IRMNR,
     1                 L, M, INDSIN, INDCOS, INDLOC, IH3, IR
      DOUBLE PRECISION ZERO, DZCOEF, COEF, DCOEF, RAD
      PARAMETER ( ZERO = 0.0d0 )
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
        PRINT *,' Subroutine SF2VGB'
        PRINT *,' NCMX = ', NCMX
        PRINT *,' ICMR = ', ICMR
        PRINT *,' ICMT = ', ICMT
        PRINT *,' ICMP = ', ICMP
        PRINT *,' Program aborted,'
        STOP
      ENDIF
C
      ISIGN = -1
C
      LENV = 2*NPHP
C ............. start to loop around theta points ...................
      DO IR = ILNR, IRNR
       IRMNR = IR - NR
       RAD = XARR( IR )
C
C Uncomment following lines for this safety check:
c      IF ( DABS( RAD ).LT.TOL ) THEN
c        PRINT *,' Subroutine SF2VGB.'
c        PRINT *,' RAD = ', RAD,' Division by zero imminent.'
c        PRINT *,' Program stopped.'
c        STOP
c      ENDIF
C ............. start to loop around theta points ...................
       DO ITHETA = 1, NTHP
        DZCOEF = ZERO
C
C ............. set FTF1, FTF2, FTF3 all to zero
        CALL DVECZ( FTF1, LENV )
        CALL DVECZ( FTF2, LENV )
        CALL DVECZ( FTF3, LENV )
C
C ............. start to loop around Harmonics ......................
C
        INDLOC = IRMNR
        DO IH3 = 1, NH3 
           INDLOC = INDLOC + NR
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
           COEF   = SV3( INDLOC )/RAD
           DCOEF  = DSV3( INDLOC )
C
C Uncomment next line if you want this "check" - I suspect
C it merely wastes time ...
c          IF ( COEF.EQ.ZERO .AND. DCOEF.EQ.ZERO ) GOTO 50
C
           M      = MM3( IH3 )
C          .
           IF ( M.GE.0 ) THEN
C            .
C            . We are dealing with a cos ( m phi ) term
C            .
             IF ( M.EQ.0 ) THEN
               INDCOS         = 1
               FTF1( INDCOS ) = FTF1( INDCOS ) +
     1                   DCOEF*VGFA1( IH3, ITHETA )
               FTF2( INDCOS ) = FTF2( INDCOS ) +
     1                   COEF*VGFA2( IH3, ITHETA )
             ELSE
               INDCOS         = 1 + 2*M/M0
               INDSIN         = INDCOS + 1
               FTF1( INDCOS ) = FTF1( INDCOS ) +
     1                   DCOEF*VGFA1( IH3, ITHETA )
               FTF2( INDCOS ) = FTF2( INDCOS ) +
     1                   COEF*VGFA2( IH3, ITHETA )
               FTF3( INDSIN ) = FTF3( INDSIN ) +
     1                   COEF*VGFA3( IH3, ITHETA )
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS         = 1 - 2*M/M0
             INDSIN         = INDCOS + 1
             FTF1( INDSIN ) = FTF1( INDSIN ) +
     1                   DCOEF*VGFA1( IH3, ITHETA )
             FTF2( INDSIN ) = FTF2( INDSIN ) +
     1                   COEF*VGFA2( IH3, ITHETA )
             FTF3( INDCOS ) = FTF3( INDCOS ) +
     1                   COEF*VGFA3( IH3, ITHETA )
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
        INDLOC = -1
        DO IPHI = 1, NPHP
          INDLOC = INDLOC + 2
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
