C*********************************************************************
C subroutine Scalar Function 2 Vector Gradient Xtra special vector D *
C            -      -        - -      -        -                   - *
C Steve Gibbons Wed Feb 21 16:12:20 WET 2001                         C
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
C  SV3( ( ih3 - 1 )*NR + IR )     for IFORMF = 4        or           C
C  SV3( ( IR  - 1 )*NH3 + ih3 )   for IFORMF = 3                     C
C                                                                    C
C The radial derivative, d psi_{ih3}/dr at grid node ir is stored    C
C in the element    DSV3( ( ih3 - 1 )*NR + IR )  for IFORMF = 4 or   C
C                   DSV3( ( IR  - 1 )*NH3 + ih3) for IFORMF = 3      C
C                                                                    C
C SF2VGD evaluates at grid nodes IR = ILNR to IRNR, the gradient of  C
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
C Requires the 2 arrays VGFA( 3, NH3, NTHP ) and IVGFA( 2, NH3 )     C
C which have been prepared by the routine SF2VGC.                    C
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
C     IFORMF    : Defines the arrangement of the solution vec.s      C
C                   (see above: either set to 3 or 4)                C
C                                                                    C
C     IVGFA     : Dim ( 2, NH3 ). Locations from SF2VGD.             C
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
C     VGFA      : Dim ( 3, NH3, NTHP ). Coeff.s from SF2VGC          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SF2VGD( ILNR, IRNR, NR, NH3, ML3, NTHP,
     1                   NPHP, NCMX, ICMR, ICMT, ICMP, SV3, DSV3,
     2                   XARR, XSV, FTF1, FTF2, FTF3,
     3                   IVGFA, VGFA, IFORMF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          ILNR, IRNR, NR, NH3, ML3( NH3 ), ICMP, NTHP,
     1                 NPHP, NCMX, ICMR, ICMT, IVGFA( 2, NH3 ), IFORMF
      DOUBLE PRECISION SV3( * ), DSV3( * ), XARR( NR ),
     1                 XSV( NCMX, NPHP, NTHP, NR )
      DOUBLE PRECISION FTF1( 2*NPHP ), FTF2( 2*NPHP ), FTF3( 2*NPHP )
      DOUBLE PRECISION VGFA( 3, NH3, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ISIGN, ITHETA, LENV, IPHI,
     1                 L, INDLOC, IH3, IR, INDINC, I2, I3
      DOUBLE PRECISION ZERO, DZCOEF, COEF, DCOEF, RAD, DF1, DF2, DF3
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
        PRINT *,' Subroutine XSVSDD'
        PRINT *,' IFORMF = ', IFORMF
        PRINT *,' Program aborted,'
        STOP
      ENDIF
C
C   Check the dimensions of XSV
C
      IF ( ICMR.EQ.ICMT .OR. ICMR.EQ.ICMP .OR. ICMT.EQ.ICMP .OR.
     1     ICMR.LT.1 .OR. ICMT.LT.1 .OR. ICMP.LT.1 .OR.
     2     ICMR.GT.NCMX .OR. ICMT.GT.NCMX .OR. ICMP.GT.NCMX  ) THEN
        PRINT *,' Subroutine SF2VGD'
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
       RAD = XARR( IR )
C
C Uncomment following lines for this safety check:
c      IF ( DABS( RAD ).LT.TOL ) THEN
c        PRINT *,' Subroutine SF2VGD.'
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
         IF ( IFORMF.EQ.3 ) THEN
           INDLOC = IR*NH3 - NH3
           INDINC = 1
         ENDIF
C
         IF ( IFORMF.EQ.4 ) THEN
           INDLOC = IR - NR
           INDINC = NR
         ENDIF
C
        DO IH3 = 1, NH3 
           INDLOC = INDLOC + INDINC
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
C          .
C          . Get indices for FFT arrays
C          .
           I2     = IVGFA( 1, IH3 )
           I3     = IVGFA( 2, IH3 )
C          .
           DF1    = VGFA( 1, IH3, ITHETA )
           DF2    = VGFA( 2, IH3, ITHETA )
           DF3    = VGFA( 3, IH3, ITHETA )
C          .
           FTF1( I2 ) = FTF1( I2 ) + DF1*DCOEF
           FTF2( I2 ) = FTF2( I2 ) + DF2*COEF
           FTF3( I3 ) = FTF3( I3 ) + DF3*COEF
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
