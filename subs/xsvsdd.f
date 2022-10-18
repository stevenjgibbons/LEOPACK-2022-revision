C*********************************************************************
C subroutine Xtra Special Vector function 2 Spectrally Decomp. vec D *
C            -    -       -                 -          -           - *
C Steve Gibbons Tue Feb 13 14:02:50 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C The double precision array XSV contains a vector function of       C
C dimensions ( NCNX, NPHP, NTHP, NR ). XSVSDD transforms this back   C
C into a spectral space in the following format.                     C
C                                                                    C
C The radial component of the vector is stored in                    C
C XSV( ICMR, iphi, ithe, ir )                                        C
C                                                                    C
C The theta component of the vector is stored in                     C
C XSV( ICMT, iphi, ithe, ir )                                        C
C                                                                    C
C The phi component of the vector is stored in                       C
C XSV( ICMP, iphi, ithe, ir )                                        C
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
C      QRF( ( ihp - 1 )*NR + IR )      for IFORMF = 4        or      C
C      QRF( ( IR  - 1 )*NPH + ihp )    for IFORMF = 3                C
C and similarly, S_{ihp}(ir) is stored in the element                C
C      SRF( ( ihp - 1 )*NR + IR )      for IFORMF = 4        or      C
C      SRF( ( IR  - 1 )*NPH + ihp )    for IFORMF = 3                C
C                                                                    C
C Likewise, the function T_{iht}(r) is stored in TRF with            C
C  T_{iht}(ir) in element                                            C
C          TRF( ( iht - 1 )*NR + IR )     for IFORMF = 4        or   C
C          TRF( ( IR  - 1 )*NTH + iht )   for IFORMF = 3             C
C                                                                    C
C Requires the four arrays                                           C
C                                                                    C
C       PFA( 3, NPH, NTHP ), JPFA( 2, NPH ),                         C
C       TFA( 2, NTH, NTHP ), JTFA( 2, NTH )                          C
C                                                                    C
C created by the routine XSVSDC.                                     C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NTHP      : Number of theta points.                            C
C     NPHP      : Number of phi points for perform Fourier transform C
C     NR        : Number of radial grid nodes.                       C
C     ILNR      : Lowest radial node to be acted upon.               C
C     IRNR      : Highest radial node to be acted upon.              C
C     NPH       : Number of poloidal harmonics [number of radial     C
C                  functions Q(r) and S(r)]                          C
C     NTH       : Number of toroidal harmonics [number of radial     C
C                  functions T(r)]                                   C
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
C     JPFA      : Dim(2,NPH). Locations in FTF2/3 calc. by XSVSDC    C
C     JTFA      : Dim(2,NTH). Locations in FTF2/3 calc. by XSVSDC    C
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
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHP )                         C
C                                                                    C
C     PFA       : Dim ( 3, NPH, NTHP ). Coeffs from XSVSDC.          C
C     TFA       : Dim ( 2, NTH, NTHP ). Coeffs from XSVSDC.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XSVSDD( NTHP, NPHP, NR, ILNR, IRNR, NPH,
     1                   NTH, NCMX, ICMR, ICMT, ICMP,
     2                   QRF, SRF, TRF, XSV, FTF1, FTF2, FTF3,
     3                   PFA, TFA, JPFA, JTFA, IFORMF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHP, NPHP, NR, ILNR, IRNR, NPH,
     1        NTH, NCMX, ICMR, ICMT, ICMP, IFORMF
      INTEGER JPFA( 2, NPH ), JTFA( 2, NTH )
      DOUBLE PRECISION QRF( * ), SRF( * ), TRF( * ),
     1                 XSV( NCMX, NPHP, NTHP, NR)
      DOUBLE PRECISION FTF1( 2*NPHP ),
     1                 FTF2( 2*NPHP ), FTF3( 2*NPHP ),
     2                 PFA( 3, NPH, NTHP ), TFA( 2, NTH, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ISIGN, ITHETA, IR, IHP, IHT, INDLOC,
     1        ILEN, IPHI, INDINC, I2, I3
      DOUBLE PRECISION QRFVAL, SRFVAL, TRFVAL, DF1, DF2, DF3
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
        PRINT *,' Subroutine XSVSDD'
        PRINT *,' NCMX = ', NCMX
        PRINT *,' ICMR = ', ICMR
        PRINT *,' ICMT = ', ICMT
        PRINT *,' ICMP = ', ICMP
        PRINT *,' Program aborted,'
        STOP
      ENDIF
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
      ILEN  = 2*NPHP
      ISIGN = 1
      DO IR = ILNR, IRNR
C ....................................................................
C ...... now start to loop around theta points .......................
       DO ITHETA = 1, NTHP
C        .
C        . Zero the arrays ftf1, ftf2 and ftf3
C        .
C note that ILEN = 2*NPHP
         CALL DVECZ( FTF1, ILEN )
         CALL DVECZ( FTF2, ILEN )
         CALL DVECZ( FTF3, ILEN )
C
C .................. firstly enter the RVF information into FTF1, FTF2,
C .................. FTF3 for r, theta, phi components respectively of
C .................. VF.
C
         INDLOC = -1
         DO IPHI = 1, NPHP
           INDLOC = INDLOC + 2
           FTF1( INDLOC ) = XSV( ICMR, IPHI, ITHETA, IR )
           FTF2( INDLOC ) = XSV( ICMT, IPHI, ITHETA, IR )
           FTF3( INDLOC ) = XSV( ICMP, IPHI, ITHETA, IR )
         ENDDO
C
C .................. now perform Forward Discreet Fourier Transforms
C .................. on FTF1, FTF2, FTF3
C
         CALL FFTRLV ( FTF1, NPHP, ISIGN )
         CALL FFTRLV ( FTF2, NPHP, ISIGN )
         CALL FFTRLV ( FTF3, NPHP, ISIGN )
C ...................................................................
C .                  Now let's loop around the Harmonics..          .
C ...................................................................
C ............. First do poloidal harmonics (Q and S radial func.s) .
C
         IF ( IFORMF.EQ.3 ) THEN
           INDLOC = IR*NPH - NPH
           INDINC = 1
         ENDIF
C
         IF ( IFORMF.EQ.4 ) THEN
           INDLOC = IR - NR
           INDINC = NR
         ENDIF
C
         DO IHP = 1, NPH
           INDLOC = INDLOC + INDINC
C          .
C          . Get indices for FFT arrays
C          .
           I2     = JPFA( 1, IHP )
           I3     = JPFA( 2, IHP )
C          .
           DF1    = PFA( 1, IHP, ITHETA )
           DF2    = PFA( 2, IHP, ITHETA )
           DF3    = PFA( 3, IHP, ITHETA )
C          .
           QRFVAL  = DF1*FTF1( I2 )
           SRFVAL  = DF2*FTF2( I2 ) + DF3*FTF3( I3 )
C          .
           QRF( INDLOC ) = QRF( INDLOC ) + QRFVAL
           SRF( INDLOC ) = SRF( INDLOC ) + SRFVAL
C          .
         ENDDO
C        .
C        . End loop ihp = 1, nhp
C        .
C............. Now do toroidal harmonics (T radial func.s) .
C
         IF ( IFORMF.EQ.3 ) THEN
           INDLOC = IR*NTH - NTH
           INDINC = 1
         ENDIF
C
         IF ( IFORMF.EQ.4 ) THEN
           INDLOC = IR - NR
           INDINC = NR
         ENDIF
C
         DO IHT = 1, NTH
           INDLOC = INDLOC + INDINC
C          .
C          . Get indices for FFT arrays
C          .
           I2     = JTFA( 1, IHT )
           I3     = JTFA( 2, IHT )
C          .
           DF2    = TFA( 1, IHT, ITHETA )
           DF3    = TFA( 2, IHT, ITHETA )
C          .
           TRFVAL  = DF3*FTF3( I3 ) + DF2*FTF2( I2 )
C          .
           TRF( INDLOC ) = TRF( INDLOC ) + TRFVAL
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
C ....................................................................
      ENDDO
C
C Ended looping around radial grid nodes ...
C
      RETURN
      END
C*********************************************************************
