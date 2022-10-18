C*********************************************************************
C subroutine Radial Spectrally Decomposed Vector 2 xtra special D ****
C            -      -          -          -      -              - ****
C Steve Gibbons Fri Feb 16 15:04:15 WET 2001                         C
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
C Before RSDV2D is called, the arrays FTFP, FTFT, INFP and INFT      C
C must first be calculated by a call to RSDV2C.                      C
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
C     INFP      : Dim(2,NPH). Locations in FTF2/3 calc. by RSDV2C    C
C                                                                    C
C     INFT      : Dim(2,NTH). Locations in FTF2/3 calc. by RSDV2C    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     QRF       : Dim ( NR*NPH ) Radial functions Q_{iph}(r)         C
C     SRF       : Dim ( NR*NPH ) Radial functions S_{iph}(r)         C
C     TRF       : Dim ( NR*NTH ) Radial functions T_{ith}(r)         C
C                                                                    C
C     FTFP     : Dim (3,NPH,NTHP). Coeff.s calculated by RSDV2C      C
C     FTFT     : Dim (2,NTH,NTHP). Coeff.s calculated by RSDV2C      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
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
      SUBROUTINE RSDV2D( NTHP, NPHP, NR, ILNR, IRNR, NPH,
     1                   NTH, NCMX, ICMR, ICMT, ICMP,
     2                   QRF, SRF, TRF, XSV, FTF1, FTF2, FTF3,
     3                   FTFP, FTFT, INFP, INFT, IFORMF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NTHP, NPHP, NR, ILNR, IRNR, NPH,
     1                 NTH, IFORMF, NCMX, ICMR, ICMT, ICMP
      INTEGER          INFP( 2, NPH ), INFT( 2, NTH )
      DOUBLE PRECISION QRF( * ), SRF( * ), TRF( * ),
     1                 XSV( NCMX, NPHP, NTHP, NR)
      DOUBLE PRECISION FTF1( 2*NPHP ), FTF2( 2*NPHP ),
     1                 FTF3( 2*NPHP ), FTFP( 3, NPH, NTHP ),
     2                 FTFT( 2, NTH, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER          ISIGN, ITHETA, LENV, IPHI, IR, 
     1                 INDLOC, IHP, IHT, I2, I3, INDINC
      DOUBLE PRECISION QRFVAL, SRFVAL, TRFVAL, DF1, DF2, DF3
C
      PARAMETER ( ISIGN = -1 )
c     DOUBLE PRECISION ZERO
c     PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C   
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
        PRINT *,' Subroutine RSDV2D'
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
        PRINT *,' Subroutine RSDV2D'
        PRINT *,' NCMX = ', NCMX
        PRINT *,' ICMR = ', ICMR
        PRINT *,' ICMT = ', ICMT
        PRINT *,' ICMP = ', ICMP
        PRINT *,' Program aborted,'
        STOP
      ENDIF
C   
      LENV = 2*NPHP
C     .................. start loop around radial grid nodes
      DO IR = ILNR, IRNR
C
C ............. start to loop around theta points ...................
        DO ITHETA = 1, NTHP
C
C ............. set all of FTF1, FTF2, FTF3 to zero
C
         CALL DVECZ( FTF1, LENV )
         CALL DVECZ( FTF2, LENV )
         CALL DVECZ( FTF3, LENV )
C
C ............. start to loop around Harmonics ......................
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
           QRFVAL = QRF( INDLOC )
           SRFVAL = SRF( INDLOC )
C          .
C          . Get indices for FFT arrays
C          .
           I2     = INFP( 1, IHP )
           I3     = INFP( 2, IHP )
C          .
           DF1    = FTFP( 1, IHP, ITHETA )
           DF2    = FTFP( 2, IHP, ITHETA )
           DF3    = FTFP( 3, IHP, ITHETA )
C          .
           FTF1( I2 ) = FTF1( I2 ) + DF1*QRFVAL
           FTF2( I2 ) = FTF2( I2 ) + DF2*SRFVAL
           FTF3( I3 ) = FTF3( I3 ) + DF3*SRFVAL
C          .
         ENDDO
C        .
C        . End loop ihp = 1, nhp
C
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
           TRFVAL = TRF( INDLOC )
C          .
C          . Get indices for FFT arrays
C          .
           I2     = INFT( 1, IHT )
           I3     = INFT( 2, IHT )
C          .
           DF2    = FTFT( 1, IHT, ITHETA )
           DF3    = FTFT( 2, IHT, ITHETA )
C          .
           FTF2( I2 ) = FTF2( I2 ) + DF2*TRFVAL
           FTF3( I3 ) = FTF3( I3 ) + DF3*TRFVAL
C          .
         ENDDO
C        .
C        . End loop iht = 1, nht
C        .
C ............. ended looping around Harmonics ......................
C
C ............... now perform Fourier Transforms on FTF1, FTF2, FTF3
C
         CALL FFTRLV( FTF1, NPHP, ISIGN )
         CALL FFTRLV( FTF2, NPHP, ISIGN )
         CALL FFTRLV( FTF3, NPHP, ISIGN )
C
C ...................................................................
         INDLOC = -1
         DO IPHI = 1, NPHP
           INDLOC = INDLOC + 2
           XSV( ICMR, IPHI, ITHETA, IR ) = FTF1( INDLOC )
           XSV( ICMT, IPHI, ITHETA, IR ) = FTF2( INDLOC )
           XSV( ICMP, IPHI, ITHETA, IR ) = FTF3( INDLOC )
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

