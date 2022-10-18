C*********************************************************************
C subroutine Radial Spectrally Decomposed Vector 2 xtra special B ****
C            -      -          -          -      -              - ****
C Steve Gibbons Mon Feb 12 11:39:40 WET 2001                         C
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
C Before RSDV2B is called, the arrays FTF1P, FTF2P, FTF3P,           C
C FTF2T and FTF3T must first be calculated by a call to RSDV2A.      C
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
C     NR        : Number of radial grid nodes.                       C
C     ILNR      : Lowest radial node to be acted upon.               C
C     IRNR      : Highest radial node to be acted upon.              C
C     NPH       : Number of poloidal harmonics [number of radial     C
C                  functions Q(r) and S(r)]                          C
C     MMP       : Array dim ( NPH ). Sph. harm order, M, or -M.      C
C     NTH       : Number of toroidal harmonics [number of radial     C
C                  functions T(r)]                                   C
C     MMT       : Array dim ( NTH ). Sph. harm order, M, or -M.      C
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
C     QRF       : Dim ( NR*NPH ) Radial functions Q_{iph}(r)         C
C     SRF       : Dim ( NR*NPH ) Radial functions S_{iph}(r)         C
C     TRF       : Dim ( NR*NTH ) Radial functions T_{ith}(r)         C
C                                                                    C
C     FTF1P     : Dim (NPH,NTHP). Coeff.s calculated by RSDV2A       C
C     FTF2P     : Dim (NPH,NTHP). Coeff.s calculated by RSDV2A       C
C     FTF3P     : Dim (NPH,NTHP). Coeff.s calculated by RSDV2A       C
C     FTF2T     : Dim (NTH,NTHP). Coeff.s calculated by RSDV2A       C
C     FTF3T     : Dim (NTH,NTHP). Coeff.s calculated by RSDV2A       C
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
      SUBROUTINE RSDV2B( NTHP, NPHP, M0, NR, ILNR, IRNR, NPH,
     1                   MMP, NTH, MMT, NCMX, ICMR, ICMT, ICMP,
     2                   QRF, SRF, TRF, XSV, FTF1, FTF2, FTF3,
     3                   FTF1P, FTF2P, FTF3P, FTF2T, FTF3T )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NTHP, NPHP, M0, NR, ILNR, IRNR, NPH,
     1                 MMP( NPH ), NTH,
     2                 MMT( NTH ), NCMX, ICMR, ICMT, ICMP
      DOUBLE PRECISION QRF( * ), SRF( * ), TRF( * ),
     1                 XSV( NCMX, NPHP, NTHP, NR)
      DOUBLE PRECISION FTF1( 2*NPHP ), FTF2( 2*NPHP ),
     1                 FTF3( 2*NPHP ), FTF1P( NPH, NTHP ),
     2                 FTF2P( NPH, NTHP ), FTF3P( NPH, NTHP ),
     3                 FTF2T( NTH, NTHP ), FTF3T( NTH, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER          ISIGN, ITHETA, LENV, IPHI, IR, 
     1                 M, INDLOC, INDCOS, INDSIN, IHP, IHT, IRMNR
      DOUBLE PRECISION QRFVAL, SRFVAL, TRFVAL
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
C   Check the dimensions of XSV
C   
      IF ( ICMR.EQ.ICMT .OR. ICMR.EQ.ICMP .OR. ICMT.EQ.ICMP .OR.
     1     ICMR.LT.1 .OR. ICMT.LT.1 .OR. ICMP.LT.1 .OR. 
     2     ICMR.GT.NCMX .OR. ICMT.GT.NCMX .OR. ICMP.GT.NCMX  ) THEN
        PRINT *,' Subroutine RSDV2B'
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
        IRMNR = IR - NR
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
         INDLOC = IRMNR
         DO IHP = 1, NPH
           INDLOC = INDLOC + NR
           QRFVAL = QRF( INDLOC )
           SRFVAL = SRF( INDLOC )
C Uncomment the following line for this check:
C I suspect this serves no other purpose than
C slowing the program down!
c          IF ( QRFVAL.EQ.ZERO .AND. SRFVAL.EQ.ZERO ) GOTO 70
C          .
C          . Find wavenumber, m, (actually pseudo wavenumber)
C          .
           M      = MMP( IHP )
C          .
           IF ( M.GE.0 ) THEN
C            .
C            . We are dealing with a cos ( m phi ) term
C            .
             IF ( M.EQ.0 ) THEN
               INDCOS         = 1
               FTF1( INDCOS ) = FTF1( INDCOS ) +
     1                             FTF1P( IHP, ITHETA )*QRFVAL
               FTF2( INDCOS ) = FTF2( INDCOS ) +
     1                             FTF2P( IHP, ITHETA )*SRFVAL
             ELSE
               INDCOS         = 1 + 2*M/M0
               INDSIN         = INDCOS + 1
               FTF1( INDCOS ) = FTF1( INDCOS ) +
     1                             FTF1P( IHP, ITHETA )*QRFVAL
               FTF2( INDCOS ) = FTF2( INDCOS ) +
     1                             FTF2P( IHP, ITHETA )*SRFVAL
               FTF3( INDSIN ) = FTF3( INDSIN ) +
     1                             FTF3P( IHP, ITHETA )*SRFVAL
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
     1                             FTF1P( IHP, ITHETA )*QRFVAL
             FTF2( INDSIN ) = FTF2( INDSIN ) +
     1                             FTF2P( IHP, ITHETA )*SRFVAL
             FTF3( INDCOS ) = FTF3( INDCOS ) +
     1                             FTF3P( IHP, ITHETA )*SRFVAL
C            .
           ENDIF
C          .
 70        CONTINUE
         ENDDO
C        .
C        . End loop ihp = 1, nhp
C
C............. Now do toroidal harmonics (T radial func.s) .
         INDLOC = IRMNR
         DO IHT = 1, NTH
           INDLOC = INDLOC + NR
           TRFVAL = TRF( INDLOC )
C Uncomment the following line for check:
C I don't think this is ever the case, and so
C is a waste of time!
c          IF ( TRFVAL.EQ.ZERO ) GOTO 71
C          .
C          . Find wavenumber, m, (actually pseudo wavenumber)
C          .
           M      = MMT( IHT )
C          .
           IF ( M.GE.0 ) THEN
C            .
C            . We are dealing with a cos ( m phi ) term
C            .
             IF ( M.EQ.0 ) THEN
               INDCOS         = 1
               FTF3( INDCOS ) = FTF3( INDCOS ) +
     1                              FTF3T( IHT, ITHETA )*TRFVAL
             ELSE
               INDCOS         = 1 + 2*M/M0
               INDSIN         = INDCOS + 1
               FTF2( INDSIN ) = FTF2( INDSIN ) +
     1                              FTF2T( IHT, ITHETA )*TRFVAL
               FTF3( INDCOS ) = FTF3( INDCOS ) +
     1                              FTF3T( IHT, ITHETA )*TRFVAL
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS         = 1 - 2*M/M0
             INDSIN         = INDCOS + 1
             FTF2( INDCOS ) = FTF2( INDCOS ) +
     1                              FTF2T( IHT, ITHETA )*TRFVAL
             FTF3( INDSIN ) = FTF3( INDSIN ) +
     1                              FTF3T( IHT, ITHETA )*TRFVAL
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

