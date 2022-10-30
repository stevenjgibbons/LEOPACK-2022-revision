C*********************************************************************
C                                                                    C
C Steve Gibbons - Uniform Boundary Thermal Convection Time Step      C
C Wed Nov 22 16:11:51 WET 2000                                       C
C                                                                    C
C This is an experimental development program to try and derive      C
C a subroutine to advance solution by a single time-step.            C
C                                                                    C
C O.k. - we now want to be able to output the kinetic energy         C
C at a single time-step.                                             C
C                                                                    C
C*********************************************************************
      PROGRAM o2ubtctsc2
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - contents of common blocks.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          NR, NH1, NH2, NH3, NBN, M0, IOUTF, LH,
     1                 ITMX, NOIT, NTHP, NPHP, NCMX
      COMMON  /NDIMPARS/  NR, NH1, NH2, NH3, NBN, M0, IOUTF, LH,
     1                    ITMX, NOIT, NTHP, NPHP, NCMX
      DOUBLE PRECISION CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CFAC, DELTAT, DTOL
      COMMON  /DPHYSPARS/ CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CFAC, DELTAT, DTOL
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          NRMAX, NH1MAX, NH2MAX, NH3MAX, NHMAX,
     1                 NIVMAX, NDCS, NCFM, NDRVM, NBNM, NIV1MX,
     2                 NIV2MX, NIV3MX, NCMXX, NNDM
      PARAMETER      ( NRMAX = 50, NH1MAX =  520, NH2MAX =  520,
     1                 NH3MAX =  520, NCMXX = 13, NNDM = 4,
     2                 NHMAX = NH1MAX+NH2MAX+NH3MAX,
     3                 NIVMAX = NHMAX*NRMAX )
      PARAMETER      ( NIV1MX = NH1MAX*NRMAX,
     1                 NIV2MX = NH2MAX*NRMAX,
     2                 NIV3MX = NH3MAX*NRMAX )
      PARAMETER      ( NDCS = 3, NBNM = 3, NCFM = 2*NBNM+1,
     1                 NDRVM = 4 )
      INTEGER          LHMAX, NTHMAX, NPHMAX, NPMAX
      PARAMETER      ( LHMAX =  36, NTHMAX = 100, NPHMAX =  64,
     1                 NPMAX = (LHMAX+1)*(LHMAX+2)/2 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     1                 MHP( NHMAX ), LARR( NDCS ), MHIBC( NDCS ),
     2                 MHOBC( NDCS ), IWORK( NCFM )
      INTEGER          ML1( NH1MAX ), MM1( NH1MAX ), MP1( NH1MAX )
      INTEGER          ML2( NH2MAX ), MM2( NH2MAX ), MP2( NH2MAX )
      INTEGER          ML3( NH3MAX ), MM3( NH3MAX ), MP3( NH3MAX )
      DOUBLE PRECISION SV( NIVMAX ), XARR( NRMAX ),
     1                 SVFDC( NCFM, NRMAX, NDRVM+1, NDCS ),
     2                 COEFM1( NCFM, NCFM ), W1( NCFM ),
     3                 COEFM2( NCFM, NCFM ), W2( NCFM ),
     4                 FDCM( NCFM, NRMAX, 1 )
      DOUBLE PRECISION SV1( NIV1MX ), SV2( NIV2MX ), SV3( NIV3MX )
      DOUBLE PRECISION VE1( NIV1MX ), VE2( NIV2MX ), VE3( NIV3MX )
      DOUBLE PRECISION DV1( NIV1MX ), DV2( NIV2MX ), DV3( NIV3MX )
      DOUBLE PRECISION DSV3( NIV3MX ), VQ1( NIV1MX ), VS1( NIV1MX ),
     1                 VT2( NIV2MX ), VT1( NIV1MX ), VQ2( NIV2MX ),
     2                 VS2( NIV2MX )
      DOUBLE PRECISION R1( NIV1MX ), R2( NIV2MX ), R3( NIV3MX )
C
C matrices for time-stepping
C
      DOUBLE PRECISION AM1( 3*NBNM+1, NIV1MX ),
     1                 BM1( 2*NBNM+1, NIV1MX )
      DOUBLE PRECISION AM2( 3*NBNM+1, NIV2MX ),
     1                 BM2( 2*NBNM+1, NIV2MX )
      DOUBLE PRECISION AM3( 3*NBNM+1, NIV3MX ),
     1                 BM3( 2*NBNM+1, NIV3MX )
      INTEGER          IPIV1( NIV1MX ), IPIV2( NIV2MX ),
     1                 IPIV3( NIV3MX )
C
C matrix for taking first derivatives of temperature fnctns.
C
      DOUBLE PRECISION SV3D( 2*NBNM+1, NIV3MX )
C
C matrices for forming QST decomposition of velocity
C
      DOUBLE PRECISION SV1Q(        1, NIV1MX )
      DOUBLE PRECISION SV1S( 2*NBNM+1, NIV1MX )
      DOUBLE PRECISION SV2T(        1, NIV2MX )
C
C matrices for taking curls of vel. QST decomposition
C
      DOUBLE PRECISION CVQ1(        1, NIV1MX )
      DOUBLE PRECISION CVS1( 2*NBNM+1, NIV1MX )
      DOUBLE PRECISION CVT2Q(        1, NIV2MX )
      DOUBLE PRECISION CVT2S( 2*NBNM+1, NIV2MX )
C
      DOUBLE PRECISION CQ1T(        1, NIV1MX )
      DOUBLE PRECISION CS1T( 2*NBNM+1, NIV1MX )
      DOUBLE PRECISION CT2P(        1, NIV2MX )
C
C Arrays for spherical transforms
C
      DOUBLE PRECISION GAUX( NTHMAX ), GAUW( NTHMAX ),
     1                 PA( NPMAX, NTHMAX ), DPA( NPMAX, NTHMAX )
      DOUBLE PRECISION XSV( NCMXX, NPHMAX, NTHMAX, NRMAX )
      DOUBLE PRECISION FTF1( 2*NPHMAX ), FTF2( 2*NPHMAX ),
     1                 FTF3( 2*NPHMAX )
      INTEGER          IW( NNDM )
      DOUBLE PRECISION DW1( NNDM ), DW2( NNDM ), DDW( NNDM, NNDM )
C
C Transform coefficients
C
C Velocity spectral --> real space
C
      INTEGER          INFPV( 2, NH1MAX ),
     1                 INFTV( 2, NH2MAX )
C
      DOUBLE PRECISION FTFPV( 3, NH1MAX, NTHMAX ),
     1                 FTFTV( 2, NH2MAX, NTHMAX )
C
C Curl velocity spectral --> real space
C
      INTEGER          INFPCV( 2, NH2MAX ),
     1                 INFTCV( 2, NH1MAX )
C
      DOUBLE PRECISION FTFPCV( 3, NH2MAX, NTHMAX ),
     1                 FTFTCV( 2, NH1MAX, NTHMAX )
C
C Grad. temp
C
      INTEGER          IVGFA( 2, NH3MAX )
      DOUBLE PRECISION VGFA( 3, NH3MAX, NTHMAX )
C
C Scalar real --> spectral transform
C
      INTEGER          ISF2S( NH3MAX )
      DOUBLE PRECISION SF2SA( NH3MAX, NTHMAX )
C
C Calculated by XSVSDC
C
      INTEGER          JPFAV( 2, NH1MAX ), JTFAV( 2, NH2MAX )
C
      DOUBLE PRECISION PFAV( 3, NH1MAX, NTHMAX ),
     1                 TFAV( 2, NH2MAX, NTHMAX )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER           LU, LULOG, I, ILENRT, ILEN, IWRTE, IFORM,
     1                  NCUDS, INARR( 3 ), NH, NR1, NDRVS, ILN, IRN,
     2                  ICOMP, IDIR, MMAX, IOPT
      INTEGER           ITS, NTS, NTSBB, NTSBS, N, INC, NDIG, N1, N2,
     1                  K, NTSBFE, LUEVAL
      DOUBLE PRECISION  PVLC( NBNM ), PVRC( NBNM ), DMONE, DPONE,
     1                  STIME, DTIME, RAD, THE, PHI, RTPFCE, PI, 
     2                  PVKE, TVKE, DKE( 2 ), VRAD, VPHI, TEMP, FAC
      EXTERNAL          RTPFCE
      CHARACTER *(120)  LINE, FNLOG, ROOT, FNAME, FNEVAL
      CHARACTER *(10)   CHNUM
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IWRTE = 3, DMONE = -1.0d0, DPONE = 1.0d0,
     1            INC = 1, IFORM = 1, PI=3.14159265358979312D0 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NCMX    = NCMXX
      NBN     = NBNM
      LU      = 11
      LULOG   = 12
      LUEVAL  = 13
 80   FORMAT(A)
C
      DO I = 1, 120
        LINE(I:I) = ' '
      ENDDO
C
      PRINT *,' Enter stem for output files.'
 200  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 200
C
C Extract ROOT from LINE
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILENRT = I - 1
          GOTO 199
        ENDIF
      ENDDO
 199  CONTINUE
      ROOT(1:ILENRT)  = LINE(1:ILENRT)
      FNLOG(1:ILENRT) = ROOT(1:ILENRT)
      FNLOG(ILENRT+1:ILENRT+4) = '.log'
      FNLOG = FNLOG(1:ILENRT+4)
      CALL FOPEN( LULOG, FNLOG, IWRTE )
C
C filename for evaluation file
C
      FNEVAL(1:ILENRT) = ROOT(1:ILENRT)
      FNEVAL(ILENRT+1:ILENRT+5) = '.eval'
      FNEVAL = FNEVAL(1:ILENRT+5)
      CALL FOPEN( LUEVAL, FNEVAL, IWRTE )
C
      PRINT *,' Enter name of harmonics file.'
 201  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 201
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 202
        ENDIF
      ENDDO
 202  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in harmonics file
C
      NCUDS = 0
C     (ncuds - number of diff. schemes already in use).
      CALL HMFRD( NH, NHMAX, MHT, MHL, MHM, MHP, NCUDS, NDCS,
     1            MHIBC, MHOBC, LARR, LU, FNAME )
      INARR( 3 ) = NH
      WRITE ( LULOG, * ) 'PROGRAM ubtctsc.'
      WRITE ( LULOG, * ) 
      WRITE ( LULOG, * ) 'Harmonics file read:'
      WRITE ( LULOG, * ) FNAME
      WRITE ( LULOG, * ) 'Total number of harmonics    = ',NH
      WRITE ( LULOG, * ) 'Number of difference schemes = ',NCUDS
      WRITE ( LULOG, * ) 
C
      DO I = 1, NCUDS
        WRITE ( LULOG, 203 ) I, MHIBC( I ), MHOBC( I )
      ENDDO
 203  FORMAT('Scheme(',I4,') IBC = ',I3,' OBC = ',I3)
      WRITE ( LULOG, * ) 
C
      PRINT *,' Enter name of vector file.'
 204  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 204
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 205
        ENDIF
      ENDDO
 205  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in vector file
C
      CALL SVFRD( INARR, LU, NRMAX, SV, FNAME )
      WRITE ( LULOG, * ) 'Vector file read:'
      WRITE ( LULOG, * ) FNAME
      NR1 = INARR( 2 )
      WRITE ( LULOG, * ) 'with ',NR1,' radial grid nodes.'
      WRITE ( LULOG, * ) 
C
      PRINT *,' Enter name of radial spacing file.'
 206  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 206
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 207
        ENDIF
      ENDDO
 207  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in radial spacing file
C
      CALL XARRRD( NR, NRMAX, XARR, LU, FNAME )
C
      WRITE ( LULOG, * ) 'Radial spacings file read:'
      WRITE ( LULOG, * ) FNAME
      WRITE ( LULOG, * )
C
      IF ( NR1.NE.NR ) THEN
        WRITE ( LULOG, * ) 'Solution vector and radial node '
        WRITE ( LULOG, * ) 'file claim differing numbers of '
        WRITE ( LULOG, * ) 'grid nodes. Program aborted.'
        GOTO 999
      ENDIF
C
C Calculate finite difference coefficients
C
      NDRVS = 1
      ILN   = 1
      IRN   = 1
      CALL SVFDCF( NR, NDCS, NBN, ILN, IRN, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, W1, W2 )
C
      NDRVS = 2
      ILN   = 2
      IRN   = 2
      CALL SVFDCF( NR, NDCS, NBN, ILN, IRN, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, W1, W2 )
C
      NDRVS = 4
      ILN   = 3
      IRN   = NR - 2
      CALL SVFDCF( NR, NDCS, NBN, ILN, IRN, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, W1, W2 )
C
      NDRVS = 2
      ILN   = NR - 1
      IRN   = NR - 1
      CALL SVFDCF( NR, NDCS, NBN, ILN, IRN, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, W1, W2 )
C
      NDRVS = 1
      ILN   = NR
      IRN   = NR
      CALL SVFDCF( NR, NDCS, NBN, ILN, IRN, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, W1, W2 )
C
      WRITE ( LULOG, * ) 'SVFDC calculated.'
      WRITE ( LULOG, * )
C
C Calculate the non-boundary specific f.d. coefficients
C
      NDRVS = 1
      ILN   = 2
      IRN   = NR - 1
      CALL FDCMBD( NR, NBN, ILN, IRN, ILN, IRN, NCFM,
     1             NCFM, NDRVS, NDRVS, IWORK, XARR, FDCM,
     2             COEFM1, W1, W2 )
C
      WRITE ( LULOG, * ) 'FDCM calculated.'
      WRITE ( LULOG, * )
C
C Now split harmonic sets into individual components
C First put poloidal velocity harmonics into SV1
C
      ICOMP = 1
      CALL IIASCE( NH, ICOMP, MHT, MHL, MHM, MHP, NH1, NH1MAX,
     1             ML1, MM1, MP1 )
      WRITE ( LULOG, * ) 'ML1, MM1 and MP1 filled.'
      WRITE ( LULOG, * ) 'There are ',NH1,' pol. vel. harms.'
      WRITE ( LULOG, * )
C
      IDIR  = 1
      ILN   = 3
      IRN   = NR - 2
      CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1             NH1, ML1, MM1, SV, SV1, IDIR )
C
      WRITE ( LULOG, * ) 'Initial SV1 vector filled.'
      WRITE ( LULOG, * )
C
C Calculate the completion coefficients for pol. vel.
C
      CALL PVCCF( NR, NDCS, NBN, NCFM, NDRVM, MP1(1), SVFDC,
     1            PVLC, PVRC )
      WRITE ( LULOG, * ) 'PVLC and PVRC calculated.'
      WRITE ( LULOG, * )
C
C Now put toroidal velocity harmonics into SV2
C
      ICOMP = 2
      CALL IIASCE( NH, ICOMP, MHT, MHL, MHM, MHP, NH2, NH2MAX,
     1             ML2, MM2, MP2 )
      WRITE ( LULOG, * ) 'ML2, MM2 and MP2 filled.'
      WRITE ( LULOG, * ) 'There are ',NH2,' tor. vel. harms.'
      WRITE ( LULOG, * )
C
      IDIR  = 1
      ILN   = 2
      IRN   = NR - 1
      CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1             NH2, ML2, MM2, SV, SV2, IDIR )
C
      WRITE ( LULOG, * ) 'Initial SV2 vector filled.'
      WRITE ( LULOG, * )
C
C Now put temperature harmonics into SV3
C
      ICOMP = 3
      CALL IIASCE( NH, ICOMP, MHT, MHL, MHM, MHP, NH3, NH3MAX,
     1             ML3, MM3, MP3 )
      WRITE ( LULOG, * ) 'ML3, MM3 and MP3 filled.'
      WRITE ( LULOG, * ) 'There are ',NH3,' temperature harms.'
      WRITE ( LULOG, * )
C
      IDIR  = 1
      ILN   = 2
      IRN   = NR - 1
      CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1             NH3, ML3, MM3, SV, SV3, IDIR )
C
      WRITE ( LULOG, * ) 'Initial SV3 vector filled.'
      WRITE ( LULOG, * )
C
C Calculate M0, MMAX and LH
C (We ignore negative MM1, MM2 and MM3 as each -m will
C have a corresponding +m)
C
      LH = 0
      DO I = 1, NH1
        IF ( ML1( I ).GT.LH ) LH = ML1( I )
      ENDDO
      DO I = 1, NH2
        IF ( ML2( I ).GT.LH ) LH = ML2( I )
      ENDDO
      DO I = 1, NH3
        IF ( ML3( I ).GT.LH ) LH = ML3( I )
      ENDDO
      WRITE ( LULOG, * ) 'LH (highest sph. harm. deg. = ', LH
C
      IF ( LH.GT.LHMAX ) THEN
        PRINT *,' LH = ', LH,' LHMAX = ', LHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      M0   = LH
      MMAX = 0
      DO I = 1, NH1
        IF ( MM1( I ).GT.0 .AND. MM1( I ).LT.M0 ) M0 = MM1( I )
        IF ( MM1( I ).GT.MMAX ) MMAX = MM1( I )
      ENDDO
      DO I = 1, NH2
        IF ( MM2( I ).GT.0 .AND. MM2( I ).LT.M0 ) M0 = MM2( I )
        IF ( MM2( I ).GT.MMAX ) MMAX = MM2( I )
      ENDDO
      DO I = 1, NH3
        IF ( MM3( I ).GT.0 .AND. MM3( I ).LT.M0 ) M0 = MM3( I )
        IF ( MM3( I ).GT.MMAX ) MMAX = MM3( I )
      ENDDO
      WRITE ( LULOG, * ) 'M0 (lowest non zero wavenumber. = ', M0
C
      PRINT *,' Enter number of points for theta quadrature.'
      PRINT *,' Must be atleast LH + 1. '
      PRINT *,' Type -1 for the minimum number and type '
      PRINT *,' -2 for LH+LH/2 (to avoid aliassing). '
C
 208  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 208
      READ ( LINE, * ) NTHP
C
C set NTHP to minimum value
C
      IF (  NTHP.EQ.-1  ) THEN
        NTHP = LH + 1
        IF ( NTHP/2*2.NE.NTHP ) NTHP = NTHP + 1
      ENDIF
C
C set NTHP to anti-aliassing value
C
      IF (  NTHP.EQ.-2  ) THEN
        NTHP = LH + LH/2
        IF ( NTHP/2*2.NE.NTHP ) NTHP = NTHP + 1
      ENDIF
C
C Validate if NTHP is large enough or within bounds
C
      IF ( NTHP.GE.(LH+1) .AND. NTHP.LE.NTHMAX ) GOTO 209
C
C O.K. so NTHP entered was invalid.
C
      WRITE ( LULOG, * ) 'Value NTHP = ',NTHP,' is invalid.'
      WRITE ( LULOG, * ) 'Program aborted.'
      GOTO 999
C
 209  CONTINUE
      WRITE ( LULOG, * ) 'NTHP (number of theta points) = ', NTHP
      WRITE ( LULOG, * )
C
C Calculate NPHP:
C
      CALL NPHPF( MMAX, M0, NPHP, NPHMAX )
      WRITE ( LULOG, * ) 'NPHP (number of phi points) = ', NPHP
      WRITE ( LULOG, * ) 'MMAX (largest wavenumber)   = ', MMAX
      WRITE ( LULOG, * )
C
C Calculate Gauss points and weights
C
      CALL GAUWTS ( DMONE, DPONE, GAUX, GAUW, NTHP )
      WRITE ( LULOG, * ) 'GAUWTS called.'
C
C Calculate arrays of Associated Legendre Functions
C
      CALL SCHNLA ( PA, DPA, GAUX, LH, NTHP )
      WRITE ( LULOG, * ) 'SCHNLA called.'
      WRITE ( LULOG, * )
C
C Read in physical parameters: First heat equation
C
      PRINT *,' Enter CA, CB1, CB2, CC, CD '
 210  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 210
      READ ( LINE, * ) CA, CB1, CB2, CC, CD
      WRITE ( LULOG, * ) 'CA  = ', CA
      WRITE ( LULOG, * ) 'CB1 = ', CB1
      WRITE ( LULOG, * ) 'CB2 = ', CB2
      WRITE ( LULOG, * ) 'CC  = ', CC
      WRITE ( LULOG, * ) 'CD  = ', CD
      WRITE ( LULOG, * )
C
C Read in physical parameters: vorticity equation
C
      PRINT *,' Enter CE, CF, CG, CH, CI '
 211  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 211
      READ ( LINE, * ) CE, CF, CG, CH, CI
      WRITE ( LULOG, * ) 'CE = ', CE
      WRITE ( LULOG, * ) 'CF = ', CF
      WRITE ( LULOG, * ) 'CG = ', CG
      WRITE ( LULOG, * ) 'CH = ', CH
      WRITE ( LULOG, * ) 'CI = ', CI
      WRITE ( LULOG, * )
C
C Enter CFAC and DELTAT
C
      PRINT *,' Enter CFAC, DELTAT, STIME.'
 212  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 212
      READ ( LINE, * ) CFAC, DELTAT, STIME
      WRITE ( LULOG, * ) 'CFAC   = ', CFAC
      WRITE ( LULOG, * ) 'DELTAT = ', DELTAT
      WRITE ( LULOG, * ) 'STIME  = ', STIME
      WRITE ( LULOG, * )
C
C No need to check value of CFAC, DELTAT etc. as this
C is done by PVTSMF, TVTSMF and TMTSMF
C      
      CALL PVTSMF( NR, NH1, NDCS, NBN, MHIBC, MHOBC, NCFM,
     1             NDRVM, ML1, MP1, IPIV1, XARR, SVFDC,
     2             AM1, BM1, CE, CI, CFAC, DELTAT )
      WRITE ( LULOG, * ) 'Called PVTSMF.'
C
      CALL TVTSMF( NR, NH2, NDCS, NBN, MHIBC, MHOBC, NCFM,
     1             NDRVM, ML2, MP2, IPIV2, XARR, SVFDC,
     2             AM2, BM2, CE, CI, CFAC, DELTAT )
      WRITE ( LULOG, * ) 'Called TVTSMF.'
C
      CALL TMTSMF( NR, NH3, NDCS, NBN, MHIBC, MHOBC, NCFM,
     1             NDRVM, ML3, MP3, IPIV3, XARR, SVFDC,
     2             AM3, BM3, CA, CD, CFAC, DELTAT )
      WRITE ( LULOG, * ) 'Called TMTSMF.'
      WRITE ( LULOG, * )
C
C Now need to calculate some auxiliary matrices which
C assist with the evaluation of the non-linear terms.
C First, we will build the array SV3D( 1+2*NBN, NH3*NR )
C which will multiply the solution vector SV3 to give the
C derivative vector DSV3. Use AVBMBR.
C
      N1     = 1+2*NBN
      N2     = NH3*NR
      IOPT   = 10
      ILN    = 2
      IRN    = NR-1
      K      = NBN
      CALL AVBMBR( N1, N2, K, NR, NH3, NBN, NCFM, NDRVM, ILN,
     1             IRN, ML3, MP3, IOPT, NDCS, XARR, SV3D, SVFDC )
      WRITE ( LULOG, * ) 'Formed SV3D.'
C
C Now we will build the array SV1Q( 1, NH1*NR )
C which will multiply the solution vector SV1 to give
C a vector of scaloidal radial functions. Use
C VOBMBR as there is no derivative.
C
      N1     = 1
      N2     = NH1*NR
      IOPT   = 1
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH1, NBN, NCFM, ILN, IRN,
     1             ML1, IOPT, XARR, SV1Q, FDCM )
      WRITE ( LULOG, * ) 'Formed SV1Q.'
C
C Now we will build the array SV1S( 1+2*NBN, NH1*NR )
C which will multiply the solution vector SV1 to give
C a vector of spheroidal radial functions. Use
C AVBMBR for increased accuracy.
C
      N1     = 1+2*NBN
      N2     = NH1*NR
      IOPT   = 2
      ILN    = 2
      IRN    = NR-1
      K      = NBN
      CALL AVBMBR( N1, N2, K, NR, NH1, NBN, NCFM, NDRVM, ILN,
     1             IRN, ML1, MP1, IOPT, NDCS, XARR, SV1S, SVFDC )
      WRITE ( LULOG, * ) 'Formed SV1S.'
C
C Now we form SV2T( 1, NH2*NR )
C which will multiply the solution vector SV2 to give
C a vector of toroidal radial functions. (T as opposed
C to tau - see my PhD thesis). Use
C VOBMBR as there is no derivative.
C
      N1     = 1
      N2     = NH2*NR
      IOPT   = 3
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH2, NBN, NCFM, ILN, IRN,
     1             ML2, IOPT, XARR, SV2T, FDCM )
      WRITE ( LULOG, * ) 'Formed SV2T.'
C
C We now form CVQ1( 1, NH1*NR )
C which will multiply VQ1 to give VT1 -
C a vector of T(toroidal) radial functions.
C
      N1     = 1
      N2     = NH1*NR
      IOPT   = 4
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH1, NBN, NCFM, ILN, IRN,
     1             ML1, IOPT, XARR, CVQ1, FDCM )
      WRITE ( LULOG, * ) 'Formed CVQ1.'
C
C We now form CVS1( 1+2*NBN, NH1*NR )
C which will multiply VS1 to add to VT1 -
C a vector of T(toroidal) radial functions.
C
      N1     = 1+2*NBN
      N2     = NH1*NR
      IOPT   = 5
      K      = NBN
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH1, NBN, NCFM, ILN, IRN,
     1             ML1, IOPT, XARR, CVS1, FDCM )
      WRITE ( LULOG, * ) 'Formed CVS1.'
C
C We now form CVT2Q( 1, NH2*NR )
C which will multiply VT2 to form VQ2 -
C a vector of scaloidal radial functions.
C
      N1     = 1
      N2     = NH2*NR
      IOPT   = 6
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH2, NBN, NCFM, ILN, IRN,
     1             ML2, IOPT, XARR, CVT2Q, FDCM )
      WRITE ( LULOG, * ) 'Formed CVT2Q.'
C
C We now form CVT2S( 1+2*NBN, NH2*NR )
C which will multiply VT2 to form VS2 -
C a vector of spheroidal radial functions.
C
      N1     = 1+2*NBN
      N2     = NH2*NR
      IOPT   = 7
      K      = NBN
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH2, NBN, NCFM, ILN, IRN,
     1             ML2, IOPT, XARR, CVT2S, FDCM )
      WRITE ( LULOG, * ) 'Formed CVT2S.'
C
C We now form CQ1T( 1, NH1*NR )
C which will multiply VQ1 to form R1 -
C a vector of toroidal (tau) radial functions.
C
      N1     = 1
      N2     = NH1*NR
      IOPT   = 11
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH1, NBN, NCFM, ILN, IRN,
     1             ML1, IOPT, XARR, CQ1T, FDCM )
      WRITE ( LULOG, * ) 'Formed CQ1T.'
C
C We now form CS1T( 1+2*NBN, NH1*NR )
C which will multiply VS1 to form R1 -
C a vector of toroidal (tau) radial functions.
C
      N1     = 1+2*NBN
      N2     = NH1*NR
      IOPT   = 12
      K      = NBN
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH1, NBN, NCFM, ILN, IRN,
     1             ML1, IOPT, XARR, CS1T, FDCM )
      WRITE ( LULOG, * ) 'Formed CS1T.'
C
C We now form CT2P( 1, NH2*NR )
C which will multiply VT2 to form R2 -
C a vector of poloidal radial functions.
C
      N1     = 1
      N2     = NH2*NR
      IOPT   = 13
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH2, NBN, NCFM, ILN, IRN,
     1             ML2, IOPT, XARR, CT2P, FDCM )
      WRITE ( LULOG, * ) 'Formed CT2P.'
C
C Calculate transform coefficients:
C Firstly the four arrays FTFPV( 3, NH1, NTHP ),
C FTFTV( 2, NH2, NTHP ), INFPV( 2, NH1 ) and INFTV( 2, NH2 )
C which transform the velocity into REAL space.
C
      CALL RSDV2C( NTHP, M0, LH, NH1, ML1,
     1             MM1, NH2, ML2, MM2, GAUX, PA, DPA,
     2             FTFPV, FTFTV, INFPV, INFTV )
      WRITE ( LULOG, * ) 'Called RSDV2C. (Velocity)'
C
C Now, the four arrays FTFPCV( 3, NH2, NTHP ),
C FTFTCV( 2, NH1, NTHP ), INFPCV( 2, NH2 ) and INFTCV( 2, NH1 )
C which transform the velocity curl into REAL space.
C
      CALL RSDV2C( NTHP, M0, LH, NH2, ML2,
     1             MM2, NH1, ML1, MM1, GAUX, PA, DPA,
     2             FTFPCV, FTFTCV, INFPCV, INFTCV )
      WRITE ( LULOG, * ) 'Called RSDV2C. (vel. curl)'
C
C Now, the 2 arrays VGFA( 3, NH3, NTHP )
C and IVGFA( NH3, NTHP ) which transform the gradient
C of the the temperature into real space.
C
      CALL SF2VGC( LH, NH3, ML3, MM3, M0, NTHP,
     1             GAUX, PA, DPA, VGFA, IVGFA )
      WRITE ( LULOG, * ) 'Called SF2VGC. (grad. theta).'
C
C Now the array SF2SA( NH3, NTHP ) which allows us
C to transform back the scalar function into spectral ...
C
      FAC  = -1.0d0*CC
      CALL SF2SDC( LH, NH3, ML3, MM3, M0, NTHP, NPHP, PA, GAUW,
     1             FAC, SF2SA, ISF2S )
      WRITE ( LULOG, * ) 'Called SF2SDC. (scalar trans.)'
C
C Now the arrays PFAV( 3, NH1, NTHP ), TFAV( 2, NH2, NTHP ),
C JPFAV( 2, NH1 ) and JTFAV( 2, NH2 )
C
      CALL XSVSDC( NTHP, M0, LH, NH1, ML1, MM1, NH2, ML2, MM2,
     1             GAUX, GAUW, PA, DPA, PFAV, TFAV, JPFAV, JTFAV )
C
      WRITE ( LULOG, * ) 'Called XSVSDC. (momentum eqn.)'
C
C Obtain information output flag along
C with number of time-steps, number of time-steps
C between backup of solution, and number of time-steps
C between solution snapshots.
C
      WRITE ( LULOG, * )
C
      PRINT *,' Enter IOUTF, NTS, NTSBB, NTSBS.'
      PRINT *,' ------------------------------ '
      PRINT *,' ioutf = 0 for no output,'
      PRINT *,' ioutf = 1 for output to .log file and'
      PRINT *,' ioutf = 6 for output to screen.'
      PRINT *,' nts is the total number of time-steps taken.'
      PRINT *,' ntsbb is no. of time-steps between backups.'
      PRINT *,' ntsbs is no. of time-steps between snapshots.'
 213  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 213
      READ ( LINE, * ) IOUTF, NTS, NTSBB, NTSBS
      IF (  IOUTF.NE.0 .AND. IOUTF.NE.1 .AND. IOUTF.NE.6 ) THEN
        WRITE ( LULOG, * ) ' IOUTF = ', IOUTF,' is invalid.'
        WRITE ( LULOG, * ) 'Program aborted.'
        GOTO 999
      ENDIF
      IF ( IOUTF.EQ.1 ) IOUTF = LULOG
      WRITE ( LULOG, * ) 'IOUTF  = ', IOUTF
      WRITE ( LULOG, * ) 'NTS    = ', NTS
      WRITE ( LULOG, * ) 'NTSBB  = ', NTSBB
      WRITE ( LULOG, * ) 'NTSBS  = ', NTSBS
      WRITE ( LULOG, * )
C
C Enter the criteria for the semi-implicit time-step
C procedure
C
      PRINT *,' Please enter ITMX, DTOL '
      PRINT *,' ITMX = max. number of iterations allowed'
      PRINT *,' per time-step.'
      PRINT *,' DTOL is how close the norms of solutions'
      PRINT *,' on successive iterations need to be to'
      PRINT *,' achieve convergence.'
 214  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 214
      READ ( LINE, * )  ITMX, DTOL
      WRITE ( LULOG, * ) 'ITMX   = ', ITMX
      WRITE ( LULOG, * ) 'DTOL   = ', DTOL
      WRITE ( LULOG, * )
C
C Enter number of time-steps between energy evaluations
C
      PRINT *,' Please enter NTSBFE.'
      PRINT *,' This is the number of time-steps between'
      PRINT *,' energy evaluations.'
 215  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 215
      READ ( LINE, * )  NTSBFE
      WRITE ( LULOG, * ) 'NTSBFE = ', NTSBFE
      WRITE ( LULOG, * )
C
C Fix the values where flow components are evaluated
C
      RAD = 0.5d0*( XARR( 1 ) + XARR( NR ) )
      THE = 0.5d0*PI
      PHI = 0.0d0
C
C Begin the time-stepping procedure
C
      ITS    = 0
      DTIME   = STIME
 50   CONTINUE
      ITS    = ITS + 1
      DTIME   = DTIME + DELTAT
C
C Call main time-step routine
C Current solution is in SV1, SV2 and SV3
C
      CALL UBTCT1( XARR,ML1,MM1,IPIV1,ML2,MM2,IPIV2,ML3,MM3,IPIV3,AM1,
     1  BM1,AM2,BM2,AM3,BM3,SV1,SV2,SV3,DV1,DV2,DV3,VE1,VE2,VE3,R1,R2,
     2  R3,PVLC,PVRC,SV3D,SV1Q,SV1S,SV2T,DSV3,VQ1,VS1,VT2,PA,DPA,FTF1,
     3  FTF2,FTF3,XSV,GAUX,GAUW,CVQ1,CVS1,CVT2Q,CVT2S,VT1,VQ2,VS2,
     4  CQ1T,CS1T,CT2P,FTFPV,FTFTV,INFPV,INFTV,
     5  FTFPCV,FTFTCV,INFPCV,INFTCV,VGFA,IVGFA,SF2SA,ISF2S,
     6  PFAV,TFAV,JPFAV,JTFAV)
C 
      IF ( NOIT.EQ.-1 ) THEN
        WRITE ( LULOG, * ) ' Time-step = ', ITS
        WRITE ( LULOG, * ) ' Maximum iterations exceeded.'
        GOTO 999
      ENDIF
C
      IF ( NOIT.EQ.-2 ) THEN
        WRITE ( LULOG, * ) ' Time-step = ', ITS
        WRITE ( LULOG, * ) ' Iteration not converged.'
        GOTO 999
      ENDIF
C
C New solution is stored in VE1, VE2 and VE3
C So copy them into SV1, SV2 and SV3
C
      N = NR*NH1
      CALL DCOPY( N, VE1, INC, SV1, INC )
C
      N = NR*NH2
      CALL DCOPY( N, VE2, INC, SV2, INC )
C
      N = NR*NH3
      CALL DCOPY( N, VE3, INC, SV3, INC )
C
C See if we wish to perform a solution write out?
C
      IF ( ITS.EQ.NTS .OR. ITS/NTSBB*NTSBB.EQ.ITS .OR.
     1     ITS/NTSBS*NTSBS.EQ.ITS  ) THEN
C       .
C       . Decide on the name of the file
C       .
        IF ( ITS/NTSBS*NTSBS.EQ.ITS .OR. ITS.EQ.NTS ) THEN
C         .
C         . We append '.ts__.sv' onto ROOT
C         .
          IF ( ITS.GE.1 .AND. ITS.LT.10 ) THEN
            NDIG = 1
            WRITE ( CHNUM(1:NDIG), 301 ) ITS
          ENDIF
          IF ( ITS.GE.10 .AND. ITS.LT.100 ) THEN
            NDIG = 2
            WRITE ( CHNUM(1:NDIG), 302 ) ITS
          ENDIF
          IF ( ITS.GE.100 .AND. ITS.LT.1000 ) THEN
            NDIG = 3
            WRITE ( CHNUM(1:NDIG), 303 ) ITS
          ENDIF
          IF ( ITS.GE.1000 .AND. ITS.LT.10000 ) THEN
            NDIG = 4
            WRITE ( CHNUM(1:NDIG), 304 ) ITS
          ENDIF
          IF ( ITS.GE.10000 .AND. ITS.LT.100000 ) THEN
            NDIG = 5
            WRITE ( CHNUM(1:NDIG), 305 ) ITS
          ENDIF
          IF ( ITS.GE.100000 .AND. ITS.LT.1000000 ) THEN
            NDIG = 6
            WRITE ( CHNUM(1:NDIG), 306 ) ITS
          ENDIF
          IF ( ITS.GE.1000000 .AND. ITS.LT.10000000 ) THEN
            NDIG = 7
            WRITE ( CHNUM(1:NDIG), 307 ) ITS
          ENDIF
C         .
          FNAME(1:ILENRT)  = ROOT(1:ILENRT)
          FNAME(ILENRT+1:ILENRT+3)  = '.ts'
          FNAME(ILENRT+4:ILENRT+3+NDIG)  = CHNUM(1:NDIG)
          FNAME(ILENRT+4+NDIG:ILENRT+6+NDIG) = '.sv'
          FNAME = FNAME(1:ILENRT+6+NDIG)
C         .
 301      FORMAT(I1)
 302      FORMAT(I2)
 303      FORMAT(I3)
 304      FORMAT(I4)
 305      FORMAT(I5)
 306      FORMAT(I6)
 307      FORMAT(I7)
C         .
        ELSE
C         .
C         . We append '.bv' onto ROOT
C         .
          FNAME(1:ILENRT)  = ROOT(1:ILENRT)
          FNAME(ILENRT+1:ILENRT+3)  = '.bv'
          FNAME = FNAME(1:ILENRT+3)
C         .
        ENDIF
C       .
C       . We have chosen output file name.
C       . Now we must collate our solution vectors.
C       .
        IDIR   = -1
C       .
        ICOMP  = 1
        ILN    = 3
        IRN    = NR - 2
        CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1               NH1, ML1, MM1, SV, SV1, IDIR )
C       .
        ICOMP  = 2
        ILN    = 2
        IRN    = NR - 1
        CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1               NH2, ML2, MM2, SV, SV2, IDIR )
C       .
        ICOMP  = 3
        ILN    = 2
        IRN    = NR - 1
        CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1               NH3, ML3, MM3, SV, SV3, IDIR )
C       .
C       . SV now contains the new solution vector,
C       . but with the boundary values missing -
C       . So let's take this opportunity to complete it.
C       .
        CALL ASVCPL( SV, NR, NDCS, INARR, MHP, MHIBC, MHOBC,
     1               NCFM, NDRVM, NDRVM, NBN, SVFDC )
C       .
        CALL SVFWT( INARR, LU, IFORM, SV, FNAME )
C       .
      ENDIF
C
C See if we wish to perform an evaluation
C
      IF ( ITS/NTSBFE*NTSBFE.EQ.ITS ) THEN
C       .
C       . Now we must collate our solution vectors.
C       .
        IDIR   = -1
C       .
        ICOMP  = 1
        ILN    = 3
        IRN    = NR - 2
        CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1               NH1, ML1, MM1, SV, SV1, IDIR )
C       .
        ICOMP  = 2
        ILN    = 2
        IRN    = NR - 1
        CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1               NH2, ML2, MM2, SV, SV2, IDIR )
C       .
        ICOMP  = 3
        ILN    = 2
        IRN    = NR - 1
        CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1               NH3, ML3, MM3, SV, SV3, IDIR )
C       .
C       . SV now contains the new solution vector,
C       . but with the boundary values missing -
C       . So let's take this opportunity to complete it.
C       .
        CALL ASVCPL( SV, NR, NDCS, INARR, MHP, MHIBC, MHOBC,
     1               NCFM, NDRVM, NDRVM, NBN, SVFDC )
C       .
        PVKE   = 0.0d0
        TVKE   = 0.0d0
C       .
        DO I = 1, NH
C         .
          CALL SHKEER( I, NDCS, NR, INARR, MHT, MHL, MHP, NBN,
     1               NDRVM, NDRVM, NCFM, SV, XARR, DKE, SVFDC )
C         .
          IF ( MHT( I ).EQ.1 ) PVKE = PVKE + DKE( 1 )
          IF ( MHT( I ).EQ.2 ) TVKE = TVKE + DKE( 1 )
C         .
        ENDDO
C       .
C       . Evaluate radial velocity component
C       .
        ICOMP = 1
        VRAD  = RTPFCE( ICOMP, INARR, NNDM, IW, MHT, MHL, MHM, RAD,
     1                 THE, PHI, SV, XARR, DW1, DW2, DDW )
C       .
C       . Evaluate phi velocity component
C       .
        ICOMP = 3
        VPHI  = RTPFCE( ICOMP, INARR, NNDM, IW, MHT, MHL, MHM, RAD,
     1                 THE, PHI, SV, XARR, DW1, DW2, DDW )
C       .
C       . Evaluate temperature
C       .
        ICOMP = 7
        TEMP  = RTPFCE( ICOMP, INARR, NNDM, IW, MHT, MHL, MHM, RAD,
     1                 THE, PHI, SV, XARR, DW1, DW2, DDW )
C       .
C       . Write to LUEVAL, the values DTIME, PVKE, TVKE,
C       . VRAD, VPHI, TEMP
C       .
        WRITE ( LUEVAL, 109 ) DTIME, PVKE, TVKE, VRAD, VPHI, TEMP
 109    FORMAT (6(1PD16.7))
C       .
      ENDIF
C
C Return to perform next time-step
C
      IF ( ITS.EQ.NTS ) GOTO 999
      CALL FLUSH( LUEVAL )
      CALL FLUSH( LULOG  )
      GOTO 50
C
 999  CONTINUE
      CALL FCLOSE( LULOG, FNLOG, 'Error' )
      CALL FCLOSE( LUEVAL, FNEVAL, 'Error' )
      STOP
      END
C*********************************************************************

C*********************************************************************
C Subroutine Uniform Boundary Thermal Convection Time step 1 *********
C            -       -        -       -          -         - *********
C Steve Gibbons Wed Nov 22 16:18:57 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C We are trying to time-step the heat and momentum equations only.   C
C The formulation necessary is:-                                     C
C                                                                    C
C  The heat equation ...                                             C
C                                                                    C
C    CA d Theta / dt  =   CD Lap Theta                               C
C                       + CC v . Grad Theta                          C
C                       + CB1 vrad * R                               C
C                       + CB2 vrad * R^{-3}                          C
C                                                                    C
C where R is the radial distance, vrad is the radial component of    C
C the velocity vector, v, and Theta is the temperature perturbation  C
C from the profile imposed by CB1 and CB2.                           C
C                                                                    C
C  The vorticity equation is the curl of ...                         C
C                                                                    C
C    CE d v / dt      =   CI Lap v                                   C
C                       + CH Theta RADVEC                            C
C                       - CG ( k x v )                               C
C                       + CF [ v x (curl v) ]                        C
C                                                                    C
C  where RADVEC is the vector ( R, 0, 0 ) and k is the vector        C
C  ( cos theta, -sin theta, 0 ).                                     C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Common Block                                                      C
C  ------------                                                      C
C                                                                    C
C   /NDIMPARS/                                                       C
C                                                                    C
C  Contains the integer numbers NR, NH1, NH2, NH3, NBN, M0, IOUTF,   C
C                               LH, ITMX, NOIT, NTHP, NPHP, NCMX     C
C                                                                    C
C  where                                                             C
C                                                                    C
C     NR        : Total number of radial grid nodes.                 C
C     NH1       : Number of poloidal velocity harmonics in soln.     C
C     NH2       : Number of toroidal velocity harmonics in soln.     C
C     NH3       : Number of temperature harmonics in soln.           C
C     NBN       : Number of bounding nodes for finite differences.   C
C     M0        : Lowest non-zero wavenumber in solution.            C
C     IOUTF     : Logical Unit number of output info file.           C
C                 (If IOUTF = 0 then no output is written).          C
C     LH        : Highest spherical harmonic degree, l.              C
C     ITMX      : Maximum number of iterations allowed per           C
C                  time-step                                         C
C     NOIT      : Actual number of iterations taken.                 C
C                  If number of iterations required exceeds ITMX     C
C                  then NOIT is returned as -1.                      C
C                  If soln. norm appears to be increasing after      C
C                  3 iterations, NOIT is returned as -2.             C
C     NTHP      : Number of theta points                             C
C     NPHP      : Number of phi points                               C
C     NCMX      : Leading dimension of XSV array.                    C
C                                                                    C
C   /DPHYSPARS/                                                      C
C                                                                    C
C  Physical parameters for problem. All double precision scalars     C
C                                                                    C
C     CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI: see above descrptn.  C
C     of equations.                                                  C
C     CFAC      : Determines how implicit/explicit integration is.   C
C                 Must be strictly greater than 0.0 and strictly     C
C                 less than 1.0  The higher CFAC is, the more        C
C                 explicit the integration. CFAC = 0.5 --> Crank-    C
C                 Nicolson scheme.                                   C
C     DELTAT    : Time-step size.                                    C
C     DTOL      : How close the solution norms on consecutive        C
C                 iterations need to be.                             C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ML1       : Array dim ( NH1 ). Sph. harm degree, L.            C
C     MM1       : Array dim ( NH1 ). Sph. harm order, M, or -M.      C
C            ( ml1 and mm1 describe poloidal velocity )              C
C     IPIV1     : Dim (NH1*NR). Pivotting information for AM1        C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C     ML2       : Array dim ( NH2 ). Sph. harm degree, L.            C
C     MM2       : Array dim ( NH2 ). Sph. harm order, M, or -M.      C
C            ( ml2 and mm2 describe toroidal velocity )              C
C     IPIV2     : Dim (NH2*NR). Pivotting information for AM2        C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C     ML3       : Array dim ( NH3 ). Sph. harm degree, L.            C
C     MM3       : Array dim ( NH3 ). Sph. harm order, M, or -M.      C
C            ( ml3 and mm3 describe temperature )                    C
C     IPIV3     : Dim (NH3*NR). Pivotting information for AM3        C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C     INFPV     : Dim ( 2, NH1 ) {from RSDV2C}                       C
C     INFTV     : Dim ( 2, NH2 ) {from RSDV2C}                       C
C                                                                    C
C     INFPCV    : Dim ( 2, NH2 ) {from RSDV2C}                       C
C     INFTCV    : Dim ( 2, NH1 ) {from RSDV2C}                       C
C                                                                    C
C     IVGFA     : Dim ( 2, NH3 ). Locations from SF2VGC              C
C                                                                    C
C     ISF2S     : Dim ( NH3 ). Indices from SF2SDC                   C
C                                                                    C
C     JPFAV     : Dim(2,NH1). Locations in FTF2/3 calc. by XSVSDC    C
C     JTFAV     : Dim(2,NH2). Locations in FTF2/3 calc. by XSVSDC    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives node radius r_i.  C
C                                                                    C
C     AM1       : Matrix for solution of f^{i+1}. (poloidal vel.)    C
C                 Has dimensions ( 3*NBN + 1, NH1*NR )               C
C                                                                    C
C     BM1       : Matrix for multiplication of f^i. (poloidal vel.)  C
C                 Has dimensions ( 2*NBN + 1, NH1*NR )               C
C                                                                    C
C     AM2       : Matrix for solution of f^{i+1}. (toroidal vel.)    C
C                 Has dimensions ( 3*NBN + 1, NH2*NR )               C
C                                                                    C
C     BM2       : Matrix for multiplication of f^i. (toroidal vel.)  C
C                 Has dimensions ( 2*NBN + 1, NH2*NR )               C
C                                                                    C
C     AM3       : Matrix for solution of f^{i+1}. (temperature)      C
C                 Has dimensions ( 3*NBN + 1, NH3*NR )               C
C                                                                    C
C     BM3       : Matrix for multiplication of f^i. (temperature)    C
C                 Has dimensions ( 2*NBN + 1, NH3*NR )               C
C                                                                    C
C     SV1       : Dim ( NH1*NR ). Poloidal velocity initial soln.    C
C     SV2       : Dim ( NH2*NR ). Toroidal velocity initial soln.    C
C     SV3       : Dim ( NH3*NR ). Temperature initial soln.          C
C                                                                    C
C     VE1       : Dim ( NH1*NR ). Poloidal velocity final soln.      C
C     VE2       : Dim ( NH2*NR ). Toroidal velocity final soln.      C
C     VE3       : Dim ( NH3*NR ). Temperature final soln.            C
C                                                                    C
C     DV1       : Dim ( NH1*NR ). Work array.                        C
C     DV2       : Dim ( NH2*NR ). Work array.                        C
C     DV3       : Dim ( NH3*NR ). Work array.                        C
C                                                                    C
C     R1        : Dim ( NH1*NR ). Work array.                        C
C     R2        : Dim ( NH2*NR ). Work array.                        C
C     R3        : Dim ( NH3*NR ). Work array.                        C
C                                                                    C
C     PVLC      : Dim ( NBN ). Formed by PVCCF.                      C
C     PVRC      : Dim ( NBN ). Formed by PVCCF.                      C
C                                                                    C
C     SV3D      : Dim ( 2*NBN + 1, NH3*NR ) takes 1st deriv. of SV3  C
C     SV1Q      : Dim (         1, NH1*NR ) turns p. vel. to Q( r )  C
C     SV1S      : Dim ( 2*NBN + 1, NH1*NR ) turns p. vel. to S( r )  C
C     SV2T      : Dim (         1, NH2*NR ) turns t. vel. to T( r )  C
C                                                                    C
C     CVQ1      : Dim (         1, NH1*NR ) takes curl of Q( r )     C
C     CVS1      : Dim ( 2*NBN + 1, NH1*NR ) takes curl of S( r )     C
C     CVT2Q     : Dim (         1, NH2*NR ) takes curl of T( r )     C
C     CVT2S     : Dim ( 2*NBN + 1, NH2*NR ) takes curl of T( r )     C
C                                                                    C
C     CQ1T      : Dim (         1, NH1*NR ) takes curl of Q( r )     C
C     CS1T      : Dim ( 2*NBN + 1, NH1*NR ) takes curl of S( r )     C
C     CT2P      : Dim (         1, NH2*NR ) takes curl of T( r )     C
C                                                                    C
C     FTF1      : Dim (2*NPHP). Work array for fourier transforming. C
C     FTF2      : Dim (2*NPHP). Work array for fourier transforming. C
C     FTF3      : Dim (2*NPHP). Work array for fourier transforming. C
C                                                                    C
C     DSV3      : Dim ( NH3*NR ). Work array.                        C
C     VQ1       : Dim ( NH1*NR ). Work array.                        C
C     VS1       : Dim ( NH1*NR ). Work array.                        C
C     VT2       : Dim ( NH2*NR ). Work array.                        C
C     VT1       : Dim ( NH1*NR ). Work array.                        C
C     VQ2       : Dim ( NH2*NR ). Work array.                        C
C     VS2       : Dim ( NH2*NR ). Work array.                        C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                  Both PA and DPA are formed by the routine SCHNLA. C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Corresponding Gauss quadrature weights. Dim (NTHP) C
C                                                                    C
C     XSV       : Work array. Dim ( NCMX, NPHP, NTHP, NR )           C
C                                                                    C
C     FTFPV     : Transform coeffs Dim ( 3, NH1, NTHP ) {from RSDV2C}C
C     FTFTV     : Transform coeffs Dim ( 2, NH2, NTHP ) {from RSDV2C}C
C                                                                    C
C     FTFPCV    : Transform coeffs Dim ( 3, NH2, NTHP ) {from RSDV2C}C
C     FTFTCV    : Transform coeffs Dim ( 2, NH1, NTHP ) {from RSDV2C}C
C                                                                    C
C     VGFA      : Dim ( 3, NH3, NTHP ). Coeff.s from SF2VGC          C
C                                                                    C
C     SF2SA     : Dim ( NH3, NTHP ). Coeff.s from SF2SDC             C
C                                                                    C
C     PFAV      : Dim ( 3, NH1, NTHP ). Coeffs from XSVSDC.          C
C     TFAV      : Dim ( 2, NH2, NTHP ). Coeffs from XSVSDC.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE UBTCT1( XARR,ML1,MM1,IPIV1,ML2,MM2,IPIV2,ML3,MM3,
     1  IPIV3,AM1,BM1,AM2,BM2,AM3,BM3,SV1,SV2,SV3,DV1,DV2,DV3,VE1,VE2,
     2  VE3,R1,R2,R3,PVLC,PVRC,SV3D,SV1Q,SV1S,SV2T,DSV3,VQ1,VS1,VT2,
     3  PA,DPA,FTF1,FTF2,FTF3,XSV,GAUX,GAUW,CVQ1,CVS1,CVT2Q,CVT2S,
     4  VT1,VQ2,VS2,CQ1T,CS1T,CT2P,FTFPV,FTFTV,INFPV,INFTV,
     5  FTFPCV,FTFTCV,INFPCV,INFTCV,VGFA,IVGFA,SF2SA,ISF2S,
     6  PFAV,TFAV,JPFAV,JTFAV)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Common block contents .....................C
      INTEGER          NR, NH1, NH2, NH3, NBN, M0, IOUTF, LH,
     1                 ITMX, NOIT, NTHP, NPHP, NCMX
      COMMON  /NDIMPARS/  NR, NH1, NH2, NH3, NBN, M0, IOUTF, LH,
     1                    ITMX, NOIT, NTHP, NPHP, NCMX
      DOUBLE PRECISION CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CFAC, DELTAT, DTOL
      COMMON  /DPHYSPARS/ CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CFAC, DELTAT, DTOL
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION XARR( NR )
C
C Basic solution vector indices
C
      INTEGER          ML1( NH1 ), MM1( NH1 ), IPIV1( NH1*NR )
      INTEGER          ML2( NH2 ), MM2( NH2 ), IPIV2( NH2*NR )
      INTEGER          ML3( NH3 ), MM3( NH3 ), IPIV3( NH3*NR )
C
C Diffusion/time-step matrices
C
      DOUBLE PRECISION AM1( 3*NBN+1, NH1*NR ), BM1( 2*NBN+1, NH1*NR )
      DOUBLE PRECISION AM2( 3*NBN+1, NH2*NR ), BM2( 2*NBN+1, NH2*NR )
      DOUBLE PRECISION AM3( 3*NBN+1, NH3*NR ), BM3( 2*NBN+1, NH3*NR )
C
C Auxiliary matrices
C
      DOUBLE PRECISION SV3D( 2*NBN + 1, NH3*NR )
      DOUBLE PRECISION SV1Q(         1, NH1*NR )
      DOUBLE PRECISION SV1S( 2*NBN + 1, NH1*NR )
      DOUBLE PRECISION SV2T(         1, NH2*NR )
C
      DOUBLE PRECISION CVQ1(          1, NH1*NR )
      DOUBLE PRECISION CVS1(  2*NBN + 1, NH1*NR )
      DOUBLE PRECISION CVT2Q(         1, NH2*NR )
      DOUBLE PRECISION CVT2S( 2*NBN + 1, NH2*NR )
C
      DOUBLE PRECISION CQ1T(       1, NH1*NR )
      DOUBLE PRECISION CS1T( 2*NBN+1, NH1*NR )
      DOUBLE PRECISION CT2P(       1, NH2*NR )
C
C Initial solution vectors
C
      DOUBLE PRECISION SV1( NH1*NR ), SV2( NH2*NR ), SV3( NH3*NR )
C
C Final solution vectors
C
      DOUBLE PRECISION VE1( NH1*NR ), VE2( NH2*NR ), VE3( NH3*NR )
C
C Arrays necessary for transform
C
      DOUBLE PRECISION FTF1(2*NPHP), FTF2(2*NPHP), FTF3(2*NPHP),
     1                 GAUX( NTHP ), GAUW( NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP ),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP )
C
C Work arrays
C
      DOUBLE PRECISION DV1( NH1*NR ), DV2( NH2*NR ), DV3( NH3*NR )
      DOUBLE PRECISION R1( NH1*NR ),  R2( NH2*NR ),  R3( NH3*NR )
      DOUBLE PRECISION DSV3( NH3*NR ), VQ1( NH1*NR ), VS1( NH1*NR ),
     1                 VT2( NH2*NR ), XSV( NCMX, NPHP, NTHP, NR),
     2                 VT1( NH1*NR ), VQ2( NH2*NR ), VS2( NH2*NR )
C
C Poloidal velocity completeness arrays
C
      DOUBLE PRECISION PVLC( NBN ), PVRC( NBN )
C
C Transform coefficients
C (pre-calculated by RSDV2C)
C
      DOUBLE PRECISION FTFPV( 3, NH1, NTHP ),
     1                 FTFTV( 2, NH2, NTHP )
C
      INTEGER          INFPV( 2, NH1 ),
     1                 INFTV( 2, NH2 )
C
      DOUBLE PRECISION FTFPCV( 3, NH2, NTHP ),
     1                 FTFTCV( 2, NH1, NTHP )
C
      INTEGER          INFPCV( 2, NH2 ),
     1                 INFTCV( 2, NH1 )
C
C Pre-calculated by SF2VGC
C
      INTEGER          IVGFA( 2, NH3 )
      DOUBLE PRECISION VGFA( 3, NH3, NTHP )
C
C Pre-calculated by SF2SDC
C
      INTEGER          ISF2S( NH3 )
      DOUBLE PRECISION SF2SA( NH3, NTHP )
C
C Calculated by XSVSDC
C
      INTEGER          JPFAV( 2, NH1 ), JTFAV( 2, NH2 )
      DOUBLE PRECISION PFAV( 3, NH1, NTHP ), TFAV( 2, NH2, NTHP )
C
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          KL, KU, LDA, M, N, INCX, INCY, NRHS, INFO,
     1                 LLU
      LOGICAL          LW
      CHARACTER *(1)   TRANS
      DOUBLE PRECISION ALPHA, BETA, ODIFF, DIFF, AVEC1, AOLD,
     1                 DNRM2, FAC
      EXTERNAL         DNRM2
C
      PARAMETER ( TRANS = 'N', INCX   = 1, INCY   = 1, NRHS = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Odiff should begin very large
C
      ODIFF = 1.0d8
C
      KL = NBN
      KU = NBN
C
      IF ( IOUTF.EQ.0 ) THEN
        LW  = .FALSE.
        LLU = 6
      ELSE
        LW  = .TRUE.
        LLU = IOUTF
      ENDIF
C
      IF ( LW ) WRITE ( IOUTF, * ) 'Entered UBTCT1.'
C
C Calculate the contribution from the diffusive parts.
C Poloidal velocity:
C   Multiply banded matrix BM1 by vector SV1 to give DV1
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = KL + KU + 1
      M      = NR*NH1
      N      = M
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, BM1, LDA, SV1, INCX,
     1             BETA, DV1, INCY )
      IF ( LW ) WRITE ( IOUTF, * ) 
     1    'UBTCT1: Called DGBMV: dv1:= bm1 . sv1.'
C
C Toroidal velocity:
C   Multiply banded matrix BM2 by vector SV2 to give DV2
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = KL + KU + 1
      M      = NR*NH2
      N      = M
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, BM2, LDA, SV2, INCX,
     1             BETA, DV2, INCY )
      IF ( LW ) WRITE ( IOUTF, * ) 
     1    'UBTCT1: Called DGBMV: dv2:= bm2 . sv2.'
C
C Temperature:
C   Multiply banded matrix BM3 by vector SV3 to give DV3
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = KL + KU + 1
      M      = NR*NH3
      N      = M
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, BM3, LDA, SV3, INCX,
     1             BETA, DV3, INCY )
      IF ( LW ) WRITE ( IOUTF, * ) 
     1    'UBTCT1: Called DGBMV: dv3:= bm3 . sv3.'
C
C Copy the diffusion forcing terms into VE1, VE2 and VE3 resp.
C
      N = NR*NH1
      CALL DCOPY( N, DV1, INCX, VE1, INCX )
C
      N = NR*NH2
      CALL DCOPY( N, DV2, INCX, VE2, INCX )
C
      N = NR*NH3
      CALL DCOPY( N, DV3, INCX, VE3, INCX )
C
C First, in order to add the heat source terms, we
C must "complete" the poloidal velocity vector.
C
      CALL PVVCPL( NR, NH1, NBN, SV1, PVLC, PVRC )
C
C We now need to form non-linear forcing terms for step i
C These are put into the vectors R1, R2 and R3.
C
      CALL UBNLTF( XARR, ML1, MM1, ML2, MM2, ML3, MM3, SV1, SV2,
     1    SV3, R1, R2, R3, SV3D, SV1Q, SV1S, SV2T, DSV3, VQ1, VS1,
     2    VT2, PA, DPA, FTF1, FTF2, FTF3, XSV, GAUX, GAUW, CVQ1, CVS1,
     3    CVT2Q,CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,CT2P,
     4    FTFPV,FTFTV,INFPV,INFTV,FTFPCV,FTFTCV,INFPCV,INFTCV,
     5    VGFA,IVGFA,SF2SA,ISF2S,PFAV,TFAV,JPFAV,JTFAV )
C
C We now add the non-linear terms (stored in R1, R2
C and R3) to VE1, VE2 and VE3
C
      N = NR*NH1
      CALL DAXPY( N, DELTAT, R1, INCX, VE1, INCX )
C
      N = NR*NH2
      CALL DAXPY( N, DELTAT, R2, INCX, VE2, INCX )
C
      N = NR*NH3
      CALL DAXPY( N, DELTAT, R3, INCX, VE3, INCX )
C
C Solve the system of equations for to form predictor
C We use the LAPACK routine DGBTRS
C First, solve for poloidal velocity:
C
      N      = NR*NH1
      LDA    = 3*NBN + 1
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM1, LDA, IPIV1, VE1, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBTCT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for predictor solution of VE1.'
      ENDIF
C
C Now, solve for toroidal velocity:
C
      N      = NR*NH2
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM2, LDA, IPIV2, VE2, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBTCT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for predictor solution of VE2.'
      ENDIF
C
C Now, solve for temperature
C
      N      = NR*NH3
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM3, LDA, IPIV3, VE3, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBTCT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for predictor solution of VE3.'
      ENDIF
C
C The predictor is now stored in VE1, VE2 and VE3
C Calculate the norm of this solution.
C
      N     = NR*NH1
      AVEC1 = DNRM2( N, VE1, INCX )
      N     = NR*NH2
      AVEC1 = AVEC1 + DNRM2( N, VE2, INCX )
      N     = NR*NH3
      AVEC1 = AVEC1 + DNRM2( N, VE3, INCX )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBTCT1: Predictor norm = ', AVEC1
C
C We then add the step^i forcing terms (non-linear)
C to the diffusion forcing terms. This will reduce
C the number of calculations done if more than one
C iteration is required.
C
      FAC = CFAC*DELTAT
C
      N   = NR*NH1
      CALL DAXPY( N, FAC, R1, INCX, DV1, INCX )
C
      N   = NR*NH2
      CALL DAXPY( N, FAC, R2, INCX, DV2, INCX )
C
      N   = NR*NH3
      CALL DAXPY( N, FAC, R3, INCX, DV3, INCX )
C
C Now begin the loop around corrector iterations
C
      NOIT  = 0
 50   CONTINUE
      NOIT = NOIT + 1
C
      IF ( NOIT.GT.ITMX ) THEN
        WRITE ( LLU, * ) 'Subroutine UBTCT1.'
        WRITE ( LLU, * ) 'Max iterations exceeded.'       
        NOIT = -1
        RETURN
      ENDIF
C
C Complete the poloidal velocity functions
C
      CALL PVVCPL( NR, NH1, NBN, VE1, PVLC, PVRC )
C
C Calculate the non-linear forcing terms in
C vectors R1, R2 and R3.
C
      CALL UBNLTF( XARR, ML1, MM1, ML2, MM2, ML3, MM3, VE1, VE2,
     1    VE3, R1, R2, R3, SV3D, SV1Q, SV1S, SV2T, DSV3, VQ1, VS1,
     2    VT2, PA, DPA, FTF1, FTF2, FTF3, XSV, GAUX, GAUW, CVQ1, CVS1,
     3    CVT2Q,CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,CT2P,
     4    FTFPV,FTFTV,INFPV,INFTV,FTFPCV,FTFTCV,INFPCV,INFTCV,
     5    VGFA,IVGFA,SF2SA,ISF2S,PFAV,TFAV,JPFAV,JTFAV )
C
      N = NR*NH1
      CALL DCOPY( N, DV1, INCX, VE1, INCX )
C
      N = NR*NH2
      CALL DCOPY( N, DV2, INCX, VE2, INCX )
C
      N = NR*NH3
      CALL DCOPY( N, DV3, INCX, VE3, INCX ) 
C
C Now add the step^{i+1} non-lin. terms to VE1, VE2 and VE3
C
      FAC = (1.0d0-CFAC)*DELTAT
C
      N   = NR*NH1
      CALL DAXPY( N, FAC, R1, INCX, VE1, INCX )
C
      N   = NR*NH2
      CALL DAXPY( N, FAC, R2, INCX, VE2, INCX )
C
      N   = NR*NH3
      CALL DAXPY( N, FAC, R3, INCX, VE3, INCX )
C
C Solve the system of equations for iteration NOIT
C We use the LAPACK routine DGBTRS
C First, solve for poloidal velocity:
C
      N      = NR*NH1
      LDA    = 3*NBN + 1
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM1, LDA, IPIV1, VE1, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBTCT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for solution ',NOIT,' of VE1.'
      ENDIF
C
C Now, solve for toroidal velocity:
C
      N      = NR*NH2
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM2, LDA, IPIV2, VE2, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBTCT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for solution ',NOIT,' of VE2.'
      ENDIF
C
C Now, solve for temperature
C
      N      = NR*NH3
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM3, LDA, IPIV3, VE3, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBTCT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for solution ',NOIT,' of VE3.'
      ENDIF
C
      AOLD  = AVEC1
      N     = NR*NH1
      AVEC1 = DNRM2( N, VE1, INCX )
      N     = NR*NH2
      AVEC1 = AVEC1 + DNRM2( N, VE2, INCX )
      N     = NR*NH3
      AVEC1 = AVEC1 + DNRM2( N, VE3, INCX )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBTCT1: Iteration ',NOIT,' norm = ', AVEC1
C
      DIFF = DABS( AVEC1 - AOLD )
      IF ( DIFF.LT.DTOL ) THEN
        IF ( LW ) WRITE ( IOUTF, * )
     1    'UBTCT1: Soln. converged iteration ', NOIT
          RETURN
      ENDIF
C
C Check to see if our norm appears to be getting bigger
C
      IF ( DIFF.GT.ODIFF .AND. NOIT.GT.3 ) THEN
        IF ( LW ) WRITE ( IOUTF, * )
     1    'UBTCT1: Soln. norm increasing: iteration ', NOIT
        NOIT = -2
        RETURN
      ENDIF
C
      ODIFF = DIFF
      GOTO 50
C
      END
C*********************************************************************

C*********************************************************************
C subroutine Uniform Boundary Non-Linear Terms Fill ******************
C            -       -        -   -      -     -    ******************
C Steve Gibbons Wed Nov 29 09:07:25 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Fills vectors R1, R2 and R3, the following terms in the heat and   C
C momentum equations:                                                C
C                                                                    C
C F_Theta =  u . ( CB1 r + CB2 r^{-2} , 0 , 0 ) :                    C
C            - CC u . Grad ( Theta )            :                    C
C                                                                    C
C F_vort. =     CH curl ( Theta )               :                    C
C             - CG curl ( K x v )               :                    C
C             - CF curl ( v. Grad) v            :                    C
C                                                                    C
C     SV1       : Dim ( NH1*NR ). Poloidal velocity initial soln.    C
C     SV2       : Dim ( NH2*NR ). Toroidal velocity initial soln.    C
C     SV3       : Dim ( NH3*NR ). Temperature initial soln.          C
C                                                                    C
C     R1        : Dim ( NH1*NR ). Toroidal vorticity forcing term    C
C     R2        : Dim ( NH2*NR ). Poloidal vorticity forcing term    C
C     R3        : Dim ( NH3*NR ). Temperature forcing term.          C
C                                                                    C
C     SV3D      : Dim ( 2*NBN + 1, NH3*NR ) takes 1st deriv. of SV3  C
C     SV1Q      : Dim (         1, NH1*NR ) turns p. vel. to Q( r )  C
C     SV1S      : Dim ( 2*NBN + 1, NH1*NR ) turns p. vel. to S( r )  C
C     SV2T      : Dim (         1, NH2*NR ) turns t. vel. to T( r )  C
C                                                                    C
C     CVQ1      : Dim (         1, NH1*NR ) takes curl of Q( r )     C
C     CVS1      : Dim ( 2*NBN + 1, NH1*NR ) takes curl of S( r )     C
C     CVT2Q     : Dim (         1, NH2*NR ) takes curl of T( r )     C
C     CVT2S     : Dim ( 2*NBN + 1, NH2*NR ) takes curl of T( r )     C
C                                                                    C
C     CQ1T      : Dim (         1, NH1*NR ) takes curl of Q( r )     C
C     CS1T      : Dim ( 2*NBN + 1, NH1*NR ) takes curl of S( r )     C
C     CT2P      : Dim (         1, NH2*NR ) takes curl of T( r )     C
C                                                                    C
C     DSV3      : Dim ( NH3*NR ). Work array.                        C
C     VQ1       : Dim ( NH1*NR ). Work array.                        C
C     VS1       : Dim ( NH1*NR ). Work array.                        C
C     VT2       : Dim ( NH2*NR ). Work array.                        C
C     VT1       : Dim ( NH1*NR ). Work array.                        C
C     VQ2       : Dim ( NH2*NR ). Work array.                        C
C     VS2       : Dim ( NH2*NR ). Work array.                        C
C                                                                    C
C     FTF1      : Dim (2*NPHP). Work array for fourier transforming. C
C     FTF2      : Dim (2*NPHP). Work array for fourier transforming. C
C     FTF3      : Dim (2*NPHP). Work array for fourier transforming. C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                  Both PA and DPA are formed by the routine SCHNLA. C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Corresponding Gauss quadrature weights. Dim (NTHP) C
C                                                                    C
C     XSV       : Work array. Dim ( NCMX, NPHP, NTHP, NR)            C
C                                                                    C
C     FTFPV     : Transform coeffs Dim ( 3, NH1, NTHP ) {from RSDV2C}C
C     FTFTV     : Transform coeffs Dim ( 2, NH2, NTHP ) {from RSDV2C}C
C                                                                    C
C     FTFPCV    : Transform coeffs Dim ( 3, NH2, NTHP ) {from RSDV2C}C
C     FTFTCV    : Transform coeffs Dim ( 2, NH1, NTHP ) {from RSDV2C}C
C                                                                    C
C     VGFA      : Dim ( 3, NH3, NTHP ). Coeff.s from SF2VGC          C
C     IVGFA     : Dim ( 2, NH3 ). Locationss from SF2VGC             C
C                                                                    C
C     SF2SA     : Dim ( NH3, NTHP ). Coeff.s from SF2SDC             C
C     ISF2S     : Dim ( NH3 ) Indices from SF2SDC                    C
C                                                                    C
C     PFAV      : Dim ( 3, NH1, NTHP ). Coeffs from XSVSDC.          C
C     TFAV      : Dim ( 2, NH2, NTHP ). Coeffs from XSVSDC.          C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE UBNLTF( XARR, ML1, MM1, ML2, MM2, ML3, MM3, SV1, SV2,
     1    SV3, R1, R2, R3, SV3D, SV1Q, SV1S, SV2T, DSV3, VQ1, VS1,
     2    VT2, PA, DPA, FTF1, FTF2, FTF3, XSV, GAUX, GAUW, CVQ1, CVS1,
     3    CVT2Q,CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,CT2P,
     4    FTFPV,FTFTV,INFPV,INFTV,FTFPCV,FTFTCV,INFPCV,INFTCV,
     5    VGFA,IVGFA,SF2SA,ISF2S,PFAV,TFAV,JPFAV,JTFAV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Common block contents .....................C
      INTEGER          NR, NH1, NH2, NH3, NBN, M0, IOUTF, LH,
     1                 ITMX, NOIT, NTHP, NPHP, NCMX
      COMMON  /NDIMPARS/  NR, NH1, NH2, NH3, NBN, M0, IOUTF, LH,
     1                    ITMX, NOIT, NTHP, NPHP, NCMX
      DOUBLE PRECISION CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CFAC, DELTAT, DTOL
      COMMON  /DPHYSPARS/ CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CFAC, DELTAT, DTOL
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION XARR( NR )
C
C Basic solution vector indices
C
      INTEGER          ML1( NH1 ), MM1( NH1 )
      INTEGER          ML2( NH2 ), MM2( NH2 )
      INTEGER          ML3( NH3 ), MM3( NH3 )
C
C Auxiliary matrices
C
      DOUBLE PRECISION SV3D( 2*NBN + 1, NH3*NR )
      DOUBLE PRECISION SV1Q(         1, NH1*NR )
      DOUBLE PRECISION SV1S( 2*NBN + 1, NH1*NR )
      DOUBLE PRECISION SV2T(         1, NH2*NR )
C
      DOUBLE PRECISION CVQ1(         1, NH1*NR )
      DOUBLE PRECISION CVS1(  2*NBN+1, NH1*NR )
      DOUBLE PRECISION CVT2Q(        1, NH2*NR )
      DOUBLE PRECISION CVT2S( 2*NBN+1, NH2*NR )
C
      DOUBLE PRECISION CQ1T(       1, NH1*NR )
      DOUBLE PRECISION CS1T( 2*NBN+1, NH1*NR )
      DOUBLE PRECISION CT2P(       1, NH2*NR )
C
C Initial solution vectors
C
      DOUBLE PRECISION SV1( NH1*NR ), SV2( NH2*NR ), SV3( NH3*NR )
C
C Forcing term vectors
C
      DOUBLE PRECISION R1( NH1*NR ), R2( NH2*NR ), R3( NH3*NR )
C
C Work arrays
C
      DOUBLE PRECISION DSV3( NH3*NR ), VQ1( NH1*NR ), VS1( NH1*NR ),
     1                 VT2( NH2*NR ), XSV( NCMX, NPHP, NTHP, NR),
     2                 VT1( NH1*NR ), VQ2( NH2*NR ), VS2( NH2*NR )
C
C Arrays necessary for transform
C
      DOUBLE PRECISION FTF1(2*NPHP), FTF2(2*NPHP), FTF3(2*NPHP),
     1                 GAUX( NTHP ), GAUW( NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP ),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP )
C
C Transform coefficients
C (pre-calculated by RSDV2C)
C
      DOUBLE PRECISION FTFPV( 3, NH1, NTHP ),
     1                 FTFTV( 2, NH2, NTHP )
C
      INTEGER          INFPV( 2, NH1 ),
     1                 INFTV( 2, NH2 )
C
      DOUBLE PRECISION FTFPCV( 3, NH2, NTHP ),
     1                 FTFTCV( 2, NH1, NTHP )
C
      INTEGER          INFPCV( 2, NH2 ),
     1                 INFTCV( 2, NH1 )
C
C Calculated by SF2VGC
C
      INTEGER          IVGFA( 2, NH3 )
      DOUBLE PRECISION VGFA( 3, NH3, NTHP )
C
C Calculated by SF2SDC
C
      INTEGER          ISF2S( NH3 )
      DOUBLE PRECISION SF2SA( NH3, NTHP )
C
C Calculated by XSVSDC
C
      INTEGER          JPFAV( 2, NH1 ), JTFAV( 2, NH2 )
      DOUBLE PRECISION PFAV( 3, NH1, NTHP ), TFAV( 2, NH2, NTHP )
C
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          K, N, IOP, INCX, INCY, LDA, M, ILNR, IRNR,
     1                 ICMR, ICMT, ICMP, ICM, IFORMF
      CHARACTER *(1)   TRANS
      LOGICAL          LW
      DOUBLE PRECISION FAC, ALPHA, BETA, ZERO
C
      PARAMETER ( TRANS = 'N', INCX   = 1, INCY   = 1, IOP = 0,
     1            ZERO = 0.0d0, IFORMF = 4 )
C____________________________________________________________________C
C
C Early return
C
      IF ( CB1.EQ.ZERO .AND. CB2.EQ.ZERO .AND. CH.EQ.ZERO .AND.
     1     CG.EQ.ZERO .AND. CC.EQ.ZERO .AND. CF.EQ.ZERO ) RETURN
C
      IF ( IOUTF.EQ.0 ) THEN
        LW  = .FALSE.
      ELSE
        LW  = .TRUE.
      ENDIF
C
C Zero R1, R2 and R3.
C
      FAC = 0.0d0
      N   = NR*NH1
      CALL VECOP( R1, FAC, N, IOP )
C
      N   = NR*NH2
      CALL VECOP( R2, FAC, N, IOP )
C
      N   = NR*NH3
      CALL VECOP( R3, FAC, N, IOP )
C
C R1, R2 and R3 are all now zero ...
C Add the heat-source terms to R3
C
      FAC = 1.0d0
      CALL NSVHST( NR, NH1, ML1, MM1, NH3, ML3, MM3, SV1, R3,
     1             FAC, XARR, CB1, CB2 )
C
C Now add buoyancy terms
C
      CALL NSVBTA( NR, NH3, ML3, MM3, NH1, ML1, MM1, SV3, R1, CH )
C
C New early escape?
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Heat and buoyancy terms added.'
      IF ( CG.EQ.ZERO .AND. CC.EQ.ZERO .AND. CF.EQ.ZERO ) RETURN
C
C Now we take the derivative of SV3 --> DSV3
C Use DGBMV to multiply SV3 by matrix SV3D to give DSV3
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1+2*NBN
      M      = NR*NH3
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, SV3D, LDA, SV3, INCX,
     1             BETA, DSV3, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: dsv3:= sv3d . sv3.'
C
C Now calculate scaloidal part of velocity.
C Use DGBMV to multiply SV1 by matrix SV1Q to give VQ1
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1
      M      = NR*NH1
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, SV1Q, LDA, SV1, INCX,
     1             BETA, VQ1, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: vq1:= sv1q . sv1.'
C
C Now calculate spheroidal part of velocity.
C Use DGBMV to multiply SV1 by matrix SV1S to give VS1
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1+2*NBN
      M      = NR*NH1
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, SV1S, LDA, SV1, INCX,
     1             BETA, VS1, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: vs1:= sv1s . sv1.'
C
C Now calculate toroidal part of velocity.
C Use DGBMV to multiply SV1 by matrix SV2T to give VT2
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1
      M      = NR*NH2
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, SV2T, LDA, SV2, INCX,
     1             BETA, VT2, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: vt2:= sv2t . sv2.'
C
C we now have the velocity in QST format and so can transform
C into real space.
C
      ILNR = 2
      IRNR = NR - 1
C     .  put radial component in icmr = 1
      ICMR = 1
C     .  put theta component in icmt = 2
      ICMT = 2
C     .  put phi component in icmp = 3
      ICMP = 3
      CALL RSDV2D( NTHP, NPHP,         NR, ILNR, IRNR, NH1,
     1                  NH2,           NCMX, ICMR, ICMT, ICMP,
     2             VQ1, VS1, VT2, XSV, FTF1, FTF2, FTF3,
     3             FTFPV, FTFTV, INFPV, INFTV, IFORMF )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called RSDV2D: XSV now contains velocity.'
C
C Velocity is now stored in XSV components 1-3.
C The temperature and temperature radial derivatives
C are now stored in SV3 and DSV3.
C We will loop around IR from ILNR to IRNR and calculate
C Grad( theta ) in the XSV array: components 4,5 and 6.
C
      ILNR = 2
      IRNR = NR - 1
C     .  put radial component in icmr = 4
      ICMR = 4
C     .  put theta component in icmt = 5
      ICMT = 5
C     .  put phi component in icmp = 6
      ICMP = 6
C
      CALL SF2VGD( ILNR, IRNR, NR, NH3, ML3, NTHP,
     1             NPHP, NCMX, ICMR, ICMT, ICMP, SV3, DSV3,
     2                   XARR,          XSV, FTF1, FTF2, FTF3,
     3             IVGFA, VGFA, IFORMF )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Grad( theta ) stored in XSV.'
C
C Now need to calculate curls of the scaloidal,
C spheroidal parts of the velocity which are stored
C in VQ1, VS1 and VT2 respectively ...
C The matrix CVQ1 multiplies VQ1 to give VT1 - which
C is a toroidal vector
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1
      M      = NR*NH1
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CVQ1, LDA, VQ1, INCX,
     1             BETA, VT1, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: vt1:= cvq1 . vq1.'
C
C Matrix CVS1 _adds_ the curl of VS1 to VT1
C
      ALPHA  = 1.0d0
      BETA   = 1.0d0
      LDA    = 1+2*NBN
      M      = NR*NH1
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CVS1, LDA, VS1, INCX,
     1             BETA, VT1, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: vt1:= cvq1 . vq1 + vt1.'
C
C The matrix CVT2Q multiplies VT2 to give VQ2 - which
C is a scaloidal vector
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1
      M      = NR*NH2
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CVT2Q, LDA, VT2, INCX,
     1             BETA, VQ2, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: vq2:= cvt2q . vt2.'
C
C The matrix CVT2S multiplies VT2 to give VS2 - which
C is a spheroidal vector
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1+2*NBN
      M      = NR*NH2
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CVT2S, LDA, VT2, INCX,
     1             BETA, VS2, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: vs2:= cvt2s . vt2.'
C
C The scaloidal vector of (curl v) is now in VQ2
C The spheroidal vector of (curl v) is now in VS2
C The toroidal vector of (curl v) is now in VT1
C So, calculate components in real space!
C
      ILNR = 2
      IRNR = NR - 1
C     .  put radial component in icmr = 7
      ICMR = 7
C     .  put theta component in icmt = 8
      ICMT = 8
C     .  put phi component in icmp = 9
      ICMP = 9
C
      CALL RSDV2D( NTHP, NPHP,         NR, ILNR, IRNR, NH2,
     1                  NH1,           NCMX, ICMR, ICMT, ICMP,
     2             VQ2, VS2, VT1,      XSV, FTF1, FTF2, FTF3,
     3             FTFPCV, FTFTCV, INFPCV, INFTCV, IFORMF )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called RSDV2D: XSV now contains curl vel.'
C
C Now evaluate all of the non-linear operations
C within XSV. This is now a single subroutine call!
C
      CALL NMCXSE( NTHP, NPHP, NR, ILNR, IRNR, NCMX, XSV, GAUX,
     1             CG, CF )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Calculated non-linear terms in real space.'
C
C v . Grad( theta ) is stored in element 10 of XSV.
C Add back to R3 in spectral coefficients
C
      ILNR = 2
      IRNR = NR - 1
      ICM  = 10
      CALL SF2SDD( ILNR, IRNR, NR,     NH3, ML3,          NTHP,
     1             NPHP, NCMX, ICM,           XSV, R3,      FTF1,
     2             SF2SA, IFORMF, ISF2S )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Added v.Grad( theta ) terms to R3.'
C
C Transform the non-linear terms back into QST space.
C Note that we will use VQ1, VS1 and VT2
C
      ILNR = 2
      IRNR = NR - 1
C     .  radial component is stored in icmr = 11
      ICMR = 11
C     .  theta component is stored in icmt = 12
      ICMT = 12
C     .  phi component is stored in icmp = 13
      ICMP = 13
      CALL XSVSDD( NTHP, NPHP,         NR, ILNR, IRNR, NH1,
     1                  NH2,           NCMX, ICMR, ICMT, ICMP,
     2             VQ1, VS1, VT2, XSV, FTF1, FTF2, FTF3,
     3             PFAV, TFAV, JPFAV, JTFAV, IFORMF )
C
C Now add the curl of the scaloidal function to the 
C toroidal forcing term of the vorticity (R1)
C
      ALPHA  = 1.0d0
      BETA   = 1.0d0
      LDA    = 1
      M      = NR*NH1
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CQ1T, LDA, VQ1, INCX,
     1             BETA, R1, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: r1:= r1 + cq1t . vq1.'
C
C Now add the curl of the spheroidal function to the
C toroidal forcing term of the vorticity (R1)
C
      ALPHA  = 1.0d0
      BETA   = 1.0d0
      LDA    = 1+2*NBN
      M      = NR*NH1
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CS1T, LDA, VS1, INCX,
     1             BETA, R1, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: r1:= r1 + cs1t . vs1.'
C
C Now add the curl of the toroidal function to the
C poloidal forcing term of the vorticity (R2)
C
      ALPHA  = 1.0d0
      BETA   = 1.0d0
      LDA    = 1
      M      = NR*NH2
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CT2P, LDA, VT2, INCX,
     1             BETA, R2, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBNLTF: Called DGBMV: r2:= r2 + ct2p . vt2.'
C
      RETURN
      END
C*********************************************************************
