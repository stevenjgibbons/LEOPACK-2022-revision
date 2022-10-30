C*********************************************************************
C                                                                    C
C                 Optimised                                          C
C Steve Gibbons - Uniform Boundary Convection + Dynamo Time Step     C
C Thu Nov  1 08:22:14 GMT 2001                                       C
C                                                                    C
C This is an experimental development program to try and derive      C
C a subroutine to advance solution by a single time-step.            C
C                                                                    C
C Version 2 also outputs energy and single component outputs ...     C
C                                                                    C
C Reads one extra line of input from the o2ibtctsc input file.       C
C This line contains three integer values ...                        C
C NTSBSE: the number of time-steps between standard energy evaluationC
C NTSBLE: ditto for L-kinetic energy spectrum information.           C
C NTSBME: ditto for M-kinetic energy spectrum information.           C
C                                                                    C
C Writes out an additional file .comps with vrad, vphi, temp         C
C (homogeneous and inhomogeneous parts) evaluated at                 C
C rad = (ri+ro)/2, the = pi/2 and phi = 0.0                          C
C                                                                    C
C                                                                    C
C*********************************************************************
      PROGRAM o2ubcdts2
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - contents of common blocks.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          NR, NH1, NH2, NH3, NH4, NH5, NBN, M0,
     1                 IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX
      COMMON  /NDIMPARS/  NR, NH1, NH2, NH3, NH4, NH5, NBN, M0,
     1                    IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX
      DOUBLE PRECISION CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CJ, CK, CL, CM, CFAC, DELTAT, DTOL
      COMMON  /DPHYSPARS/ CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CJ, CK, CL, CM, CFAC, DELTAT, DTOL
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          LHMAX, NTHMAX, NPHMAX, NPMAX
      PARAMETER      ( LHMAX = 52, NTHMAX = 100, NPHMAX = 100,
     1                 NPMAX = (LHMAX+1)*(LHMAX+2)/2 )
      INTEGER          NRMAX, NH1MAX, NH2MAX, NH3MAX, NH4MAX, NH5MAX,
     1                 NHMAX, NIVMAX, NDCS, NCFM, NDRVM, NBNM, NNDM,
     2                 NIV1MX, NIV2MX, NIV3MX, NIV4MX, NIV5MX, NCMXX
      PARAMETER      ( NRMAX = 150, NH1MAX = 1000, NH2MAX = 1000,
     1                 NH3MAX = 1000, NH4MAX = 1000, NH5MAX = 1000,
     2                 NHMAX = NH1MAX+NH2MAX+NH3MAX+NH4MAX+NH5MAX,
     3                 NCMXX = 22, NIVMAX = NHMAX*NRMAX )
      PARAMETER      ( NIV1MX = NH1MAX*NRMAX,
     1                 NIV2MX = NH2MAX*NRMAX,
     2                 NIV3MX = NH3MAX*NRMAX,
     3                 NIV4MX = NH4MAX*NRMAX,
     4                 NIV5MX = NH5MAX*NRMAX )
      PARAMETER      ( NDCS = 3+LHMAX, NBNM = 3, NCFM = 2*NBNM+1,
     1                 NDRVM = 4, NNDM = 4 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     1                 MHP( NHMAX ), LARR( NDCS ), MHIBC( NDCS ),
     2                 MHOBC( NDCS ), IWORK( NCFM ), IWRK( NNDM )
      INTEGER          ML1( NH1MAX ), MM1( NH1MAX ), MP1( NH1MAX )
      INTEGER          ML2( NH2MAX ), MM2( NH2MAX ), MP2( NH2MAX )
      INTEGER          ML3( NH3MAX ), MM3( NH3MAX ), MP3( NH3MAX )
      INTEGER          ML4( NH4MAX ), MM4( NH4MAX ), MP4( NH4MAX )
      INTEGER          ML5( NH5MAX ), MM5( NH5MAX ), MP5( NH5MAX )
      DOUBLE PRECISION SV( NIVMAX ), XARR( NRMAX ),
     1                 SVFDC( NCFM, NRMAX, NDRVM+1, NDCS ),
     2                 COEFM1( NCFM, NCFM ), W1( NCFM ),
     3                 COEFM2( NCFM, NCFM ), W2( NCFM ),
     4                 FDCM( NCFM, NRMAX, 1 )
      DOUBLE PRECISION WORK1( NNDM ),
     1                 WORK2( NNDM ), WORKM( NNDM, NNDM )
      DOUBLE PRECISION SV1( NIV1MX ), SV2( NIV2MX ), SV3( NIV3MX ),
     1                 SV4( NIV4MX ), SV5( NIV5MX )
      DOUBLE PRECISION E1( NIV1MX ), E2( NIV2MX ), E3( NIV3MX ),
     1                 E4( NIV4MX ), E5( NIV5MX )
      DOUBLE PRECISION D1( NIV1MX ), D2( NIV2MX ), D3( NIV3MX ),
     1                 D4( NIV4MX ), D5( NIV5MX )
      DOUBLE PRECISION DSV3( NIV3MX ), VQ1( NIV1MX ), VS1( NIV1MX ),
     1                 VT2( NIV2MX ), VT1( NIV1MX ), VQ2( NIV2MX ),
     2                 VS2( NIV2MX ), VQ4( NIV4MX ), VS4( NIV4MX ),
     3                 VT5( NIV5MX ), VT4( NIV4MX ), VQ5( NIV5MX ),
     4                 VS5( NIV5MX ), DKEARR( NHMAX )
      DOUBLE PRECISION R1( NIV1MX ), R2( NIV2MX ), R3( NIV3MX ),
     1                 R4( NIV4MX ), R5( NIV5MX )
C
C matrices for time-stepping
C
      DOUBLE PRECISION AM1( 3*NBNM+1, NIV1MX ),
     1                 BM1( 2*NBNM+1, NIV1MX )
      DOUBLE PRECISION AM2( 3*NBNM+1, NIV2MX ),
     1                 BM2( 2*NBNM+1, NIV2MX )
      DOUBLE PRECISION AM3( 3*NBNM+1, NIV3MX ),
     1                 BM3( 2*NBNM+1, NIV3MX )
      DOUBLE PRECISION AM4( 3*NBNM+1, NIV4MX ),
     1                 BM4( 2*NBNM+1, NIV4MX )
      DOUBLE PRECISION AM5( 3*NBNM+1, NIV5MX ),
     1                 BM5( 2*NBNM+1, NIV5MX )
      INTEGER          IP1( NIV1MX ), IP2( NIV2MX ),
     1                 IP3( NIV3MX ), IP4( NIV4MX ),
     2                 IP5( NIV5MX )
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
C matrices for forming QST decomposition of mag. field.
C
      DOUBLE PRECISION SV4Q(        1, NIV4MX )
      DOUBLE PRECISION SV4S( 2*NBNM+1, NIV4MX )
      DOUBLE PRECISION SV5T(        1, NIV5MX )
C
C matrices for taking curls of mag. QST decomposition
C
      DOUBLE PRECISION CVQ4(        1, NIV4MX )
      DOUBLE PRECISION CVS4( 2*NBNM+1, NIV4MX )
      DOUBLE PRECISION CVT5Q(        1, NIV5MX )
      DOUBLE PRECISION CVT5S( 2*NBNM+1, NIV5MX )
C
      DOUBLE PRECISION CQ5T(        1, NIV5MX )
      DOUBLE PRECISION CS5T( 2*NBNM+1, NIV5MX )
      DOUBLE PRECISION CT4P(        1, NIV4MX )
C
C Arrays for spherical transforms
C
      DOUBLE PRECISION GAUX( NTHMAX ), GAUW( NTHMAX ),
     1                 PA( NPMAX, NTHMAX ), DPA( NPMAX, NTHMAX )
      DOUBLE PRECISION XSV( NCMXX, NPHMAX, NTHMAX, NRMAX )
      DOUBLE PRECISION F1( 2*NPHMAX ), F2( 2*NPHMAX ),
     1                 F3( 2*NPHMAX )
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
C Magnetic field spectral --> real space
C
      INTEGER          INFPM( 2, NH4MAX ),
     1                 INFTM( 2, NH5MAX )
C
      DOUBLE PRECISION FTFPM( 3, NH4MAX, NTHMAX ),
     1                 FTFTM( 2, NH5MAX, NTHMAX )
C
C Curl velocity spectral --> real space
C
      INTEGER          INFPCV( 2, NH2MAX ),
     1                 INFTCV( 2, NH1MAX )
C
      DOUBLE PRECISION FTFPCV( 3, NH2MAX, NTHMAX ),
     1                 FTFTCV( 2, NH1MAX, NTHMAX )
C
C Curl magnetic field spectral --> real space
C
      INTEGER          INFPCM( 2, NH5MAX ),
     1                 INFTCM( 2, NH4MAX )
C
      DOUBLE PRECISION FTFPCM( 3, NH5MAX, NTHMAX ),
     1                 FTFTCM( 2, NH4MAX, NTHMAX )
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
C
      INTEGER          JPFAM( 2, NH5MAX ), JTFAM( 2, NH4MAX )
C
      DOUBLE PRECISION PFAM( 3, NH5MAX, NTHMAX ),
     1                 TFAM( 2, NH4MAX, NTHMAX )
C
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER           LU, LULOG, I, ILENRT, ILEN, IWRTE, IFORM,
     1                  NCUDS, INARR( 3 ), NH, NR1, NDRVS, ILN, IRN,
     2                  ICOMP, IDIR, MMAX, IOPT, NTSBSE, NTSBLE,
     3                  NTSBME, LUNRG, LUMNRG
      INTEGER           ITS, NTS, NTSBB, NTSBS, N, INC, NDIG, N1, N2,
     1                  K, L, M, LMINM, MMVAL, LLVAL, LUCOMP
      DOUBLE PRECISION  PVLC( NBNM ), PVRC( NBNM ), DMONE, DPONE,
     1                  STIME, DTIME, FAC, DKE( 2 ),
     2                  DTOTKE, DEATOT, DTORKE, DLOW, PI,
     3                  VRAD, VPHI, TEMP_TOT, RAD, THE, PHI,
     4                  DTOTME, DEATOM, DTORME, BTHE
      CHARACTER *(120)  LINE, FNLOG, ROOT, FNAME, FNNRG, FNCOMP,
     1                  FNMNRG
      CHARACTER *(10)   CHNUM
      DOUBLE PRECISION  RTPFCE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      EXTERNAL          RTPFCE
C
      PARAMETER ( IWRTE = 3, DMONE = -1.0d0, DPONE = 1.0d0,
     1            INC = 1, IFORM = 1, DLOW = 1.0d-9,
     2            PI=3.14159265358979312D0  )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NCMX    = NCMXX
      NBN     = NBNM
      LU      = 11
      LULOG   = 12
      LUNRG   = 13
      LUMNRG  = 14
      LUCOMP  = 15
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
C
      FNLOG(1:ILENRT) = ROOT(1:ILENRT)
      FNLOG(ILENRT+1:ILENRT+4) = '.log'
      FNLOG = FNLOG(1:ILENRT+4)
      CALL FOPEN( LULOG, FNLOG, IWRTE )
C
      FNNRG(1:ILENRT) = ROOT(1:ILENRT)
      FNNRG(ILENRT+1:ILENRT+4) = '.nrg'
      FNNRG = FNNRG(1:ILENRT+4)
      CALL FOPEN( LUNRG, FNNRG, IWRTE )
C
      FNCOMP(1:ILENRT) = ROOT(1:ILENRT)
      FNCOMP(ILENRT+1:ILENRT+6) = '.comps'
      FNCOMP = FNCOMP(1:ILENRT+6)
      CALL FOPEN( LUCOMP, FNCOMP, IWRTE )
C
      FNMNRG(1:ILENRT) = ROOT(1:ILENRT)
      FNMNRG(ILENRT+1:ILENRT+5) = '.mnrg'
      FNMNRG = FNMNRG(1:ILENRT+5)
      CALL FOPEN( LUMNRG, FNMNRG, IWRTE )
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
C Now put poloidal field harmonics into SV4
C
      ICOMP = 4
      CALL IIASCE( NH, ICOMP, MHT, MHL, MHM, MHP, NH4, NH4MAX,
     1             ML4, MM4, MP4 )
      WRITE ( LULOG, * ) 'ML4, MM4 and MP4 filled.'
      WRITE ( LULOG, * ) 'There are ',NH4,' poloidal field harms.'
      WRITE ( LULOG, * )
C
      IDIR  = 1
      ILN   = 2
      IRN   = NR - 1
      CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1             NH4, ML4, MM4, SV, SV4, IDIR )
C
      WRITE ( LULOG, * ) 'Initial SV4 vector filled.'
      WRITE ( LULOG, * )
C
C Now put toroidal field harmonics into SV5
C
      ICOMP = 5
      CALL IIASCE( NH, ICOMP, MHT, MHL, MHM, MHP, NH5, NH5MAX,
     1             ML5, MM5, MP5 )
      WRITE ( LULOG, * ) 'ML5, MM5 and MP5 filled.'
      WRITE ( LULOG, * ) 'There are ',NH5,' toroidal field harms.'
      WRITE ( LULOG, * )
C
      IDIR  = 1
      ILN   = 2
      IRN   = NR - 1
      CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1             NH5, ML5, MM5, SV, SV5, IDIR )
C
      WRITE ( LULOG, * ) 'Initial SV5 vector filled.'
      WRITE ( LULOG, * )
C
C Calculate M0, MMAX and LH
C (We ignore negative MM1, MM2, MM3, MM4 and MM5 as each -m will
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
      DO I = 1, NH4
        IF ( ML4( I ).GT.LH ) LH = ML4( I )
      ENDDO
      DO I = 1, NH5
        IF ( ML5( I ).GT.LH ) LH = ML5( I )
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
      DO I = 1, NH4
        IF ( MM4( I ).GT.0 .AND. MM4( I ).LT.M0 ) M0 = MM4( I )
        IF ( MM4( I ).GT.MMAX ) MMAX = MM4( I )
      ENDDO
      DO I = 1, NH5
        IF ( MM5( I ).GT.0 .AND. MM5( I ).LT.M0 ) M0 = MM5( I )
        IF ( MM5( I ).GT.MMAX ) MMAX = MM5( I )
      ENDDO
C
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
C (Non-magnetic terms only)
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
C Read in physical parameters: magnetic terms
C
      PRINT *,' Enter CJ, CK, CL, CM '
 215  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 215
      READ ( LINE, * ) CJ, CK, CL, CM
      WRITE ( LULOG, * ) 'CJ = ', CJ
      WRITE ( LULOG, * ) 'CK = ', CK
      WRITE ( LULOG, * ) 'CL = ', CL
      WRITE ( LULOG, * ) 'CM = ', CM
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
     1             NDRVM, ML1, MP1, IP1, XARR, SVFDC,
     2             AM1, BM1, CE, CI, CFAC, DELTAT )
      WRITE ( LULOG, * ) 'Called PVTSMF.'
C
      CALL TVTSMF( NR, NH2, NDCS, NBN, MHIBC, MHOBC, NCFM,
     1             NDRVM, ML2, MP2, IP2, XARR, SVFDC,
     2             AM2, BM2, CE, CI, CFAC, DELTAT )
      WRITE ( LULOG, * ) 'Called TVTSMF.'
C
      CALL TMTSMF( NR, NH3, NDCS, NBN, MHIBC, MHOBC, NCFM,
     1             NDRVM, ML3, MP3, IP3, XARR, SVFDC,
     2             AM3, BM3, CA, CD, CFAC, DELTAT )
      WRITE ( LULOG, * ) 'Called TMTSMF.'
C
      CALL PFTSMF( NR, NH4, NDCS, NBN, MHIBC, MHOBC, NCFM,
     1             NDRVM, ML4, MP4, IP4, XARR, SVFDC,
     2             AM4, BM4, CK, CL, CFAC, DELTAT )
      WRITE ( LULOG, * ) 'Called PFTSMF.'
C
      CALL TFTSMF( NR, NH5, NDCS, NBN, MHIBC, MHOBC, NCFM,
     1             NDRVM, ML5, MP5, IP5, XARR, SVFDC,
     2             AM5, BM5, CK, CL, CFAC, DELTAT )
      WRITE ( LULOG, * ) 'Called TFTSMF.'
C
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
C Now form arrays for dealing with the magnetic field.
C We will build the array SV4Q( 1, NH4*NR )
C which will multiply the solution vector SV4 to give
C a vector of scaloidal radial functions. Use
C VOBMBR as there is no derivative.
C
      N1     = 1
      N2     = NH4*NR
      IOPT   = 1
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH4, NBN, NCFM, ILN, IRN,
     1             ML4, IOPT, XARR, SV4Q, FDCM )
      WRITE ( LULOG, * ) 'Formed SV4Q.'
C
C Now we will build the array SV4S( 1+2*NBN, NH4*NR )
C which will multiply the solution vector SV4 to give
C a vector of spheroidal radial functions. Use
C AVBMBR for increased accuracy.
C
      N1     = 1+2*NBN
      N2     = NH4*NR
      IOPT   = 2
      ILN    = 2
      IRN    = NR-1
      K      = NBN
      CALL AVBMBR( N1, N2, K, NR, NH4, NBN, NCFM, NDRVM, ILN,
     1             IRN, ML4, MP4, IOPT, NDCS, XARR, SV4S, SVFDC )
      WRITE ( LULOG, * ) 'Formed SV4S.'
C
C Now we form SV5T( 1, NH5*NR )
C which will multiply the solution vector SV5 to give
C a vector of toroidal radial functions. (T as opposed
C to tau - see my PhD thesis). Use
C VOBMBR as there is no derivative.
C
      N1     = 1
      N2     = NH5*NR
      IOPT   = 3
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH5, NBN, NCFM, ILN, IRN,
     1             ML5, IOPT, XARR, SV5T, FDCM )
      WRITE ( LULOG, * ) 'Formed SV5T.'
C
C We now form CVQ4( 1, NH4*NR )
C which will multiply VQ4 to give VT4 -
C a vector of T(toroidal) radial functions.
C
      N1     = 1
      N2     = NH4*NR
      IOPT   = 4
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH4, NBN, NCFM, ILN, IRN,
     1             ML4, IOPT, XARR, CVQ4, FDCM )
      WRITE ( LULOG, * ) 'Formed CVQ4.'
C
C We now form CVS4( 1+2*NBN, NH4*NR )
C which will multiply VS4 to add to VT4 -
C a vector of T(toroidal) radial functions.
C
      N1     = 1+2*NBN
      N2     = NH4*NR
      IOPT   = 5
      K      = NBN
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH4, NBN, NCFM, ILN, IRN,
     1             ML4, IOPT, XARR, CVS4, FDCM )
      WRITE ( LULOG, * ) 'Formed CVS4.'
C
C We now form CVT5Q( 1, NH5*NR )
C which will multiply VT5 to form VQ5 -
C a vector of scaloidal radial functions.
C
      N1     = 1
      N2     = NH5*NR
      IOPT   = 6
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH5, NBN, NCFM, ILN, IRN,
     1             ML5, IOPT, XARR, CVT5Q, FDCM )
      WRITE ( LULOG, * ) 'Formed CVT5Q.'
C
C We now form CVT5S( 1+2*NBN, NH5*NR )
C which will multiply VT5 to form VS5 -
C a vector of spheroidal radial functions.
C
      N1     = 1+2*NBN
      N2     = NH5*NR
      IOPT   = 7
      K      = NBN
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH5, NBN, NCFM, ILN, IRN,
     1             ML5, IOPT, XARR, CVT5S, FDCM )
      WRITE ( LULOG, * ) 'Formed CVT5S.'
C
C We now form CQ5T( 1, NH5*NR )
C which will multiply VQ5 to form R5 -
C a vector of toroidal (tau) radial functions.
C
      N1     = 1
      N2     = NH5*NR
      IOPT   = 11
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH5, NBN, NCFM, ILN, IRN,
     1             ML5, IOPT, XARR, CQ5T, FDCM )
      WRITE ( LULOG, * ) 'Formed CQ5T.'
C
C We now form CS5T( 1+2*NBN, NH5*NR )
C which will multiply VS5 to form R5 -
C a vector of toroidal (tau) radial functions.
C
      N1     = 1+2*NBN
      N2     = NH5*NR
      IOPT   = 12
      K      = NBN
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH5, NBN, NCFM, ILN, IRN,
     1             ML5, IOPT, XARR, CS5T, FDCM )
      WRITE ( LULOG, * ) 'Formed CS5T.'
C
C We now form CT4P( 1, NH4*NR )
C which will multiply VT4 to form R4 -
C a vector of poloidal radial functions.
C
      N1     = 1
      N2     = NH4*NR
      IOPT   = 13
      K      = 0
      ILN    = 2
      IRN    = NR-1
      CALL VOBMBR( N1, N2, K, NR, NH4, NBN, NCFM, ILN, IRN,
     1             ML4, IOPT, XARR, CT4P, FDCM )
      WRITE ( LULOG, * ) 'Formed CT4P.'
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
C Now, the four arrays FTFPM( 3, NH4, NTHP ),
C FTFTM( 2, NH5, NTHP ), INFPM( 2, NH4 ) and INFTM( 2, NH5 )
C which transform the magnetic field into REAL space.
C
      CALL RSDV2C( NTHP, M0, LH, NH4, ML4,
     1             MM4, NH5, ML5, MM5, GAUX, PA, DPA,
     2             FTFPM, FTFTM, INFPM, INFTM )
      WRITE ( LULOG, * ) 'Called RSDV2C. (Mag. field)'
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
C Now, the four arrays FTFPCM( 3, NH5, NTHP ),
C FTFTCM( 2, NH4, NTHP ), INFPCM( 2, NH5 ) and INFTCM( 2, NH4 )
C which transform the magnetic field curl into REAL space.
C
      CALL RSDV2C( NTHP, M0, LH, NH5, ML5,
     1             MM5, NH4, ML4, MM4, GAUX, PA, DPA,
     2             FTFPCM, FTFTCM, INFPCM, INFTCM )
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
C Now the arrays PFAM( 3, NH5, NTHP ), TFAM( 2, NH4, NTHP ),
C JPFAM( 2, NH5 ) and JTFAM( 2, NH4 )
C
      CALL XSVSDC( NTHP, M0, LH, NH5, ML5, MM5, NH4, ML4, MM4,
     1             GAUX, GAUW, PA, DPA, PFAM, TFAM, JPFAM, JTFAM )
C
      WRITE ( LULOG, * ) 'Called XSVSDC. (induction eqn.)'
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
C Enter info. for energy evaluations
C
      PRINT *,' Please enter NTSBSE NTSBLE NTSBME'
 315  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 315
      READ ( LINE, * )  NTSBSE, NTSBLE, NTSBME
      WRITE ( LULOG, * ) 'NTSBSE = ', NTSBSE
      WRITE ( LULOG, * ) 'NTSBLE = ', NTSBLE
      WRITE ( LULOG, * ) 'NTSBME = ', NTSBME
      WRITE ( LULOG, * )
C
C Fix the values where flow components are evaluated
C
      RAD = 0.5d0*( XARR( 1 ) + XARR( NR ) )
      THE = 0.5d0*PI
      PHI = 0.0d0
C
C
C Begin the time-stepping procedure
C
      ITS    = 0
      DTIME  = STIME
 50   CONTINUE
      ITS    = ITS + 1
      DTIME  = DTIME + DELTAT
C
C Call main time-step routine
C Current solution is in SV1, SV2, SV3, SV4 and SV5
C
      CALL UBCDT1( XARR,ML1,MM1,IP1,ML2,MM2,IP2,ML3,MM3,IP3,ML4,MM4,
     1  IP4,ML5,MM5,IP5,AM1,BM1,AM2,BM2,AM3,BM3,AM4,BM4,AM5,BM5,SV1,
     2  SV2,SV3,SV4,SV5,D1,D2,D3,D4,D5,E1,E2,E3,E4,E5,R1,R2,R3,R4,R5,
     3  PVLC,PVRC,SV3D,SV1Q,SV1S,SV2T,DSV3,VQ1,VS1,VT2,PA,DPA,F1,F2,F3,
     4  XSV,GAUX,GAUW,CVQ1,CVS1,CVT2Q,CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,CT2P,
     5  SV4Q,SV4S,SV5T,CVQ4,CVS4,CVT5Q,CVT5S,VQ4,VS4,VT5,VT4,VQ5,VS5,
     6  CT4P,CQ5T,CS5T,FTFPV,FTFTV,INFPV,INFTV,
     7  FTFPM,FTFTM,INFPM,INFTM,FTFPCV,FTFTCV,INFPCV,INFTCV,
     8  FTFPCM,FTFTCM,INFPCM,INFTCM,VGFA,IVGFA,SF2SA,ISF2S,
     9  PFAV,TFAV,JPFAV,JTFAV,PFAM,TFAM,JPFAM,JTFAM)
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
C New solution is stored in E1, E2, E3, E4 and E5
C So copy them into SV1, SV2, SV3, SV4 and SV5
C
      N = NR*NH1
      CALL DCOPY( N, E1, INC, SV1, INC )
C
      N = NR*NH2
      CALL DCOPY( N, E2, INC, SV2, INC )
C
      N = NR*NH3
      CALL DCOPY( N, E3, INC, SV3, INC )
C
      N = NR*NH4
      CALL DCOPY( N, E4, INC, SV4, INC )
C
      N = NR*NH5
      CALL DCOPY( N, E5, INC, SV5, INC )
C
C See if we need to perform an energy evaluation
C
      IF ( ITS/NTSBSE*NTSBSE.EQ.ITS             .OR.
     1     ITS/NTSBLE*NTSBLE.EQ.ITS             .OR.
     2     ITS/NTSBME*NTSBME.EQ.ITS  ) THEN
C       .
C       . First need to gather up our solution vector
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
        ICOMP  = 4
        ILN    = 2
        IRN    = NR - 1
        CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1               NH4, ML4, MM4, SV, SV4, IDIR )
C       .
        ICOMP  = 5
        ILN    = 2
        IRN    = NR - 1
        CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1               NH5, ML5, MM5, SV, SV5, IDIR )
C       .
        DO I = 1, NH
          CALL SHKEER( I, NDCS, NR, INARR, MHT, MHL, MHP, NBN,
     1             NDRVS, NDRVM, NCFM, SV, XARR, DKE, SVFDC )
C         .
          DKEARR( I ) = DKE( 1 )
C         .
        ENDDO
C       .
C       . OK: The array DKEARR now contains the kinetic energy
C       . for each harmonic, so depending upon the information
C       . required, we now calculate the desired statistic.
C       .
C       . First is the general energy evaluation:
C       . Need to add to 3 variables DTOTKE, DEATOT, DTORKE
C       .
        IF ( ITS/NTSBSE*NTSBSE.EQ.ITS ) THEN
C         .
          DTOTKE = 0.0d0
          DEATOT = 0.0d0
          DTORKE = 0.0d0
          DTOTME = 0.0d0
          DEATOM = 0.0d0
          DTORME = 0.0d0
          DO I = 1, NH
C           .
C           .        find l and m
C           .
            L = MHL( I )
            IF ( MHM( I ).LT.0 ) THEN
              M = -MHM( I )
            ELSE
              M = MHM( I )
            ENDIF
            LMINM = L - M
C           .
C           . Case of poloidal velocity harmonic
C           .
            IF ( MHT( I ).EQ.1 ) THEN
C             .
C             . Add to the total in any case:
C             .
              DTOTKE = DTOTKE + DKEARR( I )
C             .
C             . Add to the EA total if (L-M) is odd.
C             .
              IF ( LMINM/2*2.NE.LMINM ) DEATOT = DEATOT + DKEARR( I )
C             .
            ENDIF
C           .
C           . Case of toroidal velocity harmonic
C           .
            IF ( MHT( I ).EQ.2 ) THEN
C             .
C             . Add to the total in any case:
C             .
              DTOTKE = DTOTKE + DKEARR( I )
C             .
C             . Add to the toroidal total in any case:
C             .
              DTORKE = DTORKE + DKEARR( I )
C             .
C             . Add to the EA total if (L-M) is even.
C             .
              IF ( LMINM/2*2.EQ.LMINM ) DEATOT = DEATOT + DKEARR( I )
C             .
            ENDIF
C           .
C           . Case of poloidal field harmonic
C           .
            IF ( MHT( I ).EQ.4 ) THEN
C             .
C             . Add to the total in any case:
C             .
              DTOTME = DTOTME + DKEARR( I )
C             .
C             . Add to the EA total if (L-M) is odd.
C             .
              IF ( LMINM/2*2.NE.LMINM ) DEATOM = DEATOM + DKEARR( I )
C             .
            ENDIF
C           .
C           . Case of toroidal field harmonic
C           .
            IF ( MHT( I ).EQ.5 ) THEN
C             .
C             . Add to the total in any case:
C             .
              DTOTME = DTOTME + DKEARR( I )
C             .
C             . Add to the toroidal total in any case:
C             .
              DTORME = DTORME + DKEARR( I )
C             .
C             . Add to the EA total if (L-M) is even.
C             .
              IF ( LMINM/2*2.EQ.LMINM ) DEATOM = DEATOM + DKEARR( I )
C             .
            ENDIF
C           .
          ENDDO
C         .
C         . Write out total values
C         .
          WRITE ( LUNRG, 737 ) DTIME, DTOTKE, DEATOT, DTORKE
          WRITE ( LUMNRG, 737 ) DTIME, DTOTME, DEATOM, DTORME
 737      FORMAT(1PD16.7,1PD16.7,1PD16.7,1PD16.7)
          CALL FLUSH( LUNRG )
          CALL FLUSH( LUMNRG )
C
C Evaluate the components of solution
C
C         .
C         . Evaluate radial velocity component
C         .
          ICOMP = 1
          VRAD  = RTPFCE( ICOMP, INARR, NNDM, IWRK, MHT, MHL, MHM, RAD,
     1                   THE, PHI, SV, XARR, WORK1, WORK2, WORKM )
C         .
C         . Evaluate phi velocity component
C         .
          ICOMP = 3
          VPHI  = RTPFCE( ICOMP, INARR, NNDM, IWRK, MHT, MHL, MHM, RAD,
     1                   THE, PHI, SV, XARR, WORK1, WORK2, WORKM )
C         .
C         . Evaluate theta magnetic field component
C         .
          ICOMP = 5
          BTHE  = RTPFCE( ICOMP, INARR, NNDM, IWRK, MHT, MHL, MHM, RAD,
     1                   THE, PHI, SV, XARR, WORK1, WORK2, WORKM )
C         .
C         . Evaluate temperature
C         .
          ICOMP     = 7
          TEMP_TOT  = RTPFCE( ICOMP, INARR, NNDM, IWRK, MHT, MHL, MHM,
     1                   RAD, THE, PHI, SV, XARR, WORK1, WORK2, WORKM )
C         .
C         . Write to LUCOMP, the values DTIME,
C         . VRAD, VPHI, BTHE, TEMP_TOT
C         .
          WRITE ( LUCOMP, 109 ) DTIME, VRAD, VPHI, BTHE, TEMP_TOT
 109      FORMAT (5(1PD16.7))
C
          CALL FLUSH( LUCOMP )
C         .
        ENDIF
C       .
C       . Now is the L-spectrum of kinetic energy
C       .
        IF ( ITS/NTSBLE*NTSBLE.EQ.ITS ) THEN
C         .
          DO LLVAL = 1, LH
C           .
            DTOTKE = 0.0d0
            DEATOT = 0.0d0
            DTORKE = 0.0d0
            DO I = 1, NH
C             .
C             .        find l and m
C             .
              L = MHL( I )
              IF ( L.NE.LLVAL ) GOTO 676
              IF ( MHM( I ).LT.0 ) THEN
                M = -MHM( I )
              ELSE
                M = MHM( I )
              ENDIF
              LMINM = L - M
C             .
C             .
C             . Case of poloidal harmonic
C             .
              IF ( MHT( I ).EQ.1 ) THEN
C               .
C               . Add to the total in any case:
C               .
                DTOTKE = DTOTKE + DKEARR( I )
C               .
C               . Add to the EA total if (L-M) is odd.
C               .
                IF ( LMINM/2*2.NE.LMINM ) DEATOT = DEATOT + DKEARR( I )
C               .
              ENDIF
C             .
C             . Case of toroidal harmonic
C             .
              IF ( MHT( I ).EQ.2 ) THEN
C               .
C               . Add to the total in any case:
C               .
                DTOTKE = DTOTKE + DKEARR( I )
C               .
C               . Add to the toroidal total in any case:
C               .
                DTORKE = DTORKE + DKEARR( I )
C               .
C               . Add to the EA total if (L-M) is even.
C               .
                IF ( LMINM/2*2.EQ.LMINM ) DEATOT = DEATOT + DKEARR( I )
C               .
              ENDIF
C             .
 676        CONTINUE
            ENDDO
C           .   end of I = 1, NH loop
C           . Write out spectrum for this L
C           .
            IF ( DTOTKE.GT.DLOW ) WRITE ( LULOG, 747 )
     1          LLVAL, DTIME, DTOTKE, DEATOT, DTORKE
 747      FORMAT('vL= ',I5,1PD16.7,1PD16.7,1PD16.7,1PD16.7)
          CALL FLUSH( LULOG )
C           .
          ENDDO
C         .   end of LLVAL = 1, LH loop
        ENDIF
C       . endif for L-spectrum
C       .
C       . Now is the M-spectrum of kinetic energy
C       .
        IF ( ITS/NTSBME*NTSBME.EQ.ITS ) THEN
C         .
          DO MMVAL = 0, MMAX, M0
C           .
            DTOTKE = 0.0d0
            DEATOT = 0.0d0
            DTORKE = 0.0d0
            DO I = 1, NH
C             .
C             .        find l and m
C             .
              IF ( MHM( I ).LT.0 ) THEN
                M = -MHM( I )
              ELSE
                M = MHM( I )
              ENDIF
              IF ( M.NE.MMVAL ) GOTO 677
              L = MHL( I )
              LMINM = L - M
C             .
C             . Case of poloidal harmonic
C             .
              IF ( MHT( I ).EQ.1 ) THEN
C               .
C               . Add to the total in any case:
C               .
                DTOTKE = DTOTKE + DKEARR( I )
C               .
C               . Add to the EA total if (L-M) is odd.
C               .
                IF ( LMINM/2*2.NE.LMINM ) DEATOT = DEATOT + DKEARR( I )
C               .
              ENDIF
C             .
C             . Case of toroidal harmonic
C             .
              IF ( MHT( I ).EQ.2 ) THEN
C               .
C               . Add to the total in any case:
C               .
                DTOTKE = DTOTKE + DKEARR( I )
C               .
C               . Add to the toroidal total in any case:
C               .
                DTORKE = DTORKE + DKEARR( I )
C               .
C               . Add to the EA total if (L-M) is even.
C               .
                IF ( LMINM/2*2.EQ.LMINM ) DEATOT = DEATOT + DKEARR( I )
C               .
              ENDIF
C             .
 677        CONTINUE
            ENDDO
C           .   end of I = 1, NH loop
C           . Write out spectrum for this L
C           .
            IF ( DTOTKE.GT.DLOW ) WRITE ( LULOG, 757 )
     1          MMVAL, DTIME, DTOTKE, DEATOT, DTORKE
 757      FORMAT('vM= ',I5,1PD16.7,1PD16.7,1PD16.7,1PD16.7)
          CALL FLUSH( LULOG )
C           .
          ENDDO
C         .   end of MMVAL = 0, MMAX loop
        ENDIF
C       . endif for M-spectrum
C       .
C       . Now deal with the magnetic energy spectrum
C       .
        IF ( ITS/NTSBLE*NTSBLE.EQ.ITS ) THEN
C         .
          DO LLVAL = 1, LH
C           .
            DTOTKE = 0.0d0
            DEATOT = 0.0d0
            DTORKE = 0.0d0
            DO I = 1, NH
C             .
C             .        find l and m
C             .
              L = MHL( I )
              IF ( L.NE.LLVAL ) GOTO 776
              IF ( MHM( I ).LT.0 ) THEN
                M = -MHM( I )
              ELSE
                M = MHM( I )
              ENDIF
              LMINM = L - M
C             .
C             .
C             . Case of poloidal harmonic
C             .
              IF ( MHT( I ).EQ.4 ) THEN
C               .
C               . Add to the total in any case:
C               .
                DTOTKE = DTOTKE + DKEARR( I )
C               .
C               . Add to the EA total if (L-M) is odd.
C               .
                IF ( LMINM/2*2.NE.LMINM ) DEATOT = DEATOT + DKEARR( I )
C               .
              ENDIF
C             .
C             . Case of toroidal harmonic
C             .
              IF ( MHT( I ).EQ.5 ) THEN
C               .
C               . Add to the total in any case:
C               .
                DTOTKE = DTOTKE + DKEARR( I )
C               .
C               . Add to the toroidal total in any case:
C               .
                DTORKE = DTORKE + DKEARR( I )
C               .
C               . Add to the EA total if (L-M) is even.
C               .
                IF ( LMINM/2*2.EQ.LMINM ) DEATOT = DEATOT + DKEARR( I )
C               .
              ENDIF
C             .
 776        CONTINUE
            ENDDO
C           .   end of I = 1, NH loop
C           . Write out spectrum for this L
C           .
            IF ( DTOTKE.GT.DLOW ) WRITE ( LULOG, 767 )
     1          LLVAL, DTIME, DTOTKE, DEATOT, DTORKE
 767      FORMAT('mL= ',I5,1PD16.7,1PD16.7,1PD16.7,1PD16.7)
          CALL FLUSH( LULOG )
C           .
          ENDDO
C         .   end of LLVAL = 1, LH loop
        ENDIF
C       . endif for L-spectrum
C       .
C       . Now is the M-spectrum of kinetic energy
C       .
        IF ( ITS/NTSBME*NTSBME.EQ.ITS ) THEN
C         .
          DO MMVAL = 0, MMAX, M0
C           .
            DTOTKE = 0.0d0
            DEATOT = 0.0d0
            DTORKE = 0.0d0
            DO I = 1, NH
C             .
C             .        find l and m
C             .
              IF ( MHM( I ).LT.0 ) THEN
                M = -MHM( I )
              ELSE
                M = MHM( I )
              ENDIF
              IF ( M.NE.MMVAL ) GOTO 877
              L = MHL( I )
              LMINM = L - M
C             .
C             . Case of poloidal harmonic
C             .
              IF ( MHT( I ).EQ.4 ) THEN
C               .
C               . Add to the total in any case:
C               .
                DTOTKE = DTOTKE + DKEARR( I )
C               .
C               . Add to the EA total if (L-M) is odd.
C               .
                IF ( LMINM/2*2.NE.LMINM ) DEATOT = DEATOT + DKEARR( I )
C               .
              ENDIF
C             .
C             . Case of toroidal harmonic
C             .
              IF ( MHT( I ).EQ.5 ) THEN
C               .
C               . Add to the total in any case:
C               .
                DTOTKE = DTOTKE + DKEARR( I )
C               .
C               . Add to the toroidal total in any case:
C               .
                DTORKE = DTORKE + DKEARR( I )
C               .
C               . Add to the EA total if (L-M) is even.
C               .
                IF ( LMINM/2*2.EQ.LMINM ) DEATOT = DEATOT + DKEARR( I )
C               .
              ENDIF
C             .
 877        CONTINUE
            ENDDO
C           .   end of I = 1, NH loop
C           . Write out spectrum for this L
C           .
            IF ( DTOTKE.GT.DLOW ) WRITE ( LULOG, 777 )
     1          MMVAL, DTIME, DTOTKE, DEATOT, DTORKE
 777      FORMAT('mM= ',I5,1PD16.7,1PD16.7,1PD16.7,1PD16.7)
          CALL FLUSH( LULOG )
C           .
          ENDDO
C         .   end of MMVAL = 0, MMAX loop
        ENDIF
C       . endif for M-spectrum
      ENDIF
C     . endif for outputting of energy spectra
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
        ICOMP  = 4
        ILN    = 2
        IRN    = NR - 1
        CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1               NH4, ML4, MM4, SV, SV4, IDIR )
C       .
        ICOMP  = 5
        ILN    = 2
        IRN    = NR - 1
        CALL MC2SCV( INARR, MHT, MHL, MHM, ILN, IRN, ICOMP,
     1               NH5, ML5, MM5, SV, SV5, IDIR )
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
C Return to perform next time-step
C
      IF ( ITS.EQ.NTS ) GOTO 999
      GOTO 50
C
 999  CONTINUE
      CALL FCLOSE( LULOG, FNLOG, 'Error' )
      CALL FCLOSE( LUCOMP, FNCOMP, 'Error' )
      CALL FCLOSE( LUNRG, FNNRG, 'Error' )
      CALL FCLOSE( LUMNRG, FNMNRG, 'Error' )
      STOP
      END
C*********************************************************************

C*********************************************************************
C Subroutine Uniform Boundary Convection + Dynamo Time step 1 ********
C            -       -        -            -      -         - ********
C Steve Gibbons Mon Dec 11 13:06:21 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C We are trying to time-step the heat, momentum and induction eqn.s  C
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
C                       - CJ [ B x (curl B) ]                        C
C                                                                    C
C  where RADVEC is the vector ( R, 0, 0 ) and k is the vector        C
C  ( cos theta, -sin theta, 0 ).                                     C
C                                                                    C
C  The induction equation is                                         C
C                                                                    C
C    CK d B / dt      =   CL Lap B                                   C
C                       + CM curl ( v x B )                          C
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
C  Contains the integer numbers NR, NH1, NH2, NH3, NH4, NH5, NBN,    C
C                    M0, IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX     C
C                                                                    C
C  where                                                             C
C                                                                    C
C     NR        : Total number of radial grid nodes.                 C
C     NH1       : Number of poloidal velocity harmonics in soln.     C
C     NH2       : Number of toroidal velocity harmonics in soln.     C
C     NH3       : Number of temperature harmonics in soln.           C
C     NH4       : Number of poloidal field harmonics in soln.        C
C     NH5       : Number of toroidal field harmonics in soln.        C
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
C     CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,                      C
C      CJ, CK, CL, CM : see above description. of equations.         C
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
C     IP1       : Dim (NH1*NR). Pivotting information for AM1        C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C     ML2       : Array dim ( NH2 ). Sph. harm degree, L.            C
C     MM2       : Array dim ( NH2 ). Sph. harm order, M, or -M.      C
C            ( ml2 and mm2 describe toroidal velocity )              C
C     IP2       : Dim (NH2*NR). Pivotting information for AM2        C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C     ML3       : Array dim ( NH3 ). Sph. harm degree, L.            C
C     MM3       : Array dim ( NH3 ). Sph. harm order, M, or -M.      C
C            ( ml3 and mm3 describe temperature )                    C
C     IP3       : Dim (NH3*NR). Pivotting information for AM3        C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C     ML4       : Array dim ( NH4 ). Sph. harm degree, L.            C
C     MM4       : Array dim ( NH4 ). Sph. harm order, M, or -M.      C
C            ( ml4 and mm4 describe poloidal magnetic field)         C
C     IP4       : Dim (NH4*NR). Pivotting information for AM4        C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C     ML5       : Array dim ( NH5 ). Sph. harm degree, L.            C
C     MM5       : Array dim ( NH5 ). Sph. harm order, M, or -M.      C
C            ( ml5 and mm5 describe toroidal magnetic field)         C
C     IP5       : Dim (NH5*NR). Pivotting information for AM5        C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C     INFPV     : Dim ( 2, NH1 ) {from RSDV2C}                       C
C     INFTV     : Dim ( 2, NH2 ) {from RSDV2C}                       C
C                                                                    C
C     INFPM     : Dim ( 2, NH4 ) {from RSDV2C}                       C
C     INFTM     : Dim ( 2, NH5 ) {from RSDV2C}                       C
C                                                                    C
C     INFPCV    : Dim ( 2, NH2 ) {from RSDV2C}                       C
C     INFTCV    : Dim ( 2, NH1 ) {from RSDV2C}                       C
C                                                                    C
C     INFPCM    : Dim ( 2, NH5 ) {from RSDV2C}                       C
C     INFTCM    : Dim ( 2, NH4 ) {from RSDV2C}                       C
C                                                                    C
C     JPFAV     : Dim(2,NH1). Locations in FTF2/3 calc. by XSVSDC    C
C     JTFAV     : Dim(2,NH2). Locations in FTF2/3 calc. by XSVSDC    C
C                                                                    C
C     JPFAM     : Dim(2,NH5). Locations in FTF2/3 calc. by XSVSDC    C
C     JTFAM     : Dim(2,NH4). Locations in FTF2/3 calc. by XSVSDC    C
C                                                                    C
C     IVGFA     : Dim ( 2, NH3 ). Locations from SF2VGC              C
C                                                                    C
C     ISF2S     : Dim ( NH3 ). Indices from SF2SDC                   C
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
C     AM4       : Matrix for solution of f^{i+1}. (poloidal m.f.)    C
C                 Has dimensions ( 3*NBN + 1, NH4*NR )               C
C                                                                    C
C     BM4       : Matrix for multiplication of f^i. (poloidal m.f.)  C
C                 Has dimensions ( 2*NBN + 1, NH4*NR )               C
C                                                                    C
C     AM5       : Matrix for solution of f^{i+1}. (toroidal m.f.)    C
C                 Has dimensions ( 3*NBN + 1, NH5*NR )               C
C                                                                    C
C     BM5       : Matrix for multiplication of f^i. (toroidal m.f.)  C
C                 Has dimensions ( 2*NBN + 1, NH5*NR )               C
C                                                                    C
C     SV1       : Dim ( NH1*NR ). Poloidal velocity initial soln.    C
C     SV2       : Dim ( NH2*NR ). Toroidal velocity initial soln.    C
C     SV3       : Dim ( NH3*NR ). Temperature initial soln.          C
C     SV4       : Dim ( NH4*NR ). Poloidal field initial soln.       C
C     SV5       : Dim ( NH5*NR ). Toroidal field initial soln.       C
C                                                                    C
C     E1        : Dim ( NH1*NR ). Poloidal velocity final soln.      C
C     E2        : Dim ( NH2*NR ). Toroidal velocity final soln.      C
C     E3        : Dim ( NH3*NR ). Temperature final soln.            C
C     E4        : Dim ( NH4*NR ). Poloidal field final soln.         C
C     E5        : Dim ( NH5*NR ). Toroidal field final soln.         C
C                                                                    C
C     D1        : Dim ( NH1*NR ). Work array.                        C
C     D2        : Dim ( NH2*NR ). Work array.                        C
C     D3        : Dim ( NH3*NR ). Work array.                        C
C     D4        : Dim ( NH4*NR ). Work array.                        C
C     D5        : Dim ( NH5*NR ). Work array.                        C
C                                                                    C
C     R1        : Dim ( NH1*NR ). Work array.                        C
C     R2        : Dim ( NH2*NR ). Work array.                        C
C     R3        : Dim ( NH3*NR ). Work array.                        C
C     R4        : Dim ( NH4*NR ). Work array.                        C
C     R5        : Dim ( NH5*NR ). Work array.                        C
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
C     SV4Q      : Dim (         1, NH4*NR ) turns p. mag. to Q( r )  C
C     SV4S      : Dim ( 2*NBN + 1, NH4*NR ) turns p. mag. to S( r )  C
C     SV5T      : Dim (         1, NH5*NR ) turns t. mag. to T( r )  C
C                                                                    C
C     CVQ4      : Dim (         1, NH4*NR ) takes curl of Q( r )     C
C     CVS4      : Dim ( 2*NBN + 1, NH4*NR ) takes curl of S( r )     C
C     CVT5Q     : Dim (         1, NH5*NR ) takes curl of T( r )     C
C     CVT5S     : Dim ( 2*NBN + 1, NH5*NR ) takes curl of T( r )     C
C                                                                    C
C     CQ5T      : Dim (         1, NH5*NR ) takes curl of Q( r )     C
C     CS5T      : Dim ( 2*NBN + 1, NH5*NR ) takes curl of S( r )     C
C     CT4P      : Dim (         1, NH4*NR ) takes curl of T( r )     C
C                                                                    C
C     F1        : Dim (2*NPHP). Work array for fourier transforming. C
C     F2        : Dim (2*NPHP). Work array for fourier transforming. C
C     F3        : Dim (2*NPHP). Work array for fourier transforming. C
C                                                                    C
C     DSV3      : Dim ( NH3*NR ). Work array.                        C
C     VQ1       : Dim ( NH1*NR ). Work array.                        C
C     VS1       : Dim ( NH1*NR ). Work array.                        C
C     VT2       : Dim ( NH2*NR ). Work array.                        C
C     VT1       : Dim ( NH1*NR ). Work array.                        C
C     VQ2       : Dim ( NH2*NR ). Work array.                        C
C     VS2       : Dim ( NH2*NR ). Work array.                        C
C     VQ4       : Dim ( NH4*NR ). Work array.                        C
C     VS4       : Dim ( NH4*NR ). Work array.                        C
C     VT5       : Dim ( NH5*NR ). Work array.                        C
C     VT4       : Dim ( NH4*NR ). Work array.                        C
C     VQ5       : Dim ( NH5*NR ). Work array.                        C
C     VS5       : Dim ( NH5*NR ). Work array.                        C
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
C     FTFPM     : Transform coeffs Dim ( 3, NH4, NTHP ) {from RSDV2C}C
C     FTFTM     : Transform coeffs Dim ( 2, NH5, NTHP ) {from RSDV2C}C
C                                                                    C
C     FTFPCV    : Transform coeffs Dim ( 3, NH2, NTHP ) {from RSDV2C}C
C     FTFTCV    : Transform coeffs Dim ( 2, NH1, NTHP ) {from RSDV2C}C
C                                                                    C
C     FTFPCM    : Transform coeffs Dim ( 3, NH5, NTHP ) {from RSDV2C}C
C     FTFTCM    : Transform coeffs Dim ( 2, NH4, NTHP ) {from RSDV2C}C
C                                                                    C
C     VGFA      : Dim ( 3, NH3, NTHP ). Coeff.s from SF2VGC          C
C                                                                    C
C     SF2SA     : Dim ( NH3, NTHP ). Coeff.s from SF2SDC             C
C                                                                    C
C     PFAV      : Dim ( 3, NH1, NTHP ). Coeffs from XSVSDC.          C
C     TFAV      : Dim ( 2, NH2, NTHP ). Coeffs from XSVSDC.          C
C                                                                    C
C     PFAM      : Dim ( 3, NH5, NTHP ). Coeffs from XSVSDC.          C
C     TFAM      : Dim ( 2, NH4, NTHP ). Coeffs from XSVSDC.          C
C                                                                    C
C*********************************************************************
      SUBROUTINE UBCDT1( XARR,ML1,MM1,IP1,ML2,MM2,IP2,ML3,MM3,IP3,ML4,
     1 MM4,IP4,ML5,MM5,IP5,AM1,BM1,AM2,BM2,AM3,BM3,AM4,BM4,AM5,BM5,
     2 SV1,SV2,SV3,SV4,SV5,D1,D2,D3,D4,D5,E1,E2,E3,E4,E5,R1,R2,R3,R4,
     3 R5,PVLC,PVRC,SV3D,SV1Q,SV1S,SV2T,DSV3,VQ1,VS1,VT2,PA,DPA,F1,F2,
     4 F3,XSV,GAUX,GAUW,CVQ1,CVS1,CVT2Q,CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,
     5 CT2P,SV4Q,SV4S,SV5T,CVQ4,CVS4,CVT5Q,CVT5S,VQ4,VS4,VT5,VT4,VQ5,
     6 VS5,CT4P,CQ5T,CS5T,FTFPV,FTFTV,INFPV,INFTV,
     7 FTFPM,FTFTM,INFPM,INFTM,FTFPCV,FTFTCV,INFPCV,INFTCV,
     8 FTFPCM,FTFTCM,INFPCM,INFTCM,VGFA,IVGFA,SF2SA,ISF2S,
     9 PFAV,TFAV,JPFAV,JTFAV,PFAM,TFAM,JPFAM,JTFAM)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Common block contents .....................C
      INTEGER          NR, NH1, NH2, NH3, NH4, NH5, NBN, M0,
     1                 IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX
      COMMON  /NDIMPARS/  NR, NH1, NH2, NH3, NH4, NH5, NBN, M0,
     1                    IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX
      DOUBLE PRECISION CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CJ, CK, CL, CM, CFAC, DELTAT, DTOL
      COMMON  /DPHYSPARS/ CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CJ, CK, CL, CM, CFAC, DELTAT, DTOL
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION XARR( NR )
C
C Basic solution vector indices
C
      INTEGER          ML1( NH1 ), MM1( NH1 ), IP1( NH1*NR )
      INTEGER          ML2( NH2 ), MM2( NH2 ), IP2( NH2*NR )
      INTEGER          ML3( NH3 ), MM3( NH3 ), IP3( NH3*NR )
      INTEGER          ML4( NH4 ), MM4( NH4 ), IP4( NH4*NR )
      INTEGER          ML5( NH5 ), MM5( NH5 ), IP5( NH5*NR )
C
C Diffusion/time-step matrices
C
      DOUBLE PRECISION AM1( 3*NBN+1, NH1*NR ), BM1( 2*NBN+1, NH1*NR )
      DOUBLE PRECISION AM2( 3*NBN+1, NH2*NR ), BM2( 2*NBN+1, NH2*NR )
      DOUBLE PRECISION AM3( 3*NBN+1, NH3*NR ), BM3( 2*NBN+1, NH3*NR )
      DOUBLE PRECISION AM4( 3*NBN+1, NH4*NR ), BM4( 2*NBN+1, NH4*NR )
      DOUBLE PRECISION AM5( 3*NBN+1, NH5*NR ), BM5( 2*NBN+1, NH5*NR )
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
      DOUBLE PRECISION SV4Q(         1, NH4*NR )
      DOUBLE PRECISION SV4S( 2*NBN + 1, NH4*NR )
      DOUBLE PRECISION SV5T(         1, NH5*NR )
C
      DOUBLE PRECISION CVQ4(          1, NH4*NR )
      DOUBLE PRECISION CVS4(  2*NBN + 1, NH4*NR )
      DOUBLE PRECISION CVT5Q(         1, NH5*NR )
      DOUBLE PRECISION CVT5S( 2*NBN + 1, NH5*NR )
C
      DOUBLE PRECISION CQ5T(       1, NH5*NR )
      DOUBLE PRECISION CS5T( 2*NBN+1, NH5*NR )
      DOUBLE PRECISION CT4P(       1, NH4*NR )
C
C Initial solution vectors
C
      DOUBLE PRECISION SV1( NH1*NR ), SV2( NH2*NR ), SV3( NH3*NR ),
     1                 SV4( NH4*NR ), SV5( NH5*NR )
C
C Final solution vectors
C
      DOUBLE PRECISION E1( NH1*NR ), E2( NH2*NR ), E3( NH3*NR ),
     1                 E4( NH4*NR ), E5( NH5*NR )
C
C Arrays necessary for transform
C
      DOUBLE PRECISION F1(2*NPHP), F2(2*NPHP), F3(2*NPHP),
     1                 GAUX( NTHP ), GAUW( NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP ),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP )
C
C Work arrays
C
      DOUBLE PRECISION D1( NH1*NR ), D2( NH2*NR ), D3( NH3*NR ),
     1                 D4( NH4*NR ), D5( NH5*NR )
      DOUBLE PRECISION R1( NH1*NR ),  R2( NH2*NR ),  R3( NH3*NR ),
     1                 R4( NH4*NR ),  R5( NH5*NR )
      DOUBLE PRECISION DSV3( NH3*NR ), VQ1( NH1*NR ), VS1( NH1*NR ),
     1                 VT2( NH2*NR ), XSV( NCMX, NPHP, NTHP, NR),
     2                 VT1( NH1*NR ), VQ2( NH2*NR ), VS2( NH2*NR )
      DOUBLE PRECISION VQ4( NH4*NR ), VS4( NH4*NR ), VT5( NH5*NR ),
     1                 VT4( NH4*NR ), VQ5( NH5*NR ), VS5( NH5*NR )
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
      DOUBLE PRECISION FTFPM( 3, NH4, NTHP ),
     1                 FTFTM( 2, NH5, NTHP )
C
      INTEGER          INFPM( 2, NH4 ),
     1                 INFTM( 2, NH5 )
C
      DOUBLE PRECISION FTFPCV( 3, NH2, NTHP ),
     1                 FTFTCV( 2, NH1, NTHP )
C
      INTEGER          INFPCV( 2, NH2 ),
     1                 INFTCV( 2, NH1 )
C
      DOUBLE PRECISION FTFPCM( 3, NH5, NTHP ),
     1                 FTFTCM( 2, NH4, NTHP )
C
      INTEGER          INFPCM( 2, NH5 ),
     1                 INFTCM( 2, NH4 )
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
C Calculated by XSVSDA
C
      INTEGER          JPFAV( 2, NH1 ), JTFAV( 2, NH2 )
      DOUBLE PRECISION PFAV( 3, NH1, NTHP ), TFAV( 2, NH2, NTHP )
C
      INTEGER          JPFAM( 2, NH5 ), JTFAM( 2, NH4 )
      DOUBLE PRECISION PFAM( 3, NH5, NTHP ), TFAM( 2, NH4, NTHP )
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
      IF ( LW ) WRITE ( IOUTF, * ) 'Entered UBCDT1.'
C
C Calculate the contribution from the diffusive parts.
C Poloidal velocity:
C   Multiply banded matrix BM1 by vector SV1 to give D1
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = KL + KU + 1
      M      = NR*NH1
      N      = M
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, BM1, LDA, SV1, INCX,
     1             BETA, D1, INCY )
      IF ( LW ) WRITE ( IOUTF, * ) 
     1    'UBCDT1: Called DGBMV: dv1:= bm1 . sv1.'
C
C Toroidal velocity:
C   Multiply banded matrix BM2 by vector SV2 to give D2
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = KL + KU + 1
      M      = NR*NH2
      N      = M
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, BM2, LDA, SV2, INCX,
     1             BETA, D2, INCY )
      IF ( LW ) WRITE ( IOUTF, * ) 
     1    'UBCDT1: Called DGBMV: dv2:= bm2 . sv2.'
C
C Temperature:
C   Multiply banded matrix BM3 by vector SV3 to give D3
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = KL + KU + 1
      M      = NR*NH3
      N      = M
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, BM3, LDA, SV3, INCX,
     1             BETA, D3, INCY )
      IF ( LW ) WRITE ( IOUTF, * ) 
     1    'UBCDT1: Called DGBMV: dv3:= bm3 . sv3.'
C
C Poloidal magnetic field:
C   Multiply banded matrix BM4 by vector SV4 to give D4
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = KL + KU + 1
      M      = NR*NH4
      N      = M
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, BM4, LDA, SV4, INCX,
     1             BETA, D4, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBCDT1: Called DGBMV: dv4:= bm4 . sv4.'
C
C Toroidal magnetic field:
C   Multiply banded matrix BM5 by vector SV5 to give D5
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = KL + KU + 1
      M      = NR*NH5
      N      = M
      CALL DGBMV ( TRANS, M, N, KL, KU, ALPHA, BM5, LDA, SV5, INCX,
     1             BETA, D5, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBCDT1: Called DGBMV: dv5:= bm5 . sv5.'
C
C Copy the diffusion forcing terms into E1, E2, E3, E4 and E5.
C
      N = NR*NH1
      CALL DCOPY( N, D1, INCX, E1, INCX )
C
      N = NR*NH2
      CALL DCOPY( N, D2, INCX, E2, INCX )
C
      N = NR*NH3
      CALL DCOPY( N, D3, INCX, E3, INCX )
C
      N = NR*NH4
      CALL DCOPY( N, D4, INCX, E4, INCX )
C
      N = NR*NH5
      CALL DCOPY( N, D5, INCX, E5, INCX )
C
C First, in order to add the heat source terms, we
C must "complete" the poloidal velocity vector.
C
      CALL PVVCPL( NR, NH1, NBN, SV1, PVLC, PVRC )
C
C We now need to form non-linear forcing terms for step i
C These are put into the vectors R1, R2, R3, R4 and R5.
C
      CALL OUBNLTF( XARR,ML1,MM1,ML2,MM2,ML3,MM3,ML4,MM4,ML5,MM5,
     1    SV1,SV2,SV3,SV4,SV5,R1,R2,R3,R4,R5,SV3D,SV1Q,SV1S,SV2T,DSV3,
     2    VQ1,VS1,VT2,PA,DPA,F1,F2,F3,XSV,GAUX,GAUW,CVQ1,CVS1,CVT2Q,
     3    CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,CT2P,SV4Q,SV4S,SV5T,CVQ4,CVS4,
     4    CVT5Q,CVT5S,VQ4,VS4,VT5,VT4,VQ5,VS5,CT4P,CQ5T,CS5T,
     5    FTFPV,FTFTV,INFPV,INFTV,FTFPM,FTFTM,INFPM,INFTM,
     6    FTFPCV,FTFTCV,INFPCV,INFTCV,FTFPCM,FTFTCM,INFPCM,INFTCM,
     7    VGFA,IVGFA,SF2SA,ISF2S,
     8    PFAV,TFAV,JPFAV,JTFAV,PFAM,TFAM,JPFAM,JTFAM)
C
C We now add the non-linear terms (stored in R1, R2, R3, R4
C and R5) to E1, E2, E3, E4 and E5
C
      N = NR*NH1
      CALL DAXPY( N, DELTAT, R1, INCX, E1, INCX )
C
      N = NR*NH2
      CALL DAXPY( N, DELTAT, R2, INCX, E2, INCX )
C
      N = NR*NH3
      CALL DAXPY( N, DELTAT, R3, INCX, E3, INCX )
C
      N = NR*NH4
      CALL DAXPY( N, DELTAT, R4, INCX, E4, INCX )
C
      N = NR*NH5
      CALL DAXPY( N, DELTAT, R5, INCX, E5, INCX )
C
C Solve the system of equations for to form predictor
C We use the LAPACK routine DGBTRS
C First, solve for poloidal velocity:
C
      N      = NR*NH1
      LDA    = 3*NBN + 1
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM1, LDA, IP1, E1, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for predictor solution of E1.'
      ENDIF
C
C Now, solve for toroidal velocity:
C
      N      = NR*NH2
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM2, LDA, IP2, E2, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for predictor solution of E2.'
      ENDIF
C
C Now, solve for temperature
C
      N      = NR*NH3
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM3, LDA, IP3, E3, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for predictor solution of E3.'
      ENDIF
C
C Now, solve for poloidal magnetic field:
C
      N      = NR*NH4
      LDA    = 3*NBN + 1
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM4, LDA, IP4, E4, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for predictor solution of E4.'
      ENDIF
C
C Now, solve for toroidal magnetic field:
C
      N      = NR*NH5
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM5, LDA, IP5, E5, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for predictor solution of E5.'
      ENDIF
C
C The predictor is now stored in E1, E2, E3, E4 and E5
C Calculate the norm of this solution.
C
      N     = NR*NH1
      AVEC1 = DNRM2( N, E1, INCX )
      N     = NR*NH2
      AVEC1 = AVEC1 + DNRM2( N, E2, INCX )
      N     = NR*NH3
      AVEC1 = AVEC1 + DNRM2( N, E3, INCX )
      N     = NR*NH4
      AVEC1 = AVEC1 + DNRM2( N, E4, INCX )
      N     = NR*NH5
      AVEC1 = AVEC1 + DNRM2( N, E5, INCX )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBCDT1: Predictor norm = ', AVEC1
C
C We then add the step^i forcing terms (non-linear)
C to the diffusion forcing terms. This will reduce
C the number of calculations done if more than one
C iteration is required.
C
      FAC = CFAC*DELTAT
C
      N   = NR*NH1
      CALL DAXPY( N, FAC, R1, INCX, D1, INCX )
C
      N   = NR*NH2
      CALL DAXPY( N, FAC, R2, INCX, D2, INCX )
C
      N   = NR*NH3
      CALL DAXPY( N, FAC, R3, INCX, D3, INCX )
C
      N   = NR*NH4
      CALL DAXPY( N, FAC, R4, INCX, D4, INCX )
C
      N   = NR*NH5
      CALL DAXPY( N, FAC, R5, INCX, D5, INCX )
C
C Now begin the loop around corrector iterations
C
      NOIT  = 0
 50   CONTINUE
      NOIT = NOIT + 1
C
      IF ( NOIT.GT.ITMX ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'Max iterations exceeded.'       
        NOIT = -1
        RETURN
      ENDIF
C
C Complete the poloidal velocity functions
C
      CALL PVVCPL( NR, NH1, NBN, E1, PVLC, PVRC )
C
C Calculate the non-linear forcing terms in
C vectors R1, R2, R3, R4 and R5.
C
      CALL OUBNLTF( XARR,ML1,MM1,ML2,MM2,ML3,MM3,ML4,MM4,ML5,MM5,
     1    E1,E2,E3,E4,E5,R1,R2,R3,R4,R5,SV3D,SV1Q,SV1S,SV2T,DSV3,
     2    VQ1,VS1,VT2,PA,DPA,F1,F2,F3,XSV,GAUX,GAUW,CVQ1,CVS1,CVT2Q,
     3    CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,CT2P,SV4Q,SV4S,SV5T,CVQ4,CVS4,
     4    CVT5Q,CVT5S,VQ4,VS4,VT5,VT4,VQ5,VS5,CT4P,CQ5T,CS5T,
     5    FTFPV,FTFTV,INFPV,INFTV,FTFPM,FTFTM,INFPM,INFTM,
     6    FTFPCV,FTFTCV,INFPCV,INFTCV,FTFPCM,FTFTCM,INFPCM,INFTCM,
     7    VGFA,IVGFA,SF2SA,ISF2S,
     8    PFAV,TFAV,JPFAV,JTFAV,PFAM,TFAM,JPFAM,JTFAM)
C
      N = NR*NH1
      CALL DCOPY( N, D1, INCX, E1, INCX )
C
      N = NR*NH2
      CALL DCOPY( N, D2, INCX, E2, INCX )
C
      N = NR*NH3
      CALL DCOPY( N, D3, INCX, E3, INCX ) 
C
      N = NR*NH4
      CALL DCOPY( N, D4, INCX, E4, INCX ) 
C
      N = NR*NH5
      CALL DCOPY( N, D5, INCX, E5, INCX ) 
C
C Now add the step^{i+1} non-lin. terms to E1, E2, E3,
C E4 and E5
C
      FAC = (1.0d0-CFAC)*DELTAT
C
      N   = NR*NH1
      CALL DAXPY( N, FAC, R1, INCX, E1, INCX )
C
      N   = NR*NH2
      CALL DAXPY( N, FAC, R2, INCX, E2, INCX )
C
      N   = NR*NH3
      CALL DAXPY( N, FAC, R3, INCX, E3, INCX )
C
      N   = NR*NH4
      CALL DAXPY( N, FAC, R4, INCX, E4, INCX )
C
      N   = NR*NH5
      CALL DAXPY( N, FAC, R5, INCX, E5, INCX )
C
C Solve the system of equations for iteration NOIT
C We use the LAPACK routine DGBTRS
C First, solve for poloidal velocity:
C
      N      = NR*NH1
      LDA    = 3*NBN + 1
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM1, LDA, IP1, E1, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for solution ',NOIT,' of E1.'
      ENDIF
C
C Now, solve for toroidal velocity:
C
      N      = NR*NH2
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM2, LDA, IP2, E2, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for solution ',NOIT,' of E2.'
      ENDIF
C
C Now, solve for temperature
C
      N      = NR*NH3
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM3, LDA, IP3, E3, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for solution ',NOIT,' of E3.'
      ENDIF
C
C Now, solve for poloidal magnetic field:
C
      N      = NR*NH4
      LDA    = 3*NBN + 1
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM4, LDA, IP4, E4, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for solution ',NOIT,' of E4.'
      ENDIF
C
C Now, solve for toroidal magnetic field:
C
      N      = NR*NH5
      CALL DGBTRS( TRANS, N, KL, KU, NRHS, AM5, LDA, IP5, E5, N,
     1             INFO )
C
      IF ( INFO.NE.0 ) THEN
        WRITE ( LLU, * ) 'Subroutine UBCDT1.'
        WRITE ( LLU, * ) 'DGBTRS returned INFO = ', INFO
        WRITE ( LLU, * ) 'for solution ',NOIT,' of E5.'
      ENDIF
C
      AOLD  = AVEC1
      N     = NR*NH1
      AVEC1 = DNRM2( N, E1, INCX )
      N     = NR*NH2
      AVEC1 = AVEC1 + DNRM2( N, E2, INCX )
      N     = NR*NH3
      AVEC1 = AVEC1 + DNRM2( N, E3, INCX )
      N     = NR*NH4
      AVEC1 = AVEC1 + DNRM2( N, E4, INCX )
      N     = NR*NH5
      AVEC1 = AVEC1 + DNRM2( N, E5, INCX )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'UBCDT1: Iteration ',NOIT,' norm = ', AVEC1
C
      DIFF = DABS( AVEC1 - AOLD )
      IF ( DIFF.LT.DTOL ) THEN
        IF ( LW ) WRITE ( IOUTF, * )
     1    'UBCDT1: Soln. converged iteration ', NOIT
          RETURN
      ENDIF
C
C Check to see if our norm appears to be getting bigger
C
      IF ( DIFF.GT.ODIFF .AND. NOIT.GT.3 ) THEN
        IF ( LW ) WRITE ( IOUTF, * )
     1    'UBCDT1: Soln. norm increasing: iteration ', NOIT
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
C Fills vectors R1, R2, R3, R4 and R5, with the following terms in   C
C the heat, induction and momentum equations:                        C
C                                                                    C
C F_Theta =  u . ( CB1 r + CB2 r^{-2} , 0 , 0 ) :                    C
C            - CC u . Grad ( Theta )            :                    C
C                                                                    C
C F_vort. =     CH curl ( Theta )               :                    C
C             - CG curl ( K x v )               :                    C
C             - CF curl ( v. Grad) v            :                    C
C             + CJ curl ( B. Grad) B            :                    C
C                                                                    C
C                                                                    C
C F_B     =   CM curl ( v x B )                 :                    C
C                                                                    C
C     SV1       : Dim ( NH1*NR ). Poloidal velocity initial soln.    C
C     SV2       : Dim ( NH2*NR ). Toroidal velocity initial soln.    C
C     SV3       : Dim ( NH3*NR ). Temperature initial soln.          C
C     SV4       : Dim ( NH4*NR ). Poloidal field initial soln.       C
C     SV5       : Dim ( NH5*NR ). Toroidal field initial soln.       C
C                                                                    C
C     R1        : Dim ( NH1*NR ). Toroidal vorticity forcing term    C
C     R2        : Dim ( NH2*NR ). Poloidal vorticity forcing term    C
C     R3        : Dim ( NH3*NR ). Temperature forcing term.          C
C     R4        : Dim ( NH4*NR ). Toroidal field forcing term        C
C     R5        : Dim ( NH5*NR ). Poloidal field forcing term        C
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
C     SV4Q      : Dim (         1, NH4*NR ) turns p. mag. to Q( r )  C
C     SV4S      : Dim ( 2*NBN + 1, NH4*NR ) turns p. mag. to S( r )  C
C     SV5T      : Dim (         1, NH5*NR ) turns t. mag. to T( r )  C
C                                                                    C
C     CVQ4      : Dim (         1, NH4*NR ) takes curl of Q( r )     C
C     CVS4      : Dim ( 2*NBN + 1, NH4*NR ) takes curl of S( r )     C
C     CVT5Q     : Dim (         1, NH5*NR ) takes curl of T( r )     C
C     CVT5S     : Dim ( 2*NBN + 1, NH5*NR ) takes curl of T( r )     C
C                                                                    C
C     CQ5T      : Dim (         1, NH5*NR ) takes curl of Q( r )     C
C     CS5T      : Dim ( 2*NBN + 1, NH5*NR ) takes curl of S( r )     C
C     CT4P      : Dim (         1, NH4*NR ) takes curl of T( r )     C
C                                                                    C
C     DSV3      : Dim ( NH3*NR ). Work array.                        C
C     VQ1       : Dim ( NH1*NR ). Work array.                        C
C     VS1       : Dim ( NH1*NR ). Work array.                        C
C     VT2       : Dim ( NH2*NR ). Work array.                        C
C     VT1       : Dim ( NH1*NR ). Work array.                        C
C     VQ2       : Dim ( NH2*NR ). Work array.                        C
C     VS2       : Dim ( NH2*NR ). Work array.                        C
C     VQ4       : Dim ( NH4*NR ). Work array.                        C
C     VS4       : Dim ( NH4*NR ). Work array.                        C
C     VT5       : Dim ( NH5*NR ). Work array.                        C
C     VT4       : Dim ( NH4*NR ). Work array.                        C
C     VQ5       : Dim ( NH5*NR ). Work array.                        C
C     VS5       : Dim ( NH5*NR ). Work array.                        C
C                                                                    C
C     F1        : Dim (2*NPHP). Work array for fourier transforming. C
C     F2        : Dim (2*NPHP). Work array for fourier transforming. C
C     F3        : Dim (2*NPHP). Work array for fourier transforming. C
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
C     FTFPM     : Transform coeffs Dim ( 3, NH4, NTHP ) {from RSDV2C}C
C     FTFTM     : Transform coeffs Dim ( 2, NH5, NTHP ) {from RSDV2C}C
C                                                                    C
C     FTFPCV    : Transform coeffs Dim ( 3, NH2, NTHP ) {from RSDV2C}C
C     FTFTCV    : Transform coeffs Dim ( 2, NH1, NTHP ) {from RSDV2C}C
C                                                                    C
C     FTFPCM    : Transform coeffs Dim ( 3, NH5, NTHP ) {from RSDV2C}C
C     FTFTCM    : Transform coeffs Dim ( 2, NH4, NTHP ) {from RSDV2C}C
C                                                                    C
C     VGFA      : Dim ( 3, NH3, NTHP ). Coeff.s from SF2VGC          C
C     IVGFA     : Dim ( 2, NH3 ). Locationss from SF2VGC             C
C                                                                    C
C     SF2SA     : Dim ( NH3, NTHP ). Coeff.s from SF2SDC             C
C     ISF2S     : Dim ( NH3 ) Indices from SF2SDC                    C
C                                                                    C
C     PFAV      : Dim ( 3, NH1, NTHP ). Coeffs from XSVSDC.          C
C     TFAV      : Dim ( 2, NH2, NTHP ). Coeffs from XSVSDC.          C
C                                                                    C
C     PFAM      : Dim ( 3, NH5, NTHP ). Coeffs from XSVSDC.          C
C     TFAM      : Dim ( 2, NH4, NTHP ). Coeffs from XSVSDC.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE OUBNLTF( XARR,ML1,MM1,ML2,MM2,ML3,MM3,ML4,MM4,ML5,MM5,
     1    SV1,SV2,SV3,SV4,SV5,R1,R2,R3,R4,R5,SV3D,SV1Q,SV1S,SV2T,DSV3,
     2    VQ1,VS1,VT2,PA,DPA,F1,F2,F3,XSV,GAUX,GAUW,CVQ1,CVS1,CVT2Q,
     3    CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,CT2P,SV4Q,SV4S,SV5T,CVQ4,CVS4,
     4    CVT5Q,CVT5S,VQ4,VS4,VT5,VT4,VQ5,VS5,CT4P,CQ5T,CS5T,
     5    FTFPV,FTFTV,INFPV,INFTV,FTFPM,FTFTM,INFPM,INFTM,
     6    FTFPCV,FTFTCV,INFPCV,INFTCV,FTFPCM,FTFTCM,INFPCM,INFTCM,
     7    VGFA,IVGFA,SF2SA,ISF2S,
     8    PFAV,TFAV,JPFAV,JTFAV,PFAM,TFAM,JPFAM,JTFAM)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Common block contents .....................C
      INTEGER          NR, NH1, NH2, NH3, NH4, NH5, NBN, M0,
     1                 IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX
      COMMON  /NDIMPARS/  NR, NH1, NH2, NH3, NH4, NH5, NBN, M0,
     1                    IOUTF, LH, ITMX, NOIT, NTHP, NPHP, NCMX
      DOUBLE PRECISION CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CJ, CK, CL, CM, CFAC, DELTAT, DTOL
      COMMON  /DPHYSPARS/ CA, CB1, CB2, CC, CD, CE, CF, CG, CH, CI,
     1                 CJ, CK, CL, CM, CFAC, DELTAT, DTOL
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION XARR( NR )
C
C Basic solution vector indices
C
      INTEGER          ML1( NH1 ), MM1( NH1 )
      INTEGER          ML2( NH2 ), MM2( NH2 )
      INTEGER          ML3( NH3 ), MM3( NH3 )
      INTEGER          ML4( NH4 ), MM4( NH4 )
      INTEGER          ML5( NH5 ), MM5( NH5 )
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
      DOUBLE PRECISION SV4Q(         1, NH4*NR )
      DOUBLE PRECISION SV4S( 2*NBN + 1, NH4*NR )
      DOUBLE PRECISION SV5T(         1, NH5*NR )
C
      DOUBLE PRECISION CVQ4(          1, NH4*NR )
      DOUBLE PRECISION CVS4(  2*NBN + 1, NH4*NR )
      DOUBLE PRECISION CVT5Q(         1, NH5*NR )
      DOUBLE PRECISION CVT5S( 2*NBN + 1, NH5*NR )
C
      DOUBLE PRECISION CQ5T(       1, NH5*NR )
      DOUBLE PRECISION CS5T( 2*NBN+1, NH5*NR )
      DOUBLE PRECISION CT4P(       1, NH4*NR )
C
C Initial solution vectors
C
      DOUBLE PRECISION SV1( NH1*NR ), SV2( NH2*NR ), SV3( NH3*NR ),
     1                 SV4( NH4*NR ), SV5( NH5*NR )
C
C Forcing term vectors
C
      DOUBLE PRECISION R1( NH1*NR ), R2( NH2*NR ), R3( NH3*NR ),
     1                 R4( NH4*NR ), R5( NH5*NR )
C
C Work arrays
C
      DOUBLE PRECISION DSV3( NH3*NR ), VQ1( NH1*NR ), VS1( NH1*NR ),
     1                 VT2( NH2*NR ), XSV( NCMX, NPHP, NTHP, NR),
     2                 VT1( NH1*NR ), VQ2( NH2*NR ), VS2( NH2*NR )
      DOUBLE PRECISION VQ4( NH4*NR ), VS4( NH4*NR ), VT5( NH5*NR ),
     1                 VT4( NH4*NR ), VQ5( NH5*NR ), VS5( NH5*NR )
C
C Arrays necessary for transform
C
      DOUBLE PRECISION F1(2*NPHP), F2(2*NPHP), F3(2*NPHP),
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
      DOUBLE PRECISION FTFPM( 3, NH4, NTHP ),
     1                 FTFTM( 2, NH5, NTHP )
C
      INTEGER          INFPM( 2, NH4 ),
     1                 INFTM( 2, NH5 )
C
      DOUBLE PRECISION FTFPCV( 3, NH2, NTHP ),
     1                 FTFTCV( 2, NH1, NTHP )
C
      INTEGER          INFPCV( 2, NH2 ),
     1                 INFTCV( 2, NH1 )
C
      DOUBLE PRECISION FTFPCM( 3, NH5, NTHP ),
     1                 FTFTCM( 2, NH4, NTHP )
C
      INTEGER          INFPCM( 2, NH5 ),
     1                 INFTCM( 2, NH4 )
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
      INTEGER          JPFAM( 2, NH5 ), JTFAM( 2, NH4 )
      DOUBLE PRECISION PFAM( 3, NH5, NTHP ), TFAM( 2, NH4, NTHP )
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
     1     CG.EQ.ZERO .AND. CC.EQ.ZERO .AND. CF.EQ.ZERO .AND.
     2     CJ.EQ.ZERO .AND. CM.EQ.ZERO ) RETURN
C
      IF ( IOUTF.EQ.0 ) THEN
        LW  = .FALSE.
      ELSE
        LW  = .TRUE.
      ENDIF
C
C Zero R1, R2, R3, R4 and R5.
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
C No need to zero R4 and R5 as they are done by DGBMV
C
c     N   = NR*NH4
c     CALL VECOP( R4, FAC, N, IOP )
C
c     N   = NR*NH5
c     CALL VECOP( R5, FAC, N, IOP )
C
C R1, R2, R3, R4 and R5 are all now zero ...
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
     1    'OUBNLTF: Heat and buoyancy terms added.'
      IF ( CG.EQ.ZERO .AND. CC.EQ.ZERO .AND. CF.EQ.ZERO .AND.
     1     CJ.EQ.ZERO .AND. CM.EQ.ZERO ) RETURN
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
     1    'OUBNLTF: Called DGBMV: dsv3:= sv3d . sv3.'
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
     1    'OUBNLTF: Called DGBMV: vq1:= sv1q . sv1.'
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
     1    'OUBNLTF: Called DGBMV: vs1:= sv1s . sv1.'
C
C Now calculate toroidal part of velocity.
C Use DGBMV to multiply SV2 by matrix SV2T to give VT2
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
     1    'OUBNLTF: Called DGBMV: vt2:= sv2t . sv2.'
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
      CALL RSDV2D( NTHP, NPHP, NR, ILNR, IRNR, NH1,
     1                  NH2,      NCMX, ICMR, ICMT, ICMP,
     2             VQ1, VS1, VT2, XSV, F1, F2, F3,      
     3             FTFPV, FTFTV, INFPV, INFTV, IFORMF )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called RSDV2D: XSV now contains velocity.'
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
      CALL SF2VGD( ILNR, IRNR, NR, NH3, ML3,          NTHP,
     1             NPHP, NCMX, ICMR, ICMT, ICMP, SV3, DSV3,
     2             XARR, XSV, F1, F2, F3, IVGFA, VGFA, IFORMF )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Grad( theta ) stored in XSV.'
C
C Now calculate scaloidal part of magnetic field.
C Use DGBMV to multiply SV4 by matrix SV4Q to give VQ4
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1
      M      = NR*NH4
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, SV4Q, LDA, SV4, INCX,
     1             BETA, VQ4, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called DGBMV: vq4:= sv4q . sv4.'
C
C Now calculate spheroidal part of magnetic field.
C Use DGBMV to multiply SV4 by matrix SV4S to give VS4
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1+2*NBN
      M      = NR*NH4
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, SV4S, LDA, SV4, INCX,
     1             BETA, VS4, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called DGBMV: vs4:= sv4s . sv4.'
C
C Now calculate toroidal part of magnetic field.
C Use DGBMV to multiply SV5 by matrix SV5T to give VT5
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1
      M      = NR*NH5
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, SV5T, LDA, SV5, INCX,
     1             BETA, VT5, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called DGBMV: vt5:= sv5t . sv5.'
C
C we now have the mag. field in QST format and so can transform
C into real space.
C
      ILNR = 2
      IRNR = NR - 1
C     .  put radial component in icmr = 10
      ICMR = 10
C     .  put theta component in icmt = 11
      ICMT = 11
C     .  put phi component in icmp = 12
      ICMP = 12
      CALL RSDV2D( NTHP, NPHP,     NR, ILNR, IRNR, NH4,
     1                  NH5,      NCMX, ICMR, ICMT, ICMP,
     2             VQ4, VS4, VT5, XSV, F1, F2, F3,
     3             FTFPM, FTFTM, INFPM, INFTM, IFORMF )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called RSDV2D: XSV now contains mag. field.'
C
C Magnetic field is now stored in XSV components 10-12.
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
     1    'OUBNLTF: Called DGBMV: vt1:= cvq1 . vq1.'
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
     1    'OUBNLTF: Called DGBMV: vt1:= cvq1 . vq1 + vt1.'
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
     1    'OUBNLTF: Called DGBMV: vq2:= cvt2q . vt2.'
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
     1    'OUBNLTF: Called DGBMV: vs2:= cvt2s . vt2.'
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
      CALL RSDV2D( NTHP, NPHP,     NR, ILNR, IRNR, NH2,
     1                  NH1,      NCMX, ICMR, ICMT, ICMP,
     2             VQ2, VS2, VT1, XSV, F1, F2, F3,
     3             FTFPCV, FTFTCV, INFPCV, INFTCV, IFORMF )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called RSDV2D: XSV now contains curl vel.'
C
C Now need to calculate curls of the scaloidal,
C spheroidal parts of the magnetic field which are stored
C in VQ4, VS4 and VT5 respectively ...
C The matrix CVQ4 multiplies VQ4 to give VT4 - which
C is a toroidal vector
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1
      M      = NR*NH4
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CVQ4, LDA, VQ4, INCX,
     1             BETA, VT4, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called DGBMV: vt4:= cvq4 . vq4.'
C
C Matrix CVS4 _adds_ the curl of VS4 to VT4
C
      ALPHA  = 1.0d0
      BETA   = 1.0d0
      LDA    = 1+2*NBN
      M      = NR*NH4
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CVS4, LDA, VS4, INCX,
     1             BETA, VT4, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called DGBMV: vt4:= cvq4 . vs4 + vt4.'
C
C The matrix CVT5Q multiplies VT5 to give VQ5 - which
C is a scaloidal vector
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1
      M      = NR*NH5
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CVT5Q, LDA, VT5, INCX,
     1             BETA, VQ5, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called DGBMV: vq5:= cvt5q . vt5.'
C
C The matrix CVT5S multiplies VT5 to give VS5 - which
C is a spheroidal vector
C
      ALPHA  = 1.0d0
      BETA   = 0.0d0
      LDA    = 1+2*NBN
      M      = NR*NH5
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CVT5S, LDA, VT5, INCX,
     1             BETA, VS5, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called DGBMV: vs5:= cvt5s . vt5.'
C
C The scaloidal vector of (curl B) is now in VQ5
C The spheroidal vector of (curl B) is now in VS5
C The toroidal vector of (curl B) is now in VT4
C So, calculate components in real space!
C
      ILNR = 2
      IRNR = NR - 1
C     .  put radial component in icmr = 13
      ICMR = 13
C     .  put theta component in icmt = 14
      ICMT = 14
C     .  put phi component in icmp = 15
      ICMP = 15
C
      CALL RSDV2D( NTHP, NPHP,     NR, ILNR, IRNR, NH5,
     1                  NH4,      NCMX, ICMR, ICMT, ICMP,
     2             VQ5, VS5, VT4, XSV, F1, F2, F3,
     3             FTFPCM, FTFTCM, INFPCM, INFTCM, IFORMF )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called RSDV2D: XSV now contains curl mag. field.'
C
C Now evaluate all of the non-linear operations
C within XSV. This is now a single subroutine call!
C
      CALL MDCXSE( NTHP, NPHP, NR, ILNR, IRNR, NCMX, XSV, GAUX,
     1             CG, CF, CJ )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Calculated non-linear terms in real space.'
C
C v . Grad( theta ) is stored in element 16 of XSV.
C Add back to R3 in spectral coefficients
C
      ILNR = 2
      IRNR = NR - 1
      ICM  = 16
c     FAC  = -1.0d0*CC
c     CALL SF2SDX( ILNR, IRNR, NR, LH, NH3, ML3, MM3, M0, NTHP,
c    1             NPHP, NCMX, ICM, PA, GAUW, XSV, R3, FAC, F1 )
c     CALL SF2SD2( ILNR, IRNR, NR, NH3, ML3, MM3, M0, NTHP,
c    1             NPHP, NCMX, ICM, XSV, R3, F1, SF2SA )
      CALL SF2SDD( ILNR, IRNR, NR, NH3, ML3,          NTHP,
     1             NPHP, NCMX, ICM, XSV, R3, F1, SF2SA,
     2             IFORMF, ISF2S )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Added v.Grad( theta ) terms to R3.'
C
C Transform the non-linear terms for the
C vorticity equation back into QST space.
C Note that we will use VQ1, VS1 and VT2
C
      ILNR = 2
      IRNR = NR - 1
C     .  radial component is stored in icmr = 17
      ICMR = 17
C     .  theta component is stored in icmt = 18
      ICMT = 18
C     .  phi component is stored in icmp = 19
      ICMP = 19
      CALL XSVSDD( NTHP, NPHP,     NR, ILNR, IRNR, NH1,
     1                  NH2,           NCMX, ICMR, ICMT, ICMP,
     2             VQ1, VS1, VT2, XSV, F1, F2, F3,
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
     1    'OUBNLTF: Called DGBMV: r1:= r1 + cq1t . vq1.'
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
     1    'OUBNLTF: Called DGBMV: r1:= r1 + cs1t . vs1.'
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
     1    'OUBNLTF: Called DGBMV: r2:= r2 + ct2p . vt2.'
C
C Transform the non-linear terms for the
C induction equation back into QST space.
C Note that we will use VQ5, VS5 and VT4
C
      ILNR = 2
      IRNR = NR - 1
C     .  radial component is stored in icmr = 20
      ICMR = 20
C     .  theta component is stored in icmt = 21
      ICMT = 21
C     .  phi component is stored in icmp = 22
      ICMP = 22
      CALL XSVSDD( NTHP, NPHP,     NR, ILNR, IRNR, NH5,
     1                  NH4,           NCMX, ICMR, ICMT, ICMP,
     2             VQ5, VS5, VT4, XSV, F1, F2, F3,
     3             PFAM, TFAM, JPFAM, JTFAM, IFORMF )
C
C Now add the curl of the scaloidal function to the
C toroidal forcing term of the magnetic field (R5)
C
      ALPHA  = CM
      BETA   = 0.0d0
      LDA    = 1
      M      = NR*NH5
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CQ5T, LDA, VQ5, INCX,
     1             BETA, R5, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called DGBMV: r5:= CM*cq5t . vq5.'
C
C Now add the curl of the spheroidal function to the
C toroidal forcing term of the magnetic field (R5)
C
      ALPHA  = CM
      BETA   = 1.0d0
      LDA    = 1+2*NBN
      M      = NR*NH5
      N      = M
      K      = NBN
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CS5T, LDA, VS5, INCX,
     1             BETA, R5, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called DGBMV: r5:= r5 + CM*cs5t . vs5.'
C
C Now add the curl of the toroidal function to the
C poloidal forcing term of the magnetic field (R4)
C
      ALPHA  = CM
      BETA   = 0.0d0
      LDA    = 1
      M      = NR*NH4
      N      = M
      K      = 0
      CALL DGBMV ( TRANS, M, N, K, K, ALPHA, CT4P, LDA, VT4, INCX,
     1             BETA, R4, INCY )
      IF ( LW ) WRITE ( IOUTF, * )
     1    'OUBNLTF: Called DGBMV: r4:= CM*ct4p . vt4.'
C
      RETURN
      END
C*********************************************************************
