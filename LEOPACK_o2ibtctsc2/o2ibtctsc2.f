C*********************************************************************
C                                                                    C
C Steve Gibbons - Inhomog. Boundary Thermal Conv. Time Step Code 2   C
C Fri Dec 15 09:16:48 WET 2000                                       C
C                                                                    C
C We read in a temperature function in the form of a .ints           C
C .vecs and .xarr file and impose as part of our convection solution C
C                                                                    C
C Reads one extra line of input from the o2ibtctsc input file.       C
C This line contains three integer values ...                        C
C NTSBSE: the number of time-steps between standard energy evaluationC
C NTSBLE: ditto for L-kinetic energy spectrum information.           C
C NTSBME: ditto for M-kinetic energy spectrum information.           C
C                                                                    C
C                                                                    C
C*********************************************************************
      PROGRAM o2ibtctsc2
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
     2                 NIV2MX, NIV3MX, NCMXX, NRMAXC, NIVCMX,
     3                 NHCMAX, NNDM
      PARAMETER      ( NRMAX = 50, NH1MAX = 1000, NH2MAX = 1000,
     1                 NH3MAX = 1000, NCMXX = 13, NRMAXC = 100,
     2                 NHMAX = NH1MAX+NH2MAX+NH3MAX, NHCMAX = 300,
     3                 NIVMAX = NHMAX*NRMAX, NIVCMX = NHCMAX*NRMAXC )
      PARAMETER      ( NIV1MX = NH1MAX*NRMAX,
     1                 NIV2MX = NH2MAX*NRMAX,
     2                 NIV3MX = NH3MAX*NRMAX )
      PARAMETER      ( NDCS = 4, NBNM = 3, NCFM = 2*NBNM+1,
     1                 NDRVM = 4, NNDM = 6 )
      INTEGER          LHMAX, NTHMAX, NPHMAX, NPMAX
      PARAMETER      ( LHMAX = 52, NTHMAX = 100, NPHMAX = 100,
     1                 NPMAX = (LHMAX+1)*(LHMAX+2)/2 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          MHTC( NHCMAX ), MHLC( NHCMAX ), MHMC( NHCMAX ),
     1                 IWRK( NNDM )
      DOUBLE PRECISION SVC( NIVCMX ), XARRC( NRMAXC ), WORK1( NNDM ),
     1                 WORK2( NNDM ), WORKM( NNDM, NNDM ),
     2                 DV3C( NIV3MX ), SV3C( NIV3MX ),
     3                 DSV3C( NIV3MX ), DKEARR( NHMAX )
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
     2                  ICOMP, IDIR, MMAX, IOPT, INARRC( 3 ),
     3                  NRC, NHC, NTSBSE, NTSBLE, NTSBME, LUNRG
      INTEGER           ITS, NTS, NTSBB, NTSBS, N, INC, NDIG, N1, N2,
     1                  K, L, M, LMINM, MMVAL, LLVAL
      DOUBLE PRECISION  PVLC( NBNM ), PVRC( NBNM ), DMONE, DPONE,
     1                  STIME, DTIME, SCAL, DTOL2, FAC, DKE( 2 ),
     2                  DTOTKE, DEATOT, DTORKE, DLOW
      CHARACTER *(120)  LINE, FNLOG, ROOT, FNAME, FNNRG
      CHARACTER *(10)   CHNUM
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IWRTE = 3, DMONE = -1.0d0, DPONE = 1.0d0,
     1            INC = 1, IFORM = 1, DLOW = 1.0d-9 )
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
      FNNRG(1:ILENRT) = ROOT(1:ILENRT)
      FNNRG(ILENRT+1:ILENRT+4) = '.nrg'
      FNNRG = FNNRG(1:ILENRT+4)
      CALL FOPEN( LUNRG, FNNRG, IWRTE )
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
C Now we need to input details about the inhomogeneous
C temperature
C
      PRINT *,' Enter name of inhomog. temp. harmonics file.'
 406  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 406
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 407
        ENDIF
      ENDDO
 407  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in harmonics file (ignoring boundary conditions!)
C
      CALL BIHFRD( NHC, NHCMAX, MHTC, MHLC, MHMC, LU, FNAME )
      WRITE ( LULOG, * ) 'Inhomog. harms. file read:'
      WRITE ( LULOG, * ) FNAME
      WRITE ( LULOG, * ) 'Number of inhomog. temp. harm.s = ', NHC
      WRITE ( LULOG, * )
C
      PRINT *,' Enter name of inhomog. temp. vector file.'
 408  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 408
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 409
        ENDIF
      ENDDO
 409  CONTINUE
      FNAME = LINE(1:ILEN)
C
      INARRC( 3 ) = NHC
      CALL SVFRD( INARRC, LU, NRMAXC, SVC, FNAME )
      WRITE ( LULOG, * ) 'Inhomog. vector file read:'
      WRITE ( LULOG, * ) FNAME
      WRITE ( LULOG, * )
C
      PRINT *,' Enter name of inhomog. radial spacings file.'
 410  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 410
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 411
        ENDIF
      ENDDO
 411  CONTINUE
      FNAME = LINE(1:ILEN)
C
      CALL XARRRD( NRC, NRMAXC, XARRC, LU, FNAME )
      WRITE ( LULOG, * ) 'Inhomog. radial spacings file read:'
      WRITE ( LULOG, * ) FNAME
      WRITE ( LULOG, * )
C
      IF ( INARRC( 2 ).NE.NRC ) THEN
        PRINT *,' Radial vector file and spacings file'
        PRINT *,' are incompatible. Program aborted.'
        STOP
      ENDIF
C
      PRINT *,' Enter scaling factor for temperature inhomog.'
 412  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 412
C
      READ ( LINE, * ) SCAL
      WRITE ( LULOG, * ) 'SCAL  = ', SCAL
      WRITE ( LULOG, * )
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
C Prepare the arrays containing terms for the inhomogeneous
C temperature - use routine ITVAAF
C
      DTOL2 = 1.0d-5
      CALL ITVAAF( NH3, NR, ML3, MM3, INARRC, MHTC, MHLC, MHMC,
     1             NNDM, IWRK, XARR, XARRC, SVC, CD, SCAL, DTOL2,
     2             DV3C, SV3C, DSV3C, WORK1, WORK2, WORKM )
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
C Enter info. for energy evaluations
C
      PRINT *,' Please enter NTSBSE NTSBLE NTSBME'
 215  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 215
      READ ( LINE, * )  NTSBSE, NTSBLE, NTSBME
      WRITE ( LULOG, * ) 'NTSBSE = ', NTSBSE
      WRITE ( LULOG, * ) 'NTSBLE = ', NTSBLE
      WRITE ( LULOG, * ) 'NTSBME = ', NTSBME
      WRITE ( LULOG, * )
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
C Current solution is in SV1, SV2 and SV3
C
      CALL IBTCT1( XARR,ML1,MM1,IPIV1,ML2,MM2,IPIV2,ML3,MM3,IPIV3,AM1,
     1  BM1,AM2,BM2,AM3,BM3,SV1,SV2,SV3,DV1,DV2,DV3,VE1,VE2,VE3,R1,R2,
     2  R3,PVLC,PVRC,SV3D,SV1Q,SV1S,SV2T,DSV3,VQ1,VS1,VT2,PA,DPA,FTF1,
     3  FTF2,FTF3,XSV,GAUX,GAUW,CVQ1,CVS1,CVT2Q,CVT2S,VT1,VQ2,VS2,
     4  CQ1T,CS1T,CT2P,DV3C,SV3C,DSV3C,FTFPV,FTFTV,INFPV,INFTV,
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
C       . Note there is no need to add the temperature now:
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
C           . Case of poloidal harmonic
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
C           . Case of toroidal harmonic
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
          ENDDO
C         .
C         . Write out total values
C         .
          WRITE ( LUNRG, 737 ) DTIME, DTOTKE, DEATOT, DTORKE
 737      FORMAT(1PD16.7,1PD16.7,1PD16.7,1PD16.7)
          CALL FLUSH( LUNRG )
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
 747      FORMAT('L= ',I5,1PD16.7,1PD16.7,1PD16.7,1PD16.7)
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
 757      FORMAT('M= ',I5,1PD16.7,1PD16.7,1PD16.7,1PD16.7)
          CALL FLUSH( LULOG )
C           .
          ENDDO
C         .   end of MMVAL = 0, MMAX loop
        ENDIF
C       . endif for M-spectrum
      ENDIF
C     . endif for outputting of energy spectra
C
C Return to perform next time-step
C
      IF ( ITS.EQ.NTS ) GOTO 999
      GOTO 50
C
 999  CONTINUE
      CALL FCLOSE( LUNRG, FNNRG, 'Error' )
      CALL FCLOSE( LULOG, FNLOG, 'Error' )
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
C     DV3C      : Dim ( NH3*NR ). Diffusion forcing from inhomog.    C
C     SV3C      : Dim ( NH3*NR ). Inhomog. temp. zero^th deriv.      C
C     DSV3C     : Dim ( NH3*NR ). Inhomog. temp. first deriv.        C
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
      SUBROUTINE IBTCT1( XARR,ML1,MM1,IPIV1,ML2,MM2,IPIV2,ML3,MM3,
     1  IPIV3,AM1,BM1,AM2,BM2,AM3,BM3,SV1,SV2,SV3,DV1,DV2,DV3,VE1,VE2,
     2  VE3,R1,R2,R3,PVLC,PVRC,SV3D,SV1Q,SV1S,SV2T,DSV3,VQ1,VS1,VT2,
     3  PA,DPA,FTF1,FTF2,FTF3,XSV,GAUX,GAUW,CVQ1,CVS1,CVT2Q,CVT2S,
     4  VT1,VQ2,VS2,CQ1T,CS1T,CT2P,DV3C,SV3C,DSV3C,FTFPV,FTFTV,INFPV,
     5  INFTV,FTFPCV,FTFTCV,INFPCV,INFTCV,VGFA,IVGFA,SF2SA,ISF2S,
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
C Inhomogeneous temperature terms
C
      DOUBLE PRECISION DV3C( NH3*NR ), SV3C( NH3*NR ), DSV3C( NH3*NR )
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
      IF ( LW ) WRITE ( IOUTF, * ) 'Entered IBTCT1.'
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
     1    'IBTCT1: Called DGBMV: dv1:= bm1 . sv1.'
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
     1    'IBTCT1: Called DGBMV: dv2:= bm2 . sv2.'
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
     1    'IBTCT1: Called DGBMV: dv3:= bm3 . sv3.'
C
C Now add on the term from the inhomogeneous part.
C
      CALL DAXPY( N, DELTAT, DV3C, INCX, DV3, INCX )
      IF ( LW ) WRITE ( IOUTF, * ) 
     1    'IBTCT1: Called DAXPY: dv3:= dv3 + deltat . dv3c.'
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
      CALL IBNLTF( XARR, ML1, MM1, ML2, MM2, ML3, MM3, SV1, SV2,
     1    SV3, R1, R2, R3, SV3D, SV1Q, SV1S, SV2T, DSV3, VQ1, VS1,
     2    VT2, PA, DPA, FTF1, FTF2, FTF3, XSV, GAUX, GAUW, CVQ1, CVS1,
     3    CVT2Q,CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,CT2P,SV3C,DSV3C,
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
        WRITE ( LLU, * ) 'Subroutine IBTCT1.'
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
        WRITE ( LLU, * ) 'Subroutine IBTCT1.'
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
        WRITE ( LLU, * ) 'Subroutine IBTCT1.'
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
     1    'IBTCT1: Predictor norm = ', AVEC1
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
        WRITE ( LLU, * ) 'Subroutine IBTCT1.'
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
      CALL IBNLTF( XARR, ML1, MM1, ML2, MM2, ML3, MM3, VE1, VE2,
     1    VE3, R1, R2, R3, SV3D, SV1Q, SV1S, SV2T, DSV3, VQ1, VS1,
     2    VT2, PA, DPA, FTF1, FTF2, FTF3, XSV, GAUX, GAUW, CVQ1, CVS1,
     3    CVT2Q,CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,CT2P,SV3C,DSV3C,
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
        WRITE ( LLU, * ) 'Subroutine IBTCT1.'
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
        WRITE ( LLU, * ) 'Subroutine IBTCT1.'
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
        WRITE ( LLU, * ) 'Subroutine IBTCT1.'
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
     1    'IBTCT1: Iteration ',NOIT,' norm = ', AVEC1
C
      DIFF = DABS( AVEC1 - AOLD )
      IF ( DIFF.LT.DTOL ) THEN
        IF ( LW ) WRITE ( IOUTF, * )
     1    'IBTCT1: Soln. converged iteration ', NOIT
          RETURN
      ENDIF
C
C Check to see if our norm appears to be getting bigger
C
      IF ( DIFF.GT.ODIFF .AND. NOIT.GT.3 ) THEN
        IF ( LW ) WRITE ( IOUTF, * )
     1    'IBTCT1: Soln. norm increasing: iteration ', NOIT
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
C subroutine Inhomogeneous Boundary Non-Linear Terms Fill ************
C            -             -        -   -      -     -    ************
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
C     SV3C      : Dim ( NH3*NR ). Inhomog. temp. zero^th deriv.      C
C     DSV3C     : Dim ( NH3*NR ). Inhomog. temp. first deriv.        C
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
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE IBNLTF( XARR, ML1, MM1, ML2, MM2, ML3, MM3, SV1, SV2,
     1    SV3, R1, R2, R3, SV3D, SV1Q, SV1S, SV2T, DSV3, VQ1, VS1,
     2    VT2, PA, DPA, FTF1, FTF2, FTF3, XSV, GAUX, GAUW, CVQ1, CVS1,
     3    CVT2Q,CVT2S,VT1,VQ2,VS2,CQ1T,CS1T,CT2P,SV3C,DSV3C,
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
C Inhomogeneous temperature terms
C
      DOUBLE PRECISION SV3C( NH3*NR ), DSV3C( NH3*NR )
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
Calculated by SF2SDC
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
     1    'IBNLTF: Called DGBMV: dsv3:= sv3d . sv3.'
C
C Add inhomogeneous part to temperature derivative
C
      ALPHA  = 1.0d0
      CALL DAXPY( N, ALPHA, SV3C, INCX, SV3, INCX )
      CALL DAXPY( N, ALPHA, DSV3C, INCX, DSV3, INCX )
C
C Now add buoyancy terms
C
      CALL NSVBTA( NR, NH3, ML3, MM3, NH1, ML1, MM1, SV3, R1, CH )
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
     1    'IBNLTF: Called DGBMV: vq1:= sv1q . sv1.'
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
     1    'IBNLTF: Called DGBMV: vs1:= sv1s . sv1.'
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
     1    'IBNLTF: Called DGBMV: vt2:= sv2t . sv2.'
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
     1    'IBNLTF: Called RSDV2D: XSV now contains velocity.'
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
     1    'IBNLTF: Grad( theta ) stored in XSV.'
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
     1    'IBNLTF: Called DGBMV: vt1:= cvq1 . vq1.'
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
     1    'IBNLTF: Called DGBMV: vt1:= cvq1 . vq1 + vt1.'
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
     1    'IBNLTF: Called DGBMV: vq2:= cvt2q . vt2.'
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
     1    'IBNLTF: Called DGBMV: vs2:= cvt2s . vt2.'
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
     1    'IBNLTF: Called RSDV2D: XSV now contains curl vel.'
C
C Now evaluate all of the non-linear operations
C within XSV. This is now a single subroutine call!
C
      CALL NMCXSE( NTHP, NPHP, NR, ILNR, IRNR, NCMX, XSV, GAUX,
     1             CG, CF )
C
      IF ( LW ) WRITE ( IOUTF, * )
     1    'IBNLTF: Calculated non-linear terms in real space.'
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
     1    'IBNLTF: Added v.Grad( theta ) terms to R3.'
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
     1    'IBNLTF: Called DGBMV: r1:= r1 + cq1t . vq1.'
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
     1    'IBNLTF: Called DGBMV: r1:= r1 + cs1t . vs1.'
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
     1    'IBNLTF: Called DGBMV: r2:= r2 + ct2p . vt2.'
C
      RETURN
      END
C*********************************************************************
C********************************************************************
C SUBROUTINE File OPEN **********************************************
C            -    ---- **********************************************
C Steve Gibbons 14.4.97 (Adapted from Dan Gordon's Code)            C
C Routine modified 15.3.99
C___________________________________________________________________C
C Opens a file with number LU, name FNAME, access OACCES  and       C
C a flag IRW to indicate whether the file is to be read or written  C
C to. ( IRW=1 ==> read only, IRW=2 ==> write but only if the file   C
C doesn't already exist, IRW=3 ==> write regardless of whether file C
C exists or not.)						    C
C___________________________________________________________________C
C Input Variables :-						    C
C ===============   						    C
C  Integer							    C
C  -------							    C
C     LU	: File number					    C
C     IRW	: Read / Write Flag 				    C
C                  = 1 for read only		                    C
C                  = 2 for write (provided that file doesn't exist. C
C                  = 3 for write (regardless of existence of file.  C
C                  = 4 for append status.                           C
C  Character							    C
C  ---------							    C
C     FNAME	: File name					    C
C___________________________________________________________________C
C Working Variables :-						    C
C =================   						    C
C  Character							    C
C  ---------							    C
C     OACCES 	: Access flag - should be set to 'OLD' for read     C
C                          and 'UNKNOWN' for write                  C
C     CONTYN    : For a yes/no to IWR = 2 option.		    C
C     LABEL     : Null string to pass into FNAMER option.	    C
C  Logical							    C
C  -------							    C
C     LEXIST	: Existence of file. File present <==> LEXIST=.TRUE.C
C___________________________________________________________________C
C Subroutines Used :-                                               C
C ================                                                  C
C     FNAMER	: For the case of IRW = 2; trying to write to an    C
C		   existing file. Used if alternative filename is   C
C                   asked for.					    C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE FOPEN ( LU, FNAME, IRW)
      IMPLICIT NONE
C___________________________________________________________________C
C Variable declarations - Parameters ...............................C
      INTEGER LU, IRW
      CHARACTER *(*) FNAME
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      LOGICAL LEXIST
      CHARACTER *(7) OACCES
      CHARACTER *(1) CONTYN
      CHARACTER *(1) LABEL
C___________________________________________________________________C
      IF ( LU.EQ.0 ) THEN
         PRINT *,' Subroutine FOPEN'
         PRINT *,' I bet you ve forgotten to set LU ...??'
         PRINT *,' Think again and come back when you have'
         PRINT *,' remembered that LU must be a non zero integer!'
         PRINT *,' See you later. Bye for now!!'
         STOP
      ENDIF
C************************
C temporary code : SJG Thu Jun  1 08:07:06 BST 2000
C The Linux compiler will not allow file opening with
C lu.ge.100, so the following line should prevent it.
C 
      IF ( LU.GT.99 ) THEN
         PRINT *,' Subroutine FOPEN'
         PRINT *,' LU = ', LU,' too large.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C************************
C--------------
 600  CONTINUE
      INQUIRE (FILE=FNAME, EXIST=LEXIST)
C Case of read only - LEXIST must be = .TRUE.
      IF ( IRW.EQ.1 ) THEN
         OACCES = 'OLD'
         IF ( .NOT. LEXIST ) THEN
            PRINT *,' Subroutine FOPEN. You are trying to open an'
            PRINT *,' old file which does not exist.'
            PRINT *,' Filename = ',FNAME,' LU= ',LU
            PRINT *,' Program aborted.'
            STOP
         ELSE
            GOTO 500
         ENDIF
      ENDIF
C Case of write to file provided that it doesn't exist
      IF ( IRW.EQ.2 ) THEN
         OACCES = 'UNKNOWN'
         IF ( LEXIST ) THEN
            PRINT *, ' Subroutine FOPEN. You are trying to write'
            PRINT *, ' to an existing file with IRW set to 2.'
            PRINT *,' Filename = ',FNAME,' LU= ',LU
            PRINT *,' Do you wish to give an alternative FNAME?'
            PRINT *,' Type y or n.'
            READ ( 5, 267) CONTYN
 267         FORMAT (A)
            IF (CONTYN.NE.'y'.AND.CONTYN.NE.'Y') THEN
               PRINT *, ' Program Aborted.'
               STOP
            ELSE
               LABEL=' '
               CALL FNAMER ( FNAME, LABEL )
               GOTO 600
            ENDIF
         ELSE
            GOTO 500
         ENDIF
      ENDIF
C Case of write to file regardless of the existence of file.
      IF ( IRW.EQ.3 ) THEN
         OACCES = 'UNKNOWN'
         GOTO 500
      ENDIF
C Treat appendment case
      IF ( IRW.EQ.4 ) THEN
         OACCES = 'UNKNOWN'
         OPEN ( UNIT=LU , FILE=FNAME , STATUS=OACCES,
     1          ACCESS='APPEND', ERR=999 )
         RETURN
      ENDIF
C___________________________________________________________________C
C All the IRW cases as of 14.4.97 have now been covered.
      PRINT *,' Subroutine FOPEN. IRW must be set to 1, 2, 3 or 4.'
      PRINT *,' Program aborted.'
      STOP

 500  CONTINUE
      OPEN ( UNIT=LU , FILE=FNAME , STATUS=OACCES, ERR=999 )
      RETURN

 999  PRINT *,' Subroutine FOPEN. Error in opening file ',FNAME
      STOP

      END
C********************************************************************
C___________________________________________________________________C
C*********************************************************************
C subroutine Single Harmonic Kinetic Energy Evaluation Routine *******
C            -      -        -       -      -          -       *******
C Steve Gibbons Thu Oct 28 08:51:22 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Returns 0.5d0 * \int_{volume} A.A dV where A is a single vector    C
C harmonic - the type being indicated by MHT( ih ).                  C
C                                                                    C
C This can be either magnetic energy or kinetic energy.              C
C                                                                    C
C The value for A = ( A_r, A_theta, A_phi ) is returned in DKE( 1 ). C
C                                                                    C
C The value for A = ( A_r, 0 , 0 ) is returned in DKE( 2 ).          C
C                                                                    C
C If MHT( ih ) does not correspond to a poloidal or toroidal         C
C vector harmonic, DKE( 1 ) and DKE( 2 ) are both returned zero.     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IH        : Number of harmonic to be evaluated.                C
C                                                                    C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1, 3. INDFUN = ( IR - 1 )*NH + IH       C
C                   IFORMF = 2, 4. INDFUN = ( IH - 1 )*NR + IR       C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C MHT defines what each scalar function in a solution vector         C
C represents.                                                        C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NDRVS     : Highest derivative stored in SVFDC.                C
C                (Not needed if we are doing tor --> pol but         C
C                 must be atleast 2 if doing pol --> tor )           C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NR*NH                       C
C                                                                    C
C     XARR      : Array of dimension ( NR )                          C
C                 XARR( j ) = element x_j                            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DKE       : Dimension ( 2 ).                                   C
C                                                                    C
C                 DKE( 1 ) is returned with the whole contribution   C
C                 to the kinetic energy from harmonic IH             C
C                                                                    C
C                 DKE( 2 ) is returned with kinetic energy of the    C
C                 radial component of harmonic IH                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHKEER( IH, NDCS, NR, INARR, MHT, MHL, MHP, NBN,
     1                   NDRVS, NDRVM, NFDCM, SV, XARR, DKE, SVFDC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IH, NDCS, NR, INARR( * ), MHT( * ), MHL( * ), MHP( * ),
     1        NBN, NDRVS, NDRVM, NFDCM
      DOUBLE PRECISION XARR( NR ), SV( * ), DKE( 2 ),
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IRAD, IHD, ITYPE, L, NR2, IS
      DOUBLE PRECISION COEF, DK, SQRLL1, PI, DERV( 2 ), D0F, D1F,
     1                 RAD, ER, ETOT, FAC
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR2 = INARR( 2 )
      IF ( NR2.NE.NR ) THEN
        PRINT *,' Subroutine SHKEER.'
        PRINT *,' NR         = ', NR
        PRINT *,' INARR( 2 ) = ', NR2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DKE( 1 ) = 0.0d0
      DKE( 2 ) = 0.0d0
C     .
      ITYPE = MHT( IH )
      L     = MHL( IH )
      IS    = MHP( IH )
      DK    = SQRLL1( L )
C     .
C     . First do poloidal harmonics
C     .
      IF ( ITYPE.EQ.1 .OR. ITYPE.EQ.4 ) THEN
        IHD = 1
        DO IRAD = 1, NR
          RAD = XARR( IRAD )
C         .
          IF ( IRAD.EQ.1 ) THEN
            COEF = 0.5d0*(  XARR( 2 ) - XARR( 1 )  )
          ENDIF
C         .
          IF ( IRAD.GT.1 .AND. IRAD.LT.NR ) THEN
            COEF = 0.5d0*(  XARR( IRAD+1 ) - XARR( IRAD-1 )  )
          ENDIF
C         .
          IF ( IRAD.EQ.NR ) THEN
            COEF = 0.5d0*(  XARR( NR ) - XARR( NR-1 )  )
          ENDIF
C         .
          CALL ASVDR( SV, IRAD, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                NDRVM, DERV, INARR, SVFDC, NDCS )
          D0F = DERV( 1 )
          D1F = DERV( 2 )
C         .
C         . Add contribution to ER and ETOT
C         .
          ER   = D0F*D0F
          ETOT = DK*DK*D0F*D0F + (D0F + RAD*D1F)*(D0F + RAD*D1F)
C         .
C         . Add to cumulative integral
C         .
          DKE( 1 ) = DKE( 1 ) + ETOT*COEF
          DKE( 2 ) = DKE( 2 ) + ER*COEF
C         .
        ENDDO
C       .
C       . Finally, multiply by leading factors
C       .
        FAC = 2.0d0*PI*DK*DK/(2.0d0*DBLE( L ) + 1.0d0)
        DKE( 1 ) = DKE( 1 )*FAC
        DKE( 2 ) = DKE( 2 )*FAC*DK*DK
C       .
      ENDIF
C     .
C     . Now do toroidal harmonics
C     .
      IF ( ITYPE.EQ.2 .OR. ITYPE.EQ.5 ) THEN
        IHD = 0
        DO IRAD = 1, NR
          RAD = XARR( IRAD )
C         .
          IF ( IRAD.EQ.1 ) THEN
            COEF = 0.5d0*(  XARR( 2 ) - XARR( 1 )  )
          ENDIF
C         .
          IF ( IRAD.GT.1 .AND. IRAD.LT.NR ) THEN
            COEF = 0.5d0*(  XARR( IRAD+1 ) - XARR( IRAD-1 )  )
          ENDIF
C         .
          IF ( IRAD.EQ.NR ) THEN
            COEF = 0.5d0*(  XARR( NR ) - XARR( NR-1 )  )
          ENDIF
C         .
          CALL ASVDR( SV, IRAD, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                NDRVM, DERV, INARR, SVFDC, NDCS )
          D0F = DERV( 1 )
C         .
C         . Add contribution to ER and ETOT
C         .
          ETOT = RAD*RAD*D0F*D0F
C         .
C         . Add to cumulative integral
C         .
          DKE( 1 ) = DKE( 1 ) + ETOT*COEF
C         .
        ENDDO
C       .
C       . Finally, multiply by leading factors
C       .
        FAC = 2.0d0*PI*DK*DK/(2.0d0*DBLE( L ) + 1.0d0)
        DKE( 1 ) = DKE( 1 )*FAC
      ENDIF
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine HarMonic File ReaD **************************************
C            -  -     -    -  - **************************************
C Steve Gibbons Sat Nov 13 12:50:33 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Reads in the integer indices of spherical harmonic sets incl.      C
C the appropriate boundary conditions from a file.                   C
C It carefully assesses the boundary conditions and creates or       C
C appends the MHIBC, MHOBC and LARR arrays which must be sent to     C
C SVFDCF to calculate the finite difference coefficients.            C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Number of vector spherical harmonics. (Output)     C
C     NHMAX     : Maximum number of vector spherical harmonics.      C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NHMAX          C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NHMAX          C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NHMAX          C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NHMAX          C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     NCUDS     : Number of currently used finite diff. schemes.     C
C     NDCS      : Number of distinct finite difference schemes.      C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). Governs behaviour at inner     C
C                  boundary for finite diff. scheme ( is )           C
C                                                                    C
C  MHIBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( is ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                        where L = MHL( ih )                         C
C                                                                    C
C     MHOBC     : Dimension ( NDCS ). Governs behaviour at outer     C
C                  boundary for finite diff. scheme ( is )           C
C                                                                    C
C  MHOBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( is ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                        where L = MHL( ih )                         C
C                                                                    C
C     LARR      : Spherical harmonic degree, L. Dim. ( NDCS ).       C
C                 Array to be passed to SVFDCF.                      C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE HMFRD( NH, NHMAX, MHT, MHL, MHM, MHP, NCUDS, NDCS,
     1                  MHIBC, MHOBC, LARR, LU, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, NHMAX, MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     1        MHP( NHMAX ), NCUDS, NDCS,
     1        MHIBC( NDCS ), MHOBC( NDCS ), LARR( NDCS ), LU
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH, IWR, IS, IIBF, IOBF, ND1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
C Read in number of spherical harmonics
C
      READ ( LU, * ) NH
C
C Check value of NH
C
      IF ( NH.GT.NHMAX ) THEN
        PRINT *, ' Subroutine HMFRD.'
        PRINT *, ' From file, NH = ', NH
        PRINT *, ' NHMAX = ', NHMAX
        PRINT *, ' Program aborted.'
        STOP
      ENDIF
C
C     Now loop around the harmonics and read in their
C     properties
C
      DO IH = 1, NH
        READ ( LU, * ) MHT( IH ), MHL( IH ), MHM( IH ), IIBF, IOBF
C       .
C       . We now need to see whether we already have a
C       . suitable finite difference scheme for this outcome.
C       .
        ND1 = MIN( NCUDS, NDCS )
        DO IS = 1, ND1
         IF ( MHIBC( IS ).EQ.IIBF .AND. MHOBC( IS ).EQ.IOBF ) THEN
C          .
C          . Simple case where LARR needn't be applied
C          .
           IF ( IIBF.NE.7 .AND. IOBF.NE.7 ) THEN
             MHP( IH ) = IS
             GOTO 60
           ENDIF
C          .
C          . OK - check to see if LARR value is correct
C          .
           IF ( LARR( IS ).EQ.MHL( IH ) ) THEN
             MHP( IH ) = IS
             GOTO 60
           ENDIF
C          .
         ENDIF
        ENDDO
C       .
C       . OK - we don't already have this diff. scheme
C       .
        NCUDS = NCUDS + 1
        IF ( NCUDS.GT.NDCS ) THEN
          PRINT *, ' Subroutine HMFRD.'
          PRINT *, NDCS,' schemes are available.'
          PRINT *, ' and these have now been exhausted.'
          PRINT *, ' Program aborted.'
          STOP
        ENDIF
C       .
        MHIBC( NCUDS ) = IIBF
        MHOBC( NCUDS ) = IOBF
        IF ( IIBF.NE.7 .AND. IOBF.NE.7 ) THEN
          LARR( NCUDS ) = 0
        ELSE
          LARR( NCUDS ) = MHL( IH )
        ENDIF
        MHP( IH ) = NCUDS
C       .
 60     CONTINUE
C       .
      ENDDO
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
C Abort program if more finite difference schemes were
C requested by the harmonic sets than were allowed
C
      IF ( NCUDS.LT.NDCS ) THEN
        DO IS = NCUDS + 1, NDCS
          LARR( IS ) = -1
        ENDDO
      ENDIF
C
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine Solution Vector File ReaD *******************************
C            -        -      -    -  - *******************************
C Steve Gibbons Sat Nov 13 14:30:24 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Reads in a solution vector from a file.                            C
C IMPORTANT. The set of spherical harmonics must ALREADY BE KNOWN    C
C before calling SVFRD. If SVFRD reads in a file with NH different   C
C to the NH input in INARR( 3 ), an error will be reported.          C
C                                                                    C
C The only check on NR is that it is not greater than NRMAX.         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1, 3. INDFUN = ( IR - 1 )*NH + IH       C
C                   IFORMF = 2, 4. INDFUN = ( IH - 1 )*NR + IR       C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C     NRMAX     : Maximum permitted radial grid nodes.               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Dim ( * ) but length atleast NR*NH.                C
C                  Solution vector defined by INARR.                 C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVFRD( INARR, LU, NRMAX, SV, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), LU, NRMAX
      CHARACTER *(*) FNAME
      DOUBLE PRECISION SV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IWR, ILEN, IFORMF, NR, NH, NH1, IFORM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NH     = INARR( 3 )
C
C Open file for reading
C
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
C Read iformf, nr, nh, iform
C  
       READ ( LU, * ) IFORMF, NR, NH1, IFORM
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine SVFRD.'
        PRINT *,' IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check that value of NH is legal
C
      IF ( NH.NE.NH1 ) THEN
        PRINT *,' Subroutine SVFRD.'
        PRINT *,' INARR( 3 ) = ', NH
        PRINT *,' File contains ',NH1,' harmonics.'
        PRINT *,' Load in correct indices before calling SVFRD.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check that value of NR is legal
C
      IF ( NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine SVFRD.'
        PRINT *,' In file, NR = ', NR
        PRINT *,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INARR( 1 ) = IFORMF
      INARR( 2 ) = NR
C
      ILEN   = NR*NH
C
C OK, so read SV values ...
C
      IF ( IFORM.EQ.1 ) READ ( LU, 41 ) ( SV( I ), I = 1, ILEN )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 41   FORMAT(5(1PD16.7))
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine X value ARRay ReaD **************************************
C            -       ---   -  - **************************************
C Steve Gibbons Fri Nov 12 07:58:51 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Reads in the XARR array of abscissae from a file.                  C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                 Note that NR is not checked for correspondence to  C
C                 any other value - merely for being not greater     C
C                 than NRMAX.                                        C
C                                                                    C
C     NRMAX     : Maximum number of radial grid nodes.               C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim ( * ) but length atleast NR. Location of       C
C                  radial grid nodes.                                C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XARRRD( NR, NRMAX, XARR, LU, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NRMAX, LU
      CHARACTER *(*) FNAME
      DOUBLE PRECISION XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IWR, IFORM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
C Read number of radial grid nodes
C
      READ ( LU, * ) NR, IFORM
C
C Check NR
C
      IF ( NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine XARRRD.'
        PRINT *,' From file, NR = ', NR
        PRINT *,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine XARRRD.'
        PRINT *,' From file, IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C OK, so read X values ...
C
      IF ( IFORM.EQ.1 ) READ ( LU, 41 ) ( XARR( I ), I = 1, NR )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 41   FORMAT(5(1PD16.7))
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Solution Vector Finite Difference Coefficients Form *****
C            -        -      -      -          -            -    *****
C Steve Gibbons Fri Oct 22 09:33:36 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  If XARR is an array of length ( NR ) such that the j^{th}         C
C  element is the value of x_j, then SVFDCF builds an array          C
C  SVFDC of dimension ( NFDCM, NR, NDRVM+1, NDCS ) such that for a   C
C  node number, j, the ND^{th} derivative of radial function         C
C  given f_{IH} ( x ) will be given by                               C
C                                                                    C
C  f_{IH}^{ND}( x_j ) =                                              C
C         \sum_{i=LN}^{RN} SVFDC ( IRAD, j, ND+1, K ) f_{IH} ( x_i ) C
C                                                                    C
C  where LN (the left node)  = MAX( NLMR, j - NBN ) and              C
C        RN (the right node) = MIN( NRMC, j + NBN ),                 C
C                                                                    C
C  IRAD = i - j + NBN + 1 and K = MHP( ih ).                         C
C                                                                    C
C  NDCS is the number of distinct sets of coefficients required for  C
C  different types of harmonics (in generally will be considerably   C
C  smaller than the number of harmonics, NH).                        C
C                                                                    C
C  The arrays MHIBC and MHOBC instruct SVFDCF how to manipulate      C
C  the finite difference coefficients at the boundaries.             C
C                                                                    C
C  MHIBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( ih ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                        where L = LARR( ih )                        C
C                                                                    C
C  Similarly, at the outer boundary:-                                C
C                                                                    C
C  MHOBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( ih ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                        where L = LARR( ih )                        C
C                                                                    C
C  The elements of this array are filled in from j = NLMR            C
C  to j = NRMR ( number of the left most node and number of the      C
C  right most node ) - other rows are left unreferred to.            C
C  This is incase a higher order derivative is required for          C
C  central nodes than boundary nodes; in which case SVFDCF must      C
C  be called for the remaining nodes with modified parameters.       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
C                                                                    C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C     NLMR      : This is the lowest j for which the terms are       C
C                  calculated for SVFDC( i, j, ND+1, K )             C
C     NRMR      : This is the highest j for which the terms are      C
C                  calculated for SVFDC( i, j, ND+1, K )             C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). See above.                     C
C     MHOBC     : Dimension ( NDCS ). See above.                     C
C     LARR      : Spherical harmonic degree, L. Dim. ( NDCS ).       C
C                 This value is only useful for harmonics            C
C                 when calculating derivatives for magnetic fields.  C
C                 If LARR( ih ) = -1, the harmonic is ignored        C
C                 completely.                                        C
C                                                                    C
C     NCFM      : Leading order of working coefficient matrix.       C
C                 Must be atleast (2*NBN + 1) where NBN is the       C
C                 maximum number of nodes on either side of the      C
C                 central node.                                      C
C     NFDCM     : Leading order of the array SVFDC.                  C
C                 This must be atleast (2*NBN + 1)                   C
C     NDRVS     : Number of derivatives required.                    C
C                  This will be limited by the available bandwidth.  C
C                                                                    C
C                  Let NLCS = NLMR - 1    and let                    C
C                      NRCS = NR - NRMR                              C
C                                                                    C
C                  Now, let I = MIN( NLCS, NRCS) + NBN               C
C                                                                    C
C                  then NDRVS must be no greater than I.             C
C                  This is checked for.                              C
C     NDRVM     : Maximum number of derivatives allowed.             C
C                                                                    C
C     IWORK     : Integer work array. Dimension ( NCFM )             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                                                                    C
C     COEFM1    : Coefficient work array. Dimension ( NCFM, NCFM )   C
C     COEFM2    : Coefficient work array. Dimension ( NCFM, NCFM )   C
C     WORK1     : Coefficient work array. Dimension ( NCFM )         C
C     WORK2     : Coefficient work array. Dimension ( NCFM )         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVFDCF( NR, NDCS, NBN, NLMR, NRMR, MHIBC, MHOBC,
     1                   LARR, NCFM, NFDCM, NDRVS, NDRVM, XARR,
     2                   IWORK, SVFDC, COEFM1, COEFM2, WORK1, WORK2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NDCS, NBN, NLMR, NRMR, MHIBC( NDCS ), MHOBC( NDCS ),
     1        LARR( NDCS ), NCFM, NFDCM, NDRVS, NDRVM, 
     1        IWORK( NCFM )
      DOUBLE PRECISION XARR( NR ), SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                 COEFM1( NCFM, NCFM ), COEFM2( NCFM, NCFM ),
     2                 WORK1( NCFM ), WORK2( NCFM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NDER, IRAD, NLCS, NRCS, NLN, NRN, IDCS, L,
     1        NNDS, INDS, I, INODE, NSNIB, NSNOB, NALF, NARF,
     2        IIBC, IOBC, ND1
C
C nsnib is the number of special nodes on the inner boundary
C nsnob is the number of special nodes on the outer boundary
C
      DOUBLE PRECISION DZERO, X0, EMMULT, FAC
      PARAMETER ( DZERO = 0.0d0 )
C
      LOGICAL OCHNGE
C
C ochnge is .TRUE. when the boundary conditions
C play a part in the finite difference coefficients
C and .FALSE. otherwise
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     .
C     . Check the values of integers ...
C     .
      IF ( NDRVS.GT.NDRVM ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NDRVS  = ',NDRVS
         PRINT *,' NDRVM  = ',NDRVM
         STOP
      ENDIF
C     .
      NLCS = NLMR - 1
      NRCS = NR - NRMR
C     . 
C     . Check that sufficient points are allowed
C     . for the derivatives ...
C     . 
      IF ( (NR-1).LT.(NBN+1) ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NBN  = ', NBN
         PRINT *,' NLMR = ', NLMR
         PRINT *,' NRMR = ', NRMR
         PRINT *,' Insufficient nodes for differencing.'
         STOP
      ENDIF
C     . 
      IF ( NLCS.LT.0 ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NLMR = ', NLMR
         STOP
      ENDIF
C     . 
      IF ( NRCS.LT.0 ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NRMR = ', NRMR
         STOP
      ENDIF
C     .
      I = MIN( NLCS, NRCS) + NBN
C     .
      IF ( NDRVS.GT.I ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' You have requested deriv.s to order ',NDRVS
         PRINT *,' At one node, you have only', I
         PRINT *,' side nodes to use for differencing.'
         STOP
      ENDIF
C     .
      I = 2*NBN + 1
C     .
      IF ( NCFM.LT.I ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NCFM = ', NCFM
         PRINT *,' NBN  = ', NBN
         STOP
      ENDIF
C     .
      IF ( NFDCM.LT.I ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NFDCM = ', NFDCM
         PRINT *,' NBN  = ', NBN 
         STOP
      ENDIF
C     .
      IF ( NLMR.GT.NRMR ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NLMR = ', NLMR
         PRINT *,' NRMR = ', NRMR
         STOP
      ENDIF
C     .
C     . Whew ... all input parameters seem to be o.k.
C     . Loop around the requested nodes
C     .
      DO IRAD = NLMR, NRMR
C      .
C      . Loop around the different 'harmonic forms'
C      .
       DO IDCS = 1, NDCS
C
C Check for flag to ignore this harmonic
C
        IF ( LARR( IDCS ).EQ.-1 ) GOTO 50
C
C We now need to check which boundary condition is
C required. Make sure that it is valid.
C
        IF ( MHIBC( IDCS ).LT.1 .AND. MHIBC( IDCS ).GT.7 ) THEN
          PRINT *,' Subroutine SVFDCF.'
          PRINT *,' MHIBC(',IDCS,') = ', MHIBC( IDCS )
          STOP
        ENDIF
C
C O.k. inner b.c. is fine. Now need to 
C see how many points this effects.
C
        IF ( MHIBC( IDCS ).EQ.1 ) THEN
          NSNIB = 0
        ENDIF
C
        IF ( MHIBC( IDCS ).EQ.2 .OR. MHIBC( IDCS ).EQ.3 .OR.
     1       MHIBC( IDCS ).EQ.6 .OR. MHIBC( IDCS ).EQ.7     ) THEN
          NSNIB = 1
        ENDIF
C
        IF ( MHIBC( IDCS ).EQ.4 .OR. MHIBC( IDCS ).EQ.5 ) THEN
          NSNIB = 2
        ENDIF
C
        IF ( MHOBC( IDCS ).LT.1 .AND. MHOBC( IDCS ).GT.7 ) THEN
          PRINT *,' Subroutine SVFDCF.'
          PRINT *,' MHOBC(',IDCS,') = ', MHOBC( IDCS )
          STOP
        ENDIF
C
C O.k. outer b.c. is fine. Now need to
C see how many points this effects.
C
        IF ( MHOBC( IDCS ).EQ.1 ) THEN
          NSNOB = 0
        ENDIF
C
        IF ( MHOBC( IDCS ).EQ.2 .OR. MHOBC( IDCS ).EQ.3 .OR.
     1       MHOBC( IDCS ).EQ.6 .OR. MHOBC( IDCS ).EQ.7     ) THEN
          NSNOB = 1 
        ENDIF
C
        IF ( MHOBC( IDCS ).EQ.4 .OR. MHOBC( IDCS ).EQ.5 ) THEN
          NSNOB = 2 
        ENDIF
C
        DO NDER = 0, NDRVS
          ND1 = NDER + 1
          DO I = 1, NFDCM
            SVFDC( I, IRAD, ND1, IDCS ) = DZERO
          ENDDO
        ENDDO
C
C irad is the node for which we want to calculate
C our coefficients
C
        NLCS = IRAD - 1
        NRCS = NR - IRAD
C
C we wish to calculate NLN (number of left nodes)
C and NRN ( number of right nodes )
C NNDS ( total number of nodes) is then NLN + NRN + 1 ...
C
        NLN = MIN( NBN, NLCS )
        NRN = MIN( NBN, NRCS )
C
        NNDS = NLN + NRN + 1
C
C We must work out how many nodes are affected
C to the left and the right. NALF and NARF
C are respectively the number of affected nodes
C to the left and right.
C
        IF ( (IRAD-NLN).GT.NSNIB ) NALF = 0
        IF ( (IRAD-NLN).EQ.NSNIB ) NALF = 1
        IF ( (IRAD-NLN).LT.NSNIB ) NALF = 2
C
        IF ( (IRAD+NRN).LT.(NR+1-NSNOB) ) NARF = 0
        IF ( (IRAD+NRN).EQ.(NR+1-NSNOB) ) NARF = 1
        IF ( (IRAD+NRN).GT.(NR+1-NSNOB) ) NARF = 2
C
        IF ( NALF.EQ.0 .AND. NARF.EQ.0 ) THEN
          OCHNGE = .FALSE.
        ELSE
          OCHNGE = .TRUE.
        ENDIF
        IF ( .NOT. OCHNGE ) GOTO 51
C       .
C       . OK - we need to form a matrix COEFM2 such that
C       . the correct coeffcients are given when
C       . COEFM1 is multiplied by COEFM2
C       .
        L    = LARR( IDCS )
        IIBC = MHIBC( IDCS )
        IOBC = MHOBC( IDCS )
C       .
        CALL LDGNMF( NR, NNDS, NALF, NARF, L, IIBC, IOBC, NCFM,
     1             XARR, COEFM2, COEFM1, WORK1, WORK2, IWORK )
C       .
 51     CONTINUE
        X0 = XARR( IRAD )
        DO INDS = 1, NNDS
          INODE = IRAD - NLN - 1 + INDS
          WORK1( INDS ) = XARR( INODE )
        ENDDO
C
C Now ready to calculate the coefficients
C
        CALL GFDCFD( X0, WORK1, NNDS, COEFM1, NCFM, 
     1               IWORK, WORK2 )
C
C coefm matrix should now contain the coeff.s
C
        IF ( OCHNGE ) THEN
C        .
C        . Our coefficients are modified
C        . by the boundary conditions
C        .
         DO NDER = 0, NDRVS
          ND1 = NDER + 1
          DO INDS = 1, NNDS
            INODE = INDS - NLN + NBN
            FAC = EMMULT( ND1, INDS, NCFM, NCFM, NNDS,
     1                      COEFM1, COEFM2 )
            SVFDC( INODE, IRAD, ND1, IDCS ) = FAC
          ENDDO
         ENDDO
        ELSE
C        .
C        . Our coefficients are not modified
C        . by the boundary conditions
C        .
         DO NDER = 0, NDRVS
          ND1 = NDER + 1
          DO INDS = 1, NNDS
            INODE = INDS - NLN + NBN
            SVFDC( INODE, IRAD, ND1, IDCS ) = 
     1                      COEFM1( ND1, INDS )
          ENDDO
         ENDDO
C        .
        ENDIF
C
 50    CONTINUE
       ENDDO
C      .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Finite Difference Coefficient Matrix BuilD **************
C            -      -          -           -      -   - **************
C Steve Gibbons Tue Sep 21 09:25:54 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  If XARR is an array of length ( NR ) such that the j^{th}         C
C  element is the value of x_j, then FDCMBD builds an array          C
C  FDCM of dimension ( NFDCM, NR, NDRVS ) such that for a given      C
C  node number, j, the ND^{th} derivative of a function f( x )       C
C  will be given by                                                  C
C                                                                    C
C  f^{ND}( x_j ) = \sum_{i=LN}^{RN} FDCM( IRAD, j, ND ) f( x_i )     C
C                                                                    C
C  where LN (the left node)  = MAX( NLMC, j - NBN ) and              C
C        RN (the right node) = MIN( NRMC, j + NBN )                  C
C                                                                    C
C  and IRAD = i - j + NBN + 1                                        C
C                                                                    C
C  NLMC and NRMC are respectively the left most and right most       C
C  nodes (columns) which may be used to obtain a difference formula. C
C  In most matrix applications NLMC = 1 and NRMC = NR, although      C
C  when differentiating a vector it may be necessary to omit an      C
C  extreme point; for instance when this would result in a division  C
C  by zero.                                                          C
C                                                                    C
C  The elements of this array are filled in from j = NLMN            C
C  to j = NRMN ( number of the left most node and number of the      C
C  right most node ) - other rows are left unreferred to.            C
C  This is incase a higher order derivative is required for          C
C  central nodes than boundary nodes; in which case FDCMBD must      C
C  be called for the remaining nodes with modified parameters.       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     FDCM      : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM ).                   C
C                                                                    C
C     COEFM     : Coefficient work array. Dimension ( NCFM, NCFM )   C
C     WORK1     : Coefficient work array. Dimension ( NCFM )         C
C     WORK2     : Coefficient work array. Dimension ( NCFM )         C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C     NLMN      : This is the lowest j for which the terms are       C
C                  calculated for FDCM( i, j, ND )                   C
C     NRMN      : This is the highest j for which the terms are      C
C                  calculated for FDCM( i, j, ND )                   C
C     NLMC      : This is the lowest i for which the terms are       C
C                  calculated for FDCM( i, j, ND )                   C
C     NRMC      : This is the highest i for which the terms are      C
C                  calculated for FDCM( i, j, ND )                   C
C                                                                    C
C     NCFM      : Leading order of working coefficient matrix.       C
C                 Must be atleast (2*NBN + 1) where NBN is the       C
C                 maximum number of nodes on either side of the      C
C                 central node.                                      C
C     NFDCM     : Leading order of the array FDCM.                   C
C                 This must be atleast (2*NBN + 1)                   C
C     NDRVS     : Number of derivatives required.                    C
C                  This will be limited by the available bandwidth.  C
C                                                                    C
C                  Let NLCS = NLMN - NLMC and let                    C
C                      NRCS = NRMC - NRMN                            C
C                                                                    C
C                  Now, let I = MIN( NLCS, NRCS) + NBN               C
C                                                                    C
C                  then NDRVS must be no greater than I.             C
C                  This is checked for.                              C
C     NDRVM     : Maximum number of derivatives required.            C
C                                                                    C
C     IWORK     : Integer work array. Dimension ( NCFM )             C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FDCMBD( NR, NBN, NLMN, NRMN, NLMC, NRMC, NCFM,
     1                   NFDCM, NDRVS, NDRVM, IWORK, XARR, FDCM,
     2                   COEFM, WORK1, WORK2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NBN, NLMN, NRMN, NCFM, NFDCM, NDRVS, NDRVM, 
     1        IWORK( NCFM ), NLMC, NRMC
      DOUBLE PRECISION XARR( NR ), FDCM( NFDCM, NR, NDRVM ),
     1                 COEFM( NCFM, NCFM ), WORK1( NCFM ),
     2                 WORK2( NCFM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NDER, IRAD, NLCS, NRCS, NLN, NRN,
     1        NNDS, INDS, I, INODE
      DOUBLE PRECISION DZERO, X0
      PARAMETER ( DZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     .
C     . Check the values of integers ...
C     .
      IF ( NDRVS.GT.NDRVM ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NDRVS  = ',NDRVS
         PRINT *,' NDRVM  = ',NDRVM
         STOP
      ENDIF
C     .
      NLCS = NLMN - NLMC
      NRCS = NRMC - NRMN
C     . 
C     . Check that sufficient points are allowed
C     . for the derivatives ...
C     . 
      IF ( (NRMC-NLMC).LT.(NBN+1) ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NBN  = ', NBN
         PRINT *,' NLMC = ', NLMC
         PRINT *,' NRMC = ', NRMC
         PRINT *,' Insufficient nodes for differencing.'
         STOP
      ENDIF
C     . 
      IF ( NLCS.LT.0 ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NLMN = ', NLMN
         PRINT *,' NLMC = ', NLMC
         STOP
      ENDIF
C     . 
      IF ( NRCS.LT.0 ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NRMN = ', NRMN
         PRINT *,' NRMC = ', NRMC
         STOP
      ENDIF
C     .
      I = MIN( NLCS, NRCS) + NBN
C     .
      IF ( NDRVS.GT.I ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' You have requested deriv.s to order ',NDRVS
         PRINT *,' At one node, you have only', I
         PRINT *,' side nodes to use for differencing.'
         STOP
      ENDIF
C     .
      I = 2*NBN + 1
C     .
      IF ( NCFM.LT.I ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NCFM = ', NCFM
         PRINT *,' NBN  = ', NBN
         STOP
      ENDIF
C     .
      IF ( NFDCM.LT.I ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NFDCM = ', NFDCM
         PRINT *,' NBN  = ', NBN 
         STOP
      ENDIF
C     .
      IF ( NLMN.GT.NRMN ) THEN
         PRINT *,' Subroutine FDCMBD.'
         PRINT *,' NLMN = ', NLMN
         PRINT *,' NRMN = ', NRMN
         STOP
      ENDIF
C     .
C     . Whew ... all input parameters seem to be o.k.
C     . Now loop around the requested rows
C     .
      DO IRAD = NLMN, NRMN
C
        DO NDER = 1, NDRVS
          DO I = 1, NFDCM
            FDCM( I, IRAD, NDER ) = DZERO
          ENDDO
        ENDDO
C
C irad is the node for which we want to calculate
C our coefficients
C
        NLCS = IRAD - NLMC
        NRCS = NRMC - IRAD
C
C we wish to calculate NLN (number of left nodes)
C and NRN ( number of right nodes )
C NNDS ( total number of nodes) is then NLN + NRN + 1 ...
C
        NLN = MIN( NBN, NLCS )
        NRN = MIN( NBN, NRCS )
C
        NNDS = NLN + NRN + 1
C
        X0 = XARR( IRAD )
        DO INDS = 1, NNDS
          INODE = IRAD - NLN - 1 + INDS
          WORK1( INDS ) = XARR( INODE )
        ENDDO
C
C Now ready to calculate the coefficients
C
        CALL GFDCFD( X0, WORK1, NNDS, COEFM, NCFM, 
     1               IWORK, WORK2 )
C
C coefm matrix should now contain the coeff.s
C
        DO NDER = 1, NDRVS
          DO INDS = 1, NNDS
            INODE = INDS - NLN + NBN
            FDCM( INODE, IRAD, NDER ) = COEFM( NDER+1, INDS )
          ENDDO
        ENDDO
C
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Integer Index Array Single Component Extract ************
C            -       -     -     -      -         -       ************
C Steve Gibbons Wed Nov  1 08:51:08 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Takes as input the large integer arrays defining spherical         C
C harmonic properties and fills equivalent arrays with indices       C
C corresponding to a single component (1,2,3,4 or 5)                 C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Total number of harmonics in full harm. set.       C
C     ICOMP     : Single component identifier for comparison with    C
C                  MHT array. Must be one of the following:-         C
C                                                                    C
C         ICOMP = 1 for a poloidal velocity harmonic.                C
C         ICOMP = 2 for a toroidal velocity harmonic.                C
C         ICOMP = 3 for a temperature harmonic.                      C
C         ICOMP = 4 for a poloidal magnetic field harmonic.          C
C         ICOMP = 5 for a toroidal magnetic field harmonic.          C
C                                                                    C
C     MHT       : MHT( ih ) contains itype for harmonic 'ih'         C
C                                                                    C
C     MHL       : MHL( ih ) contains degree, l, for harmonic 'ih'    C
C                                                                    C
C     MHM       : MHM( ih ) contains order, m, for harmonic 'ih' if  C
C                  ih has cos m phi dependency and (-m) if ih has    C
C                   sin m phi dependency.                            C
C                                                                    C
C     MHP       : MHP( ih ) contains is, finite difference scheme    C
C                  identifier.                                       C
C                                                                    C
C     NCH       : (Output) number of harmonics in single component   C
C                   set,                                             C
C                                                                    C
C     NCHMAX    : Maximum value for NCH.                             C
C                                                                    C
C     MLC       : Equiv. of MHL for single component set.            C
C     MMC       : Equiv. of MHL for single component set.            C
C     MPC       : Equiv. of MHL for single component set.            C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE IIASCE( NH, ICOMP, MHT, MHL, MHM, MHP, NCH, NCHMAX,
     1                   MLC, MMC, MPC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, ICOMP, MHT( * ), MHL( * ), MHM( * ), MHP( * ), NCH,
     1        NCHMAX, MLC( * ), MMC( * ), MPC( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH
      LOGICAL OK, KO
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      IF ( ICOMP.LT.1 .OR. ICOMP.GT.5 ) THEN
        PRINT *,' Subroutine IIASCE.'
        PRINT *,' ICOMP = ', ICOMP,' : Invalid value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
c
      OK  = .TRUE.
      NCH = 0
      DO IH = 1, NH
        KO = .FALSE.
        IF ( MHT( IH ).EQ.ICOMP ) THEN
          NCH = NCH + 1
          KO  = .TRUE.
        ENDIF
        IF ( NCH.GT.NCHMAX ) OK = .FALSE.
        IF ( KO .AND. OK ) THEN
          MLC( NCH ) = MHL( IH )
          MMC( NCH ) = MHM( IH )
          MPC( NCH ) = MHP( IH )
        ENDIF
      ENDDO
c
      IF ( OK ) RETURN
c
      PRINT *,' Subroutine IIASCE. NCHMAX = ', NCHMAX
      PRINT *,' You require ', NCH,' harmonics for '
      PRINT *,' component ', ICOMP,'. Program aborted.'
      STOP
      END
C*********************************************************************
C*********************************************************************
C subroutine Multiple Component 2 Single Component Vector ************
C            -        -         - -      -         -      ************
C Steve Gibbons Wed Nov  1 11:08:11 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C   IDIR = 1:                                                        C
C Takes a single component vector VECC defined by MLC, MMC, with     C
C NHC harmonics and fills it up, from grid node ILNR to IRNR, with   C
C values from the corresponding harmonics in the vector VEC which    C
C is defined by the arrays INARR, MHT, MHL and MHM.                  C
C                                                                    C
C   IDIR = -1:                                                       C
C Fills in values of VEC from VECC.                                  C
C                                                                    C
C Note that VECC must be arranged with ind = (ih-1)*nr + ir          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Int. parameter array corresponding to vectors.     C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NR      See INDFUN for details        C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     MHT       : MHT( ih ) contains itype for harmonic 'ih'         C
C                  (see ICOMP)                                       C
C                                                                    C
C     MHL       : MHL( ih ) contains degree, l, for harmonic 'ih'    C
C                                                                    C
C     MHM       : MHM( ih ) contains order, m, for harmonic 'ih' if  C
C                  ih has cos m phi dependency and (-m) if ih has    C
C                   sin m phi dependency.                            C
C                                                                    C
C N.B. MHT, MHL, MHM and INARR all correspond to the array, VEC.     C
C                                                                    C
C     ILNR      : Lowest radial node to transfer.                    C
C     IRNR      : Highest radial node to transfer.                   C
C                                                                    C
C     ICOMP     : Number of component (1,2,3,4 or 5) -               C
C                 corresponds to MHT( ih ) and                       C
C                                                                    C
C         ICOMP = 1 for a poloidal velocity harmonic.                C
C         ICOMP = 2 for a toroidal velocity harmonic.                C
C         ICOMP = 3 for a temperature harmonic.                      C
C         ICOMP = 4 for a poloidal magnetic field harmonic.          C
C         ICOMP = 5 for a toroidal magnetic field harmonic.          C
C                                                                    C
C     NHC       : Number of harmonic in vector VECC                  C
C                                                                    C
C     MLC       : Dim ( * ). Equivalent of MHL, for VECC.            C
C                                                                    C
C     MMC       : Dim ( * ). Equivalent of MHM, for VECC.            C
C                                                                    C
C     IDIR      : Set to 1 to put values from VEC into VECC          C
C                 Set to -1 to put values from VECC into VEC         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Dim ( * ). (Multiple component vector).            C
C     VECC      : Dim ( * ). (Single component vector).              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MC2SCV( INARR, MHT, MHL, MHM, ILNR, IRNR, ICOMP,
     1                   NHC, MLC, MMC, VEC, VECC, IDIR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * ), ILNR, IRNR,
     1        ICOMP, NHC, MLC( * ), MMC( * ), IDIR
      DOUBLE PRECISION VEC( * ), VECC( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IND, INDC, INDFUN, NH, NR, IR, IH, IHC,
     1        IBEGIN
      EXTERNAL INDFUN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      NR = INARR( 2 )
      NH = INARR( 3 )
c
      IF ( IDIR.NE.1 .AND. IDIR.NE.-1 ) THEN
        PRINT *,' Subroutine MC2SCV.'
        PRINT *,' IDIR = ', IDIR,' : Invalid value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
c
      IF ( ICOMP.LT.1 .OR. ICOMP.GT.5 ) THEN
        PRINT *,' Subroutine MC2SCV.'
        PRINT *,' ICOMP = ', ICOMP,' : Invalid value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
c
      DO IHC = 1, NHC
        IBEGIN = (IHC-1)*NR 
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.ICOMP   .AND.   MHL( IH ).EQ.MLC( IHC )
     1           .AND.  MHM( IH ).EQ.MMC( IHC )    ) THEN
            DO IR = ILNR, IRNR
              IND  = INDFUN( IR, IH, INARR )
              INDC = IBEGIN + IR
              IF ( IDIR.EQ.1 ) VECC( INDC ) = VEC( IND )
              IF ( IDIR.EQ.-1 ) VEC( IND ) = VECC( INDC )
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Poloidal Velocity Completeness Coefficients Find ********
C            -        -        -            -            -    ********
C Steve Gibbons Thu Nov  9 10:10:59 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Because of the additional boundary conditions on the poloidal      C
C velocity radial function, the values at grid nodes 2 and NR-1 are  C
C not solved for in the solution vector.                             C
C                                                                    C
C This routine extracts from the array SVFDC, two vectors of length  C
C NBN, called PVLC (poloidal velocity left coefficients) and PVRC    C
C (poloidal velocity right coefficients) such that for poloidal      C
C harmonic IH, the vectors can be completed with the following       C
C segments of code:-                                                 C
C                                                                    C
C     IRND = 2                                                       C
C     IND2 = INDFUN( IRND, IH, INARR )                               C
C     TEMP = 0.0d0                                                   C
C     DO I = 1, NBN                                                  C
C       IR   = IRND + I                                              C
C       IND  = INDFUN( IR, IH, INARR )                               C
C       TEMP = TEMP + PVLC( I )*SV( IND )                            C
C     ENDDO                                                          C
C     SV( IND2 ) = TEMP                                              C
C                                                                    C
C     IRND = NR - 1                                                  C
C     IND2 = INDFUN( IRND, IH, INARR )                               C
C     TEMP = 0.0d0                                                   C
C     DO I = 1, NBN                                                  C
C       IR   = IRND - I                                              C
C       IND  = INDFUN( IR, IH, INARR )                               C
C       TEMP = TEMP + PVRC( I )*SV( IND )                            C
C     ENDDO                                                          C
C     SV( IND2 ) = TEMP                                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C                                                                    C
C     NDCS      : Number of distinct differencing coeff.s            C
C                  represented in SVFDC.                             C
C                                                                    C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C     NFDCM     : Leading order of the array SVFDC.                  C
C                 This must be atleast (2*NBN + 1)                   C
C                                                                    C
C     NDRVM     : Maximum number of derivatives allowed.             C
C                                                                    C
C     IS        : Number of difference scheme for poloidal velocity. c
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                                                                    C
C     PVLC      : Dim (NBN). See above.                              C
C     PVRC      : Dim (NBN). See above.                              C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PVCCF( NR, NDCS, NBN, NFDCM, NDRVM, IS, SVFDC,
     1                  PVLC, PVRC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NR, NDCS, NBN, NFDCM, NDRVM, IS
      DOUBLE PRECISION SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                 PVLC( NBN ), PVRC( NBN )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IR, IND
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IR  = 2
      DO I = 1, NBN
        IND       = NBN + 1 + I
        PVLC( I ) = SVFDC( IND, IR, 1, IS )
      ENDDO
C     .
      IR  = NR - 1
      DO I = 1, NBN
        IND       = NBN + 1 - I
        PVRC( I ) = SVFDC( IND, IR, 1, IS )
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Number of PHi Points Find *******************************
C            -         --  -      -    *******************************
C Steve Gibbons Mon Nov  6 11:36:13 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Chooses an adequate number of phi points to perform a fast         C
C Fourier transform. Now a function f( phi ) is such that            C
C                                                                    C
C f( phi ) = f( phi + 2.pi/M0 )                                      C
C                                                                    C
C If f requires a maximum wavenumber MMAX then NPHPF finds a number  C
C of points in PHI ( NPHP ) such that  NPHP > 2*MMAX/M0 and NPHP     C
C is a power of 2.                                                   C
C                                                                    C
C The i^{th} point, phi_{i} = (i-1)*deltap                           C
C                                                                    C
C where deltap = 2*pi/(NPHP*M0)                                      C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     MMAX      : Highest wavenumber, m, in spectral expansion.      C
C     M0        : Integer, which divides MMAX such that              C
C                        f( phi ) = f( phi + 2.pi/M0 )               C
C                                                                    C
C     NPHP      : Number of phi points required.                     C
C     NPHMAX    : Maximum number of phi points allowed.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NPHPF( MMAX, M0, NPHP, NPHMAX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER MMAX, M0, NPHP, NPHMAX
C____________________________________________________________________C
C Variable declarations - working variables .........................C
      INTEGER NPHL
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( M0.LE.0 .OR. MMAX.LT.0 ) THEN
        PRINT *,' Subroutine NPHPF.'
        PRINT *,' MMAX = ', MMAX,', M0 = ', M0
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( MMAX/M0*M0.NE.MMAX ) THEN
        PRINT *,' Subroutine NPHPF.'
        PRINT *,' MMAX = ', MMAX,', M0 = ', M0
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NPHL = 2*MMAX/M0
      NPHP = 2
 500  CONTINUE
      IF ( NPHL.GE.NPHP ) THEN
        NPHP = NPHP*2
        GOTO 500
      ENDIF
C
      IF ( NPHP.GT.NPHMAX ) THEN
        PRINT *,' Subroutine NPHPF.'
        PRINT *,' MPHP   = ', NPHP
        PRINT *,' MPHMAX = ', NPHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine GAUWTS **************************************************
C Adapted 22.4.97 from Numerical Recipes routine GAULEG              C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHPTS	: Number of theta points.                            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X1	: Starting value for integration.                    C
C     X2	: Ending value for integration.                      C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX	: Array containing abscissae of the Gauss-Legendre   C
C                  NTHPTS-points quadrature formula.                 C
C     GAUW      : Array containing the weights for the above points. C
C                                                                    C
C ( Both GAUX and GAUW have dimension NTHPTS ).                      C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE GAUWTS ( X1, X2, GAUX, GAUW, NTHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHPTS
      DOUBLE PRECISION X1, X2, GAUX( NTHPTS ), GAUW( NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER M,I,J
      DOUBLE PRECISION XM,XL,P1,P2,P3,EPS,PP,Z,Z1
      PARAMETER (EPS=1.0d-13)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C 
      IF ( ABS( X1 - X2 ).LT.EPS ) THEN
        PRINT *,' Subroutine GAUWTS,'
        PRINT *,' X1 = ', X1
        PRINT *,' X2 = ', X2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C ....... the roots are symmetric in the interval so need only find
C half of them.
      M = ( NTHPTS + 1 )/2
      XM = 0.5d0 * ( X2 + X1 )
      XL = 0.5d0 * ( X2 - X1 )
C ........start looping over the desired roots .......................
      DO I = 1, M
         Z = DCOS( PI*( DBLE(I) - 0.25d0)/( DBLE(NTHPTS) + 0.5d0 ))
C           ..... starting with this approximation to the Ith root, we
C                enter the main loop of refinement by Newton's method.
 100     CONTINUE
            P1 = 1.0D0
            P2 = 0.0D0
C           ........... Loop up the recurrence relation to get the
C                      legendre Polynomial evaluated at Z.
            DO J = 1, NTHPTS
               P3 = P2
               P2 = P1
               P1 = ((2.0d0*J-1.0d0)*Z*P2 - (J-1.0d0)*P3)/DBLE( J )
            ENDDO
C           ..................... finish recurrence relation loop ...
C ... P1 is now the desired Legendre Polynomial. We now compute PP,
C    its derivative by a standard relation involving also P2, the 
C    polynomial of one order lower.
            PP = NTHPTS*(Z*P1-P2)/(Z*Z-1.0d0)
            Z1 = Z
            Z = Z1 - P1/PP
         IF ( ABS(Z-Z1).GT.EPS ) GOTO 100
C ...........scale the root to the desired interval .................
         GAUX( I ) = XM - XL*Z
C ...........and add its symmetric counterpart ......................
         GAUX( NTHPTS+1-I ) = XM + XL*Z
C ...........calculate the weight ...................................
         GAUW( I ) = 2.0d0*XL/((1.0d0-Z*Z)*PP*PP)
C ...........and add its symmetric counterpart ......................
         GAUW( NTHPTS + 1 - I ) = GAUW( I )
      ENDDO
C ......... end looping over the desired roots .......................
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine SCHmidt Normalised Legendre function Array **************
C            ---     -          -                 -     **************
C Steve Gibbons 22.4.97                                              C
C____________________________________________________________________C
C Does the same as SCHNLF except that instead of a single valued X   C
C for one theta point, it fills arrays PA and DPA with the           C
C legendre Functions etc. for each of the NTHPTS values of cos(theta)C
C in the array GAUX.                                                 C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH	: Highest degree, l, of spherical harmonic.          C
C     NTHPTS	: Number of theta points.                            C
C  Double Precision                                                  C
C  ----------------                                                  C
C     PA	: Schmidt Normalised Legendre Functions Dimension.   C
C		   {  ( LH + 1 )*( LH + 2 )/2 , NTHPTS }             C
C		   P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA	: Derivatives of the above.                          C
C     GAUX	: Array of cosines to the NTHPTS angles.             C
C                  Dimension ( NTHPTS ).                             C
C____________________________________________________________________C
C Functions ... Calling Proceedures :-                               C
C  Double Precision                                                  C
C  ----------------                                                  C
C PMM ( M, S )					                     C
C DPMM ( M, C, S )				                     C
C PMM1 ( M, X, PMM0 )                                                C
C PLM ( L, M, X, PLMIN1, PLMIN2 )				     C
C DPMM1 ( M , X , S, PMM , DPMM)                                     C
C DPLM ( L, M , X , S, PMM1 , DPMM1, DPMM2 )                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SCHNLA ( PA, DPA, GAUX, LH, NTHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS
      DOUBLE PRECISION GAUX( NTHPTS ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER INDEX,L,M,IOLD1,IOLD2,NTHETA
      DOUBLE PRECISION SINE,PMIN1,PMIN2,TOL,DPMIN1,
     1                 DPMIN2,X
      PARAMETER (TOL=1.0d-6)
C____________________________________________________________________C
C Variable declarations - Functions called ..........................C
      DOUBLE PRECISION PMM,PMM1,PLM,DPMM,DPMM1,DPLM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check validity of arguments ....
C
C ........................... now loop around theta points
      DO NTHETA = 1, NTHPTS
C
        X = GAUX ( NTHETA )
        SINE = X*X
        IF ( SINE.GT.1.0d0 ) THEN
           PRINT *,' Subroutine SCHNLA.'
           PRINT *,' Illegal Cos(theta) has been entered.'
           PRINT *,' ( For NTHETA = ',NTHETA,' )'
           PRINT *,' Program stopped.'
           STOP
        ENDIF
C
C Set SINE (theta) in terms of X
        SINE = DSQRT ( (1.0d0 + X)*(1.0d0 - X) )
        IF ( SINE.LT.TOL ) THEN
           PRINT *,' Subroutine SCHNLA.'
           PRINT *,' SINE is too small. Division by zero imminent.'
           PRINT *,' ( For NTHETA = ',NTHETA,' )'
           PRINT *,' Program stopped.'
           STOP
        ENDIF
C..................... first calculate the P_l^m ..........
        DO M = 0, LH - 2
C                        ............. Calculate P_M^M .....
           L = M
           INDEX = L*(L+1)/2+M+1
           PA ( INDEX , NTHETA) = PMM ( M , SINE )
           DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
C                         ............. Calculate P_(M+1)^M .
           PMIN1 = PA ( INDEX , NTHETA)
           DPMIN1 = DPA ( INDEX , NTHETA)
           IOLD2 = INDEX
           L = L + 1
           INDEX = L*(L+1)/2+M+1
           PA (INDEX , NTHETA) = PMM1 ( M , X , PMIN1 )
           DPA (INDEX , NTHETA) = DPMM1 (M,X,SINE,PMIN1, DPMIN1)
           IOLD1 = INDEX
C                         ......... Calculate P_L^M general .
           DO L = M + 2, LH
              PMIN2 = PA ( IOLD2 , NTHETA)
              PMIN1 = PA ( IOLD1 , NTHETA)
              DPMIN2 = DPA ( IOLD2 , NTHETA)
              DPMIN1 = DPA ( IOLD1 , NTHETA)
              INDEX = L*(L+1)/2+M+1
              PA ( INDEX , NTHETA) = PLM ( L,M,X,PMIN1,PMIN2 )
              DPA ( INDEX , NTHETA) = DPLM (L,M,X,SINE , PMIN1, 
     1                              DPMIN1, DPMIN2 )
              IOLD2 = IOLD1
              IOLD1 = INDEX
           ENDDO
        ENDDO
        M = LH - 1
        L = M
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM ( M , SINE )
        DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
        PMIN1 = PA( INDEX , NTHETA)
        DPMIN1 = DPA( INDEX , NTHETA)
        L = LH
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM1 ( M , X , PMIN1 )
        DPA(INDEX ,NTHETA) = DPMM1 (M ,X ,SINE,PMIN1,DPMIN1 )
        M = LH
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM ( M , SINE )
        DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
C......................finished calculating P_l^m .........
      ENDDO
C......................finished looping around theta points

      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Poloidal Velocity Time-Step Matrices Form ***************
C            -        -        -    -    -        -    ***************
C Steve Gibbons Wed Nov  1 15:04:22 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we are time-stepping, solving for a poloidal velocity, p^{i+1}, C
C such that                                                          C
C                                                                    C
C           [ p^{i+1} - p^i ]            [ (1-c) \nabla^2 p^{i+1} ]  C
C  c_e curl [ ------------- ] = c_i curl [    +  c \nabla^2 p^i   ]  C
C           [  delta t      ]    + Additional_forcing_terms          C
C                                                                    C
C then PVTSMF builds the two matrices, AM1 and BM1, such that        C
C                                                                    C
C AM1 p^{i+1} = BM1 p^{i} + Forcing terms.                           C
C                                                                    C
C The solution vector must consist only of poloidal velocity         C
C harmonics, of which there are NH1, and the values of these         C
C radial functions, at NR grid nodes, must be given by the index     C
C                                                                    C
C   ind = ( ih - 1 )*nr + ir                                         C
C                                                                    C
C where ih is the number of the radial function and ir is the        C
C number of the radial grid node.                                    C
C                                                                    C
C The radial functions are characterised by the integer arrays       C
C ML1, MM1 (not referred to here) and MP1.                           C
C                                                                    C
C ML1( ih ) gives the spherical harmonic degree, l.                  C
C MM1( ih ) gives the spherical harmonic order, m, when ih has a     C
C (cos m phi) dependence and -m for a (sin m phi) dependence.        C
C                                                                    C
C NBN is the number of bounding nodes.                               C
C Rows are added to the matrices between grid nodes 3 and NR - 2.    C
C The matrix AM1 is LU decomposed and so additional diagonal         C
C elements must be added using the routine AMSDEA. Due to the        C
C LU decomposition, AM1 must have the dimensions                     C
C  ( 3*NBN+1, NH1*NR )                                               C
C whereas BM1 is only required for a matrix-vector multiplication    C
C and so only the dimensions ( 2*NBN+1, NH1*NR ) are needed.         C
C                                                                    C
C It zeroes both matrices on input and assumes a LAPACK format for   C
C the banded matrix ( IMF = 1 in MATIND ).                           C
C                                                                    C
C The boundary conditions imposed in the different finite diff.      C
C schemes are given by the arrays MHIBC and MHOBC. Now MP1(ih) = is  C
C where MHIBC( is ) is one of the following:-                        C
C                                                                    C
C  MHIBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C                                                                    C
C Similarly for MHOBC. MHIBC and MHOBC need not be identical.        C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NH1       : Number of poloidal velocity harmonics.             C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
C                                                                    C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). See above.                     C
C     MHOBC     : Dimension ( NDCS ). See above.                     C
C     NFDCM     : Leading order of the array SVFDC.                  C
C                 This must be atleast (2*NBN + 1)                   C
C     NDRVM     : Maximum number of derivatives in SVFDC.            C
C                                                                    C
C     ML1       : Dim ( NH1 ). See above.                            C
C     MP1       : Dim ( NH1 ). See above.                            C
C                                                                    C
C     IPIV1     : Dim (NH1*NR). Pivotting information for            C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                                                                    C
C     AM1      : Matrix for solution of p^{i+1}. Has dimensions      C
C                 ( 3*NBN + 1, NH1*NR )                              C
C                                                                    C
C     BM1      : Matrix for multiplication of p^i. Has dimensions    C
C                 ( 2*NBN + 1, NH1*NR )                              C
C                                                                    C
C     CE        : Coefficient of time-derivative.                    C
C     CI        : Coefficient of viscous diffusion.                  C
C     CFAC      : Determines how implicit/explicit integration is.   C
C                 Must be strictly greater than 0.0 and strictly     C
C                 less than 1.0  The higher CFAC is, the more        C
C                 explicit the integration. CFAC = 0.5 --> Crank-    C
C                 Nicolson scheme.                                   C
C     DELTAT    : Time-step size.                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PVTSMF( NR, NH1, NDCS, NBN, MHIBC, MHOBC, NFDCM,
     1                   NDRVM, ML1, MP1, IPIV1, XARR, SVFDC,
     2                   AM1, BM1, CE, CI, CFAC, DELTAT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH1, NDCS, NBN, MHIBC( NDCS ), MHOBC( NDCS ),
     1        NFDCM, NDRVM, ML1( * ), MP1( * ), IPIV1( * )
      DOUBLE PRECISION XARR( * ), SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                AM1( 3*NBN+1, NH1*NR ), BM1( 2*NBN+1, NH1*NR ),
     2                CE, CI, CFAC, DELTAT 
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, IMF, INARR( 3 ), ILNR, IRNR, INFO, N1, N2, KL,
     1        IH1, IHD, IPARS( 2 ), IS, KLE
      DOUBLE PRECISION ZERO, ONE, FAC, DPARS( 1 ), WORK( 5 )
      PARAMETER ( ZERO = 0.0d0, ONE = 1.0d0, IOP = 0 )
      EXTERNAL AMDLT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check input parameters
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH1
      DPARS( 1 ) = ZERO
C (actually dpars is not referred to by AMDLT)
C
C Enforce matrix is built in LAPACK format
C
      IMF = 1
      KL  = NBN
C
C Check CFAC
C
      IF ( CFAC.LE.ZERO .OR. CFAC.GE.ONE ) THEN
        PRINT *,' Subroutine PVTSMF.'
        PRINT *,' CFAC = ', CFAC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DELTAT.EQ.ZERO .OR. CE.EQ.ZERO .OR. CI.EQ.ZERO ) THEN
        PRINT *,' Subroutine PVTSMF.'
        PRINT *,' CE = ', CE,' CI = ', CI,' DELTAT = ', DELTAT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Zero the matrix AM1
C
      N2 = NH1*NR
      N1 = 3*NBN + 1
      CALL MATOP( AM1, ZERO, N1, N2, IOP )
C
      KLE  = KL
      ILNR = 3
      IRNR = NR - 2
C
C Add the curl of \nabla^2 parts onto AM1 matrix
C Curl of a poloidal harmonic is the -D_l^2 operator
C
      FAC  = CI*DELTAT*(1.0d0 - CFAC)
C
      IHD  = 4
      DO IH1 = 1, NH1
        IS         = MP1( IH1 )
        IPARS( 1 ) = 3
        IPARS( 2 ) = ML1( IH1 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH1, IH1, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM1, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the curl parts of time derivative onto AM1 matrix
C Curl is -D_l for poloidal harmonics and so we multiply
C CE by -1.0d0 to make FAC
C
      FAC  = (-1.0d0)*CE
      IHD  = 2
      DO IH1 = 1, NH1
        IS         = MP1( IH1 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML1( IH1 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH1, IH1, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM1, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the diagonal elements to AM1 matrix.
C
      FAC = 1.0d0
      CALL AMSDEA( AM1, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP1, MHIBC, 'Inner Boundary', FAC, NDCS )
      CALL AMSDEA( AM1, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP1, MHOBC, 'Outer Boundary', FAC, NDCS )
C
C Now we must attempt to do an LU decomposition of the
C DMAT matrix ... for this, we need the LAPACK routine DGBTRF
C
      CALL DGBTRF( N2, N2, KL, KL, AM1, N1, IPIV1, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine PVTSMF.'
        PRINT *,' The LAPACK subroutine DGBTRF has been called'
        PRINT *,' and has returned ',INFO,' as a value of INFO. '
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now zero the matrix BM1 ...
C
      N1 = 2*NBN + 1
      CALL MATOP( BM1, ZERO, N1, N2, IOP )
      KLE = 0
C
C Add the curl of \nabla^2 parts onto BM1 matrix.
C Curl of a poloidal harmonic is the -D_l^2 operator
C
      FAC  = CI*DELTAT*CFAC*(-1.0d0)
      IHD  = 4
      DO IH1 = 1, NH1
        IS         = MP1( IH1 )
        IPARS( 1 ) = 3
        IPARS( 2 ) = ML1( IH1 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH1, IH1, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM1, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the curl parts of time derivative onto BM1 matrix
C Curl is -D_l for poloidal harmonics and so we multiply
C CE by -1.0d0 to make FAC
C
      FAC  = (-1.0d0)*CE
      IHD  = 2
      DO IH1 = 1, NH1
        IS         = MP1( IH1 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML1( IH1 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH1, IH1, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM1, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Toroidal Velocity Time-Step Matrices Form ***************
C            -        -        -    -    -        -    ***************
C Steve Gibbons Fri Nov  3 09:22:18 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we are time-stepping, solving for a toroidal velocity, t^{i+1}, C
C such that                                                          C
C                                                                    C
C           [ t^{i+1} - t^i ]            [ (1-c) \nabla^2 t^{i+1} ]  C
C  c_e curl [ ------------- ] = c_i curl [    +  c \nabla^2 t^i   ]  C
C           [  delta t      ]    + Additional_forcing_terms          C
C                                                                    C
C then TVTSMF builds the two matrices, AM2 and BM2, such that        C
C                                                                    C
C AM2 t^{i+1} = BM2 t^{i} + Forcing terms.                           C
C                                                                    C
C The solution vector must consist only of toroidal velocity         C
C harmonics, of which there are NH2, and the values of these         C
C radial functions, at NR grid nodes, must be given by the index     C
C                                                                    C
C   ind = ( ih - 1 )*nr + ir                                         C
C                                                                    C
C where ih is the number of the radial function and ir is the        C
C number of the radial grid node.                                    C
C                                                                    C
C The radial functions are characterised by the integer arrays       C
C ML2, MM2 (not referred to here) and MP2.                           C
C                                                                    C
C ML2( ih ) gives the spherical harmonic degree, l.                  C
C MM2( ih ) gives the spherical harmonic order, m, when ih has a     C
C (cos m phi) dependence and -m for a (sin m phi) dependence.        C
C                                                                    C
C NBN is the number of bounding nodes.                               C
C Rows are added to the matrices between grid nodes 2 and NR - 1.    C
C The matrix AM2 is LU decomposed and so additional diagonal         C
C elements must be added using the routine AMSDEA. Due to the        C
C LU decomposition, AM2 must have the dimensions                     C
C  ( 3*NBN+1, NH2*NR )                                               C
C whereas BM2 is only required for a matrix-vector multiplication    C
C and so only the dimensions ( 2*NBN+1, NH2*NR ) are needed.         C
C                                                                    C
C It zeroes both matrices on input and assumes a LAPACK format for   C
C the banded matrix ( IMF = 1 in MATIND ).                           C
C                                                                    C
C The boundary conditions imposed in the different finite diff.      C
C schemes are given by the arrays MHIBC and MHOBC. Now MP2(ih) = is  C
C where MHIBC( is ) is one of the following:-                        C
C                                                                    C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C                                                                    C
C Similarly for MHOBC. MHIBC and MHOBC need not be identical.        C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NH2       : Number of toroidal velocity harmonics.             C
C     NDCS      : Number of distinct differencing coeff.s            C
C                 represented in SVFDC.                              C
C                                                                    C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). See above.                     C
C     MHOBC     : Dimension ( NDCS ). See above.                     C
C     NFDCM     : Leading order of the array SVFDC.                  C
C                 This must be atleast (2*NBN + 1)                   C
C     NDRVM     : Maximum number of derivatives in SVFDC.            C
C                                                                    C
C     ML2       : Dim ( NH2 ). See above.                            C
C     MP2       : Dim ( NH2 ). See above.                            C
C                                                                    C
C     IPIV2     : Dim (NH2*NR). Pivotting information for            C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                                                                    C
C     AM2      : Matrix for solution of p^{i+1}. Has dimensions      C
C                 ( 3*NBN + 1, NH2*NR )                              C
C                                                                    C
C     BM2      : Matrix for multiplication of p^i. Has dimensions    C
C                 ( 2*NBN + 1, NH2*NR )                              C
C                                                                    C
C     CE        : Coefficient of time-derivative.                    C
C     CI        : Coefficient of viscous diffusion.                  C
C     CFAC      : Determines how implicit/explicit integration is.   C
C                 Must be strictly greater than 0.0 and strictly     C
C                 less than 1.0  The higher CFAC is, the more        C
C                 explicit the integration. CFAC = 0.5 --> Crank-    C
C                 Nicolson scheme.                                   C
C     DELTAT    : Time-step size.                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE TVTSMF( NR, NH2, NDCS, NBN, MHIBC, MHOBC, NFDCM,
     1                   NDRVM, ML2, MP2, IPIV2, XARR, SVFDC,
     2                   AM2, BM2, CE, CI, CFAC, DELTAT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH2, NDCS, NBN, MHIBC( NDCS ), MHOBC( NDCS ),
     1        NFDCM, NDRVM, ML2( * ), MP2( * ), IPIV2( * )
      DOUBLE PRECISION XARR( * ), SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                AM2( 3*NBN+1, NH2*NR ), BM2( 2*NBN+1, NH2*NR ),
     2                CE, CI, CFAC, DELTAT 
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, IMF, INARR( 3 ), ILNR, IRNR, INFO, N1, N2, KL,
     1        IH2, IHD, IPARS( 2 ), IS, KLE
      DOUBLE PRECISION ZERO, ONE, FAC, DPARS( 1 ), WORK( 5 )
      PARAMETER ( ZERO = 0.0d0, ONE = 1.0d0, IOP = 0 )
      EXTERNAL AMDLT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check input parameters
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH2
      DPARS( 1 ) = ZERO
C (actually dpars is not referred to by AMDLT)
C
C Enforce matrix is built in LAPACK format
C
      IMF = 1
      KL  = NBN
C
C Check CFAC
C
      IF ( CFAC.LE.ZERO .OR. CFAC.GE.ONE ) THEN
        PRINT *,' Subroutine TVTSMF.'
        PRINT *,' CFAC = ', CFAC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DELTAT.EQ.ZERO .OR. CE.EQ.ZERO .OR. CI.EQ.ZERO ) THEN
        PRINT *,' Subroutine TVTSMF.'
        PRINT *,' CE = ', CE,' CI = ', CI,' DELTAT = ', DELTAT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Zero the matrix AM2
C
      N2 = NH2*NR
      N1 = 3*NBN + 1
      CALL MATOP( AM2, ZERO, N1, N2, IOP )
C
      KLE  = KL
      ILNR = 2
      IRNR = NR - 1
C
C Add the curl of \nabla^2 parts onto AM2 matrix
C Curl Lap of a toroidal vel. harmonic is the D_l operator
C
      FAC  = CI*DELTAT*(CFAC - 1.0d0)
C
      IHD  = 2
      DO IH2 = 1, NH2
        IS         = MP2( IH2 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML2( IH2 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH2, IH2, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM2, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the curl parts of time derivative onto AM2 matrix
C Curl is D_l^0(!) for toroidal harmonics and so we multiply
C CE by +1.0d0 to make FAC
C
      FAC  = CE
      IHD  = 0
      DO IH2 = 1, NH2
        IS         = MP2( IH2 )
        IPARS( 1 ) = 1
        IPARS( 2 ) = ML2( IH2 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH2, IH2, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM2, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the diagonal elements to AM2 matrix.
C
      FAC = 1.0d0
      CALL AMSDEA( AM2, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP2, MHIBC, 'Inner Boundary', FAC, NDCS )
      CALL AMSDEA( AM2, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP2, MHOBC, 'Outer Boundary', FAC, NDCS )
C
C Now we must attempt to do an LU decomposition of the
C AM2 matrix ... for this, we need the LAPACK routine DGBTRF
C
      CALL DGBTRF( N2, N2, KL, KL, AM2, N1, IPIV2, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine TVTSMF.'
        PRINT *,' The LAPACK subroutine DGBTRF has been called'
        PRINT *,' and has returned ',INFO,' as a value of INFO. '
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now zero the matrix BM2 ...
C
      N1 = 2*NBN + 1
      CALL MATOP( BM2, ZERO, N1, N2, IOP )
      KLE = 0
C
C Add the curl of \nabla^2 parts onto BM2 matrix.
C Curl Lap of a toroidal harmonic is the D_l operator
C
      FAC  = CI*DELTAT*CFAC
      IHD  = 2
      DO IH2 = 1, NH2
        IS         = MP2( IH2 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML2( IH2 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH2, IH2, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM2, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the curl parts of time derivative onto BM2 matrix
C Curl is D_l^0(!) for toroidal harmonics and so we multiply
C CE by +1.0d0 to make FAC
C
      FAC  = CE
      IHD  = 0
      DO IH2 = 1, NH2
        IS         = MP2( IH2 )
        IPARS( 1 ) = 1
        IPARS( 2 ) = ML2( IH2 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH2, IH2, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM2, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine TeMperature Time-Step Matrices Form *********************
C            - -         -    -    -        -    *********************
C Steve Gibbons Fri Nov  3 11:00:51 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we are time-stepping, solving for a temperature, T^{i+1},       C
C such that                                                          C
C                                                                    C
C       [ T^{i+1} - T^i ]         [ (1-c) \nabla^2 T^{i+1} ]         C
C  c_a  [ ------------- ] =  c_d  [    +  c \nabla^2 T^i   ]         C
C       [    delta t    ]      +    Additional_forcing_terms         C
C                                                                    C
C then TMTSMF builds the two matrices, AM3 and BM3, such that        C
C                                                                    C
C AM3 t^{i+1} = BM3 t^{i} + Forcing terms.                           C
C                                                                    C
C The solution vector must consist only of temperature               C
C harmonics, of which there are NH3, and the values of these         C
C radial functions, at NR grid nodes, must be given by the index     C
C                                                                    C
C   ind = ( ih - 1 )*nr + ir                                         C
C                                                                    C
C where ih is the number of the radial function and ir is the        C
C number of the radial grid node.                                    C
C                                                                    C
C The radial functions are characterised by the integer arrays       C
C ML3, MM3 (not referred to here) and MP3.                           C
C                                                                    C
C ML3( ih ) gives the spherical harmonic degree, l.                  C
C MM3( ih ) gives the spherical harmonic order, m, when ih has a     C
C (cos m phi) dependence and -m for a (sin m phi) dependence.        C
C                                                                    C
C NBN is the number of bounding nodes.                               C
C Rows are added to the matrices between grid nodes 2 and NR - 1.    C
C The matrix AM3 is LU decomposed and so additional diagonal         C
C elements must be added using the routine AMSDEA. Due to the        C
C LU decomposition, AM3 must have the dimensions                     C
C  ( 3*NBN+1, NH3*NR )                                               C
C whereas BM3 is only required for a matrix-vector multiplication    C
C and so only the dimensions ( 2*NBN+1, NH3*NR ) are needed.         C
C                                                                    C
C It zeroes both matrices on input and assumes a LAPACK format for   C
C the banded matrix ( IMF = 1 in MATIND ).                           C
C                                                                    C
C The boundary conditions imposed in the different finite diff.      C
C schemes are given by the arrays MHIBC and MHOBC. Now MP3(ih) = is  C
C where MHIBC( is ) is one of the following:-                        C
C                                                                    C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C                                                                    C
C Similarly for MHOBC. MHIBC and MHOBC need not be identical.        C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NH3       : Number of toroidal velocity harmonics.             C
C     NDCS      : Number of distinct differencing coeff.s            C
C                 represented in SVFDC.                              C
C                                                                    C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). See above.                     C
C     MHOBC     : Dimension ( NDCS ). See above.                     C
C     NFDCM     : Leading order of the array SVFDC.                  C
C                 This must be atleast (2*NBN + 1)                   C
C     NDRVM     : Maximum number of derivatives in SVFDC.            C
C                                                                    C
C     ML3       : Dim ( NH3 ). See above.                            C
C     MP3       : Dim ( NH3 ). See above.                            C
C                                                                    C
C     IPIV3     : Dim (NH3*NR). Pivotting information for            C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                                                                    C
C     AM3      : Matrix for solution of p^{i+1}. Has dimensions      C
C                 ( 3*NBN + 1, NH3*NR )                              C
C                                                                    C
C     BM3      : Matrix for multiplication of p^i. Has dimensions    C
C                 ( 2*NBN + 1, NH3*NR )                              C
C                                                                    C
C     CA        : Coefficient of time-derivative.                    C
C     CD        : Coefficient of thermal diffusion.                  C
C     CFAC      : Determines how implicit/explicit integration is.   C
C                 Must be strictly greater than 0.0 and strictly     C
C                 less than 1.0  The higher CFAC is, the more        C
C                 explicit the integration. CFAC = 0.5 --> Crank-    C
C                 Nicolson scheme.                                   C
C     DELTAT    : Time-step size.                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE TMTSMF( NR, NH3, NDCS, NBN, MHIBC, MHOBC, NFDCM,
     1                   NDRVM, ML3, MP3, IPIV3, XARR, SVFDC,
     2                   AM3, BM3, CA, CD, CFAC, DELTAT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH3, NDCS, NBN, MHIBC( NDCS ), MHOBC( NDCS ),
     1        NFDCM, NDRVM, ML3( * ), MP3( * ), IPIV3( * )
      DOUBLE PRECISION XARR( * ), SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                AM3( 3*NBN+1, NH3*NR ), BM3( 2*NBN+1, NH3*NR ),
     2                CA, CD, CFAC, DELTAT 
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, IMF, INARR( 3 ), ILNR, IRNR, INFO, N1, N2, KL,
     1        IH3, IHD, IPARS( 2 ), IS, KLE
      DOUBLE PRECISION ZERO, ONE, FAC, DPARS( 1 ), WORK( 5 )
      PARAMETER ( ZERO = 0.0d0, ONE = 1.0d0, IOP = 0 )
      EXTERNAL AMDLT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check input parameters
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH3
      DPARS( 1 ) = ZERO
C (actually dpars is not referred to by AMDLT)
C
C Enforce matrix is built in LAPACK format
C
      IMF = 1
      KL  = NBN
C
C Check CFAC
C
      IF ( CFAC.LE.ZERO .OR. CFAC.GE.ONE ) THEN
        PRINT *,' Subroutine TMTSMF.'
        PRINT *,' CFAC = ', CFAC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DELTAT.EQ.ZERO .OR. CA.EQ.ZERO .OR. CD.EQ.ZERO ) THEN
        PRINT *,' Subroutine TMTSMF.'
        PRINT *,' CA = ', CA,' CD = ', CD,' DELTAT = ', DELTAT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Zero the matrix AM3
C
      N2 = NH3*NR
      N1 = 3*NBN + 1
      CALL MATOP( AM3, ZERO, N1, N2, IOP )
C
      KLE  = KL
      ILNR = 2
      IRNR = NR - 1
C
C Add \nabla^2 parts onto AM3 matrix
C Lap of a temperature harmonic is the D_l operator
C
      FAC  = CD*DELTAT*(CFAC - 1.0d0)
C
      IHD  = 2
      DO IH3 = 1, NH3
        IS         = MP3( IH3 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML3( IH3 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH3, IH3, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM3, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the time derivative onto AM3 matrix
C
      FAC  = CA
      IHD  = 0
      DO IH3 = 1, NH3
        IS         = MP3( IH3 )
        IPARS( 1 ) = 1
        IPARS( 2 ) = ML3( IH3 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH3, IH3, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM3, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the diagonal elements to AM3 matrix.
C
      FAC = 1.0d0
      CALL AMSDEA( AM3, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP3, MHIBC, 'Inner Boundary', FAC, NDCS )
      CALL AMSDEA( AM3, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP3, MHOBC, 'Outer Boundary', FAC, NDCS )
C
C Now we must attempt to do an LU decomposition of the
C AM3 matrix ... for this, we need the LAPACK routine DGBTRF
C
      CALL DGBTRF( N2, N2, KL, KL, AM3, N1, IPIV3, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine TMTSMF.'
        PRINT *,' The LAPACK subroutine DGBTRF has been called'
        PRINT *,' and has returned ',INFO,' as a value of INFO. '
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now zero the matrix BM3 ...
C
      N1 = 2*NBN + 1
      CALL MATOP( BM3, ZERO, N1, N2, IOP )
      KLE = 0
C
C Add the \nabla^2 parts onto BM3 matrix.
C Lap of a temperature harmonic is the D_l operator
C
      FAC  = CD*DELTAT*CFAC
      IHD  = 2
      DO IH3 = 1, NH3
        IS         = MP3( IH3 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML3( IH3 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH3, IH3, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM3, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the time derivative onto BM3 matrix
C
      FAC  = CA
      IHD  = 0
      DO IH3 = 1, NH3
        IS         = MP3( IH3 )
        IPARS( 1 ) = 1
        IPARS( 2 ) = ML3( IH3 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH3, IH3, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM3, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Adapted Vector operation Banded Matrix Building Routine *
C            -       -                       -      -        -       *
C Steve Gibbons Tue Nov 28 14:54:39 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Forms banded matrices to perform operations on a vector of         C
C spherical harmonic radial functions. The operation required is     C
C specified by the integer flag IOPT:                                C
C                                                                    C
C Possible values for IOPT are                                       C
C ----------------------------                                       C
C                                                                    C
C  1: multiply p(r) to get Q(r)                                      C
C                                                                    C
C        now Q( r ) = L(L+1)/RAD p( r )                              C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  2: multiply p(r) to get S(r)                                      C
C                                                                    C
C        now S( r ) = DSQRT( L(L+1) )( p/RAD + dp/dr )               C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C  3: multiply tau(r) to get T(r)                                    C
C                                                                    C
C        now T( r ) = -DSQRT( L(L+1) ) tau( r )                      C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C In the following q, s and t refer to the vectors                   C
C                                                                    C
C  4: curl of scaloidal vector Q(r) q                                C
C                                                                    C
C       curl[ Q(r) q ] = -DSQRT( L(L+1) ) Q( r )/RAD  t              C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  5: curl of spheroidal vector S(r) s                               C
C                                                                    C
C       curl[ S(r) s ] = t[ dS/dr + S(r)/RAD ]                       C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C  6: curl of toroidal vector T(r) t: (scaloidal component)          C
C                                                                    C
C       curl[ T(r) t ] = - q[ T(r) DSQRT( L(L+1) )/RAD ]             C
C                        - s[ dT/dr + T(r)/RAD ]                     C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  7: curl of toroidal vector T(r) t: (spheroidal component)         C
C                                                                    C
C       curl[ T(r) t ] = - q[ T(r) DSQRT( L(L+1) )/RAD ]             C
C                        - s[ dT/dr + T(r)/RAD ]                     C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C  8: Mulitply Q(r) to get p(r)                                      C
C                                                                    C
C       p( r ) = RAD/( L*L + L ) Q( r )                              C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  9: Mulitply T(r) to get tau(r)                                    C
C                                                                    C
C       tau( r ) = - T( r ) / SQRT( L*L + L )                        C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C 10: Calculate a pure first derivative.                             C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C 11: Curl of scaloidal vector Q(r) q to give tau( r ) radial        C
C      function                                                      C
C                                                                    C
C        curl[ Q(r) q ] = T_c( r ) t                                 C
C                                                                    C
C     where T_c( r ) = -DSQRT( L(L+1) ) Q( r )/RAD                   C
C                                                                    C
C     and so tau( r ) = Q( r )/RAD                                   C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C 12: Curl of spheroidal vector S(r) s to give tau( r ) radial       C
C      function                                                      C
C                                                                    C
C        curl[ S(r) s ] = T_c( r ) t                                 C
C                                                                    C
C     where T_c( r ) = dS/dr + S(r)/RAD                              C
C                                                                    C
C       and so tau( r ) =  -( dS/dr + S(r)/RAD )/DSQRT( L(L+1) )     C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C 13: Curl of toroidal vector T( r ) t to give p( r )                C
C                                                                    C
C        curl[ T(r) t ] = Q_c( r ) q   +   S_c( r ) s                C
C                                                                    C
C         where  Q_c( r ) = -T(r) DSQRT( L(L+1) )/RAD   and          C
C                S_c( r ) = -[ dT/dr + T(r)/RAD ]                    C
C                                                                    C
C   hence p( r ) = -T(r)/DSQRT( L(L+1) )                             C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C Unlike AVBMBR, this routine uses the SVFDC finite derivative       C
C array and so will in general assume a particular kind of boundary  C
C condition, or will require values from points 1 to NR in order to  C
C take derivatives.                                                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C The matrix DPOM has dimensions ( N1, N2 ). The vector which is     C
C multiplied by DPOM contains NH radial functions, each with NR      C
C radial grid nodes. The radius at grid node IR is XARR( IR ).       C
C The IRth value of function IH is stored in element                 C
C                                                                    C
C   ( IH - 1 )*NR + IR                                               C
C                                                                    C
C The spherical harmonic degree of the radial function is stored     C
C in MLT( IH ). The finite difference scheme required for this       C
C radial function is given by MPT( IH )                              C
C                                                                    C
C N2 must be equal to NR*NH                                          C
C N1 must be equal to 1 + 2*K where K is the number of diagonals.    C
C                                                                    C
C If IOPT is such that IHD = 0, then K = 0.                          C
C If IOPT is such that IHD = 1, then K = NBN.                        C
C                                                                    C
C The matrix is zero-ed on entry.                                    C
C The terms are added to the matrix between nodes ILNR and IRNR.     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N1        : Leading dimension of DPOM.                         C
C     N2        : Second dimension of DPOM.                          C
C     K         : Number of diagonals in banded matrix.              C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of radial functions.                        C
C     NBN       : Number of side nodes for radial derivatives.       C
C     NCFM      : Leading coefficient of SVFDC.                      C
C     NDRVM     : Maximum number of derivatives stored in SVFDC.     C
C     ILNR      : First grid node to operate upon.                   C
C     IRNR      : Last grid node to operate upon.                    C
C     MLT       : Dim (NH). Element IH contains degree L.            C
C     MPT       : Dim (NH). Element IH contains degree IS.           C
C     IOPT      : Choice of matrix operation. See above list.        C
C     NDCS      : Number of finite difference schemes.               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim (NR). Radial node values.                      C
C     DPOM      : Dim (N1, N2). Double precision operations matrix.  C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AVBMBR( N1, N2, K, NR, NH, NBN, NCFM, NDRVM, ILNR,
     1             IRNR, MLT, MPT, IOPT, NDCS, XARR, DPOM, SVFDC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          N1, N2, K, NR, NH, NBN, NCFM, NDRVM, ILNR,
     1                 IRNR, MLT( NH ), MPT( NH ), IOPT, NDCS
      DOUBLE PRECISION XARR( NR ), DPOM( N1, N2 ),
     1                 SVFDC( NCFM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IARR( 13 ), IHD, IOP, IPARS( 2 ), IMF,
     1                 INARR( 3 ), IH, IS
      DOUBLE PRECISION CVEC( 2 ), ZERO, DPARS( 1 ), FAC
      PARAMETER        ( IOP = 0, ZERO = 0.0d0 )
      EXTERNAL         VOBMAR
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IARR(  1 ) = 0
      IARR(  2 ) = 1
      IARR(  3 ) = 0
      IARR(  4 ) = 0
      IARR(  5 ) = 1
      IARR(  6 ) = 0
      IARR(  7 ) = 1
      IARR(  8 ) = 0
      IARR(  9 ) = 0
      IARR( 10 ) = 1
      IARR( 11 ) = 0
      IARR( 12 ) = 1
      IARR( 13 ) = 0
C
      IHD = IARR( IOPT )
C
      IF ( N1.NE.(2*K+1) .OR. N2.NE.(NR*NH) ) THEN
        PRINT *,' Subroutine AVBMBR.'
        PRINT *,' N1 = ', N1,' N2 = ', N2
        PRINT *,' NR = ', NR,' NH = ', NH
        PRINT *,' K  = ', K
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Clear the matrix
C
      CALL MATOP( DPOM, ZERO, N1, N2, IOP )
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH
      FAC        = 1.0d0
      IMF        = 1
      IPARS( 1 ) = IOPT
      DO IH = 1, NH
        IPARS( 2 ) = MLT( IH )
        IS         = MPT( IH )
        CALL AMLICA( N1, N2, K, K, IOP, IMF, IH, IH, INARR,
     1               IHD, NBN, ILNR, IRNR, NCFM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, VOBMAR, DPOM, FAC,
     3               XARR, CVEC, DPARS, SVFDC )
      ENDDO
C
      END
C*********************************************************************
C*********************************************************************
C subroutine Vector Operations Banded Matrix Building Routine ********
C            -      -          -      -      -        -       ********
C Steve Gibbons Tue Nov 28 13:43:55 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Forms banded matrices to perform operations on a vector of         C
C spherical harmonic radial functions. The operation required is     C
C specified by the integer flag IOPT:                                C
C                                                                    C
C Possible values for IOPT are                                       C
C ----------------------------                                       C
C                                                                    C
C  1: multiply p(r) to get Q(r)                                      C
C                                                                    C
C        now Q( r ) = L(L+1)/RAD p( r )                              C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  2: multiply p(r) to get S(r)                                      C
C                                                                    C
C        now S( r ) = DSQRT( L(L+1) )( p/RAD + dp/dr )               C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C  3: multiply tau(r) to get T(r)                                    C
C                                                                    C
C        now T( r ) = -DSQRT( L(L+1) ) tau( r )                      C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C In the following q, s and t refer to the vectors                   C
C                                                                    C
C  4: curl of scaloidal vector Q(r) q                                C
C                                                                    C
C       curl[ Q(r) q ] = -DSQRT( L(L+1) ) Q( r )/RAD  t              C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  5: curl of spheroidal vector S(r) s                               C
C                                                                    C
C       curl[ S(r) s ] = t[ dS/dr + S(r)/RAD ]                       C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C  6: curl of toroidal vector T(r) t: (scaloidal component)          C
C                                                                    C
C       curl[ T(r) t ] = - q[ T(r) DSQRT( L(L+1) )/RAD ]             C
C                        - s[ dT/dr + T(r)/RAD ]                     C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  7: curl of toroidal vector T(r) t: (spheroidal component)         C
C                                                                    C
C       curl[ T(r) t ] = - q[ T(r) DSQRT( L(L+1) )/RAD ]             C
C                        - s[ dT/dr + T(r)/RAD ]                     C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C  8: Mulitply Q(r) to get p(r)                                      C
C                                                                    C
C       p( r ) = RAD/( L*L + L ) Q( r )                              C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  9: Mulitply T(r) to get tau(r)                                    C
C                                                                    C
C       tau( r ) = - T( r ) / SQRT( L*L + L )                        C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C 10: Calculate a pure first derivative.                             C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C 11: Curl of scaloidal vector Q(r) q to give tau( r ) radial        C
C      function                                                      C
C                                                                    C
C        curl[ Q(r) q ] = T_c( r ) t                                 C
C                                                                    C
C     where T_c( r ) = -DSQRT( L(L+1) ) Q( r )/RAD                   C
C                                                                    C
C     and so tau( r ) = Q( r )/RAD                                   C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C 12: Curl of spheroidal vector S(r) s to give tau( r ) radial       C
C      function                                                      C
C                                                                    C
C        curl[ S(r) s ] = T_c( r ) t                                 C
C                                                                    C
C     where T_c( r ) = dS/dr + S(r)/RAD                              C
C                                                                    C
C       and so tau( r ) =  -( dS/dr + S(r)/RAD )/DSQRT( L(L+1) )     C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C 13: Curl of toroidal vector T( r ) t to give p( r )                C
C                                                                    C
C        curl[ T(r) t ] = Q_c( r ) q   +   S_c( r ) s                C
C                                                                    C
C         where  Q_c( r ) = -T(r) DSQRT( L(L+1) )/RAD   and          C
C                S_c( r ) = -[ dT/dr + T(r)/RAD ]                    C
C                                                                    C
C   hence p( r ) = -T(r)/DSQRT( L(L+1) )                             C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C The matrix DPOM has dimensions ( N1, N2 ). The vector which is     C
C multiplied by DPOM contains NH radial functions, each with NR      C
C radial grid nodes. The radius at grid node IR is XARR( IR ).       C
C The IRth value of function IH is stored in element                 C
C                                                                    C
C   ( IH - 1 )*NR + IR                                               C
C                                                                    C
C The spherical harmonic degree of the radial function is stored     C
C in MLT( IH ).                                                      C
C                                                                    C
C N2 must be equal to NR*NH                                          C
C N1 must be equal to 1 + 2*K where K is the number of diagonals.    C
C                                                                    C
C If IOPT is such that IHD = 0, then K = 0.                          C
C If IOPT is such that IHD = 1, then K = NBN.                        C
C                                                                    C
C The matrix is zero-ed on entry.                                    C
C The terms are added to the matrix between nodes ILNR and IRNR.     C
C                                                                    C
C The finite difference coefficients are stored in the array FDCM    C
C which has dimensions ( NCFM, NR, 1 ) (NCFM=2*NBN) and this         C
C matrix is formed by the subroutine FDCMBD with the following       C
C parameters ...                                                     C
C                                                                    C
C                     NDRVM = 1                                      C
C                     NLMN  = ILNR                                   C
C                     NRMN  = IRNR                                   C
C                     NLMC  = ILNR                                   C
C                     NRMC  = IRNR                                   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N1        : Leading dimension of DPOM.                         C
C     N2        : Second dimension of DPOM.                          C
C     K         : Number of diagonals in banded matrix.              C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of radial functions.                        C
C     NBN       : Number of side nodes for radial derivatives.       C
C     NCFM      : Leading coefficient of FDCM.                       C
C     ILNR      : First grid node to operate upon.                   C
C     IRNR      : Last grid node to operate upon.                    C
C     MLT       : Dim (NH). Element IH contains degree L.            C
C     IOPT      : Choice of matrix operation. See above list.        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim (NR). Radial node values.                      C
C     DPOM      : Dim (N1, N2). Double precision operations matrix.  C
C     FDCM      : Dim (NCFM,NR,1). Finite difference coefficients.   C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VOBMBR( N1, N2, K, NR, NH, NBN, NCFM, ILNR, IRNR,
     1                   MLT, IOPT, XARR, DPOM, FDCM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          N1, N2, K, NR, NH, NBN, NCFM, ILNR, IRNR,
     1                 MLT( NH ), IOPT
      DOUBLE PRECISION XARR( NR ), DPOM( N1, N2 ), FDCM(NCFM,NR,1)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          KARR( 13 ), IHD, IOP, IPARS( 2 ), IMF,
     1                 INARR( 3 ), IH
      DOUBLE PRECISION CVEC( 2 ), ZERO, DPARS( 1 ), FAC
      PARAMETER        ( IOP = 0, ZERO = 0.0d0 )
      EXTERNAL         VOBMAR
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      KARR(  1 ) = 0
      KARR(  2 ) = NBN
      KARR(  3 ) = 0
      KARR(  4 ) = 0
      KARR(  5 ) = NBN
      KARR(  6 ) = 0
      KARR(  7 ) = NBN
      KARR(  8 ) = 0
      KARR(  9 ) = 0
      KARR( 10 ) = NBN
      KARR( 11 ) = 0
      KARR( 12 ) = NBN
      KARR( 13 ) = 0
C
      IF ( K.NE.KARR( IOPT ) ) THEN
        PRINT *,' Subroutine VOBMBR.'
        PRINT *,' IOPT = ', IOPT,' K = ',K
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IHD = 0
      IF ( K.EQ.NBN ) IHD = 1
C
      IF ( N1.NE.(2*K+1) .OR. N2.NE.(NR*NH) ) THEN
        PRINT *,' Subroutine VOBMBR.'
        PRINT *,' N1 = ', N1,' N2 = ', N2
        PRINT *,' NR = ', NR,' NH = ', NH
        PRINT *,' K  = ', K
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Clear the matrix
C
      CALL MATOP( DPOM, ZERO, N1, N2, IOP )
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH
      FAC        = 1.0d0
      IMF        = 1
      IPARS( 1 ) = IOPT
      DO IH = 1, NH
        IPARS( 2 ) = MLT( IH )
        CALL NMLICA( N1, N2, K, K, IOP, IMF, IH, IH, INARR,
     1               IHD, K, ILNR, IRNR, ILNR, IRNR, NCFM,
     2               IMF, VOBMAR, DPOM, FAC, XARR, FDCM, DPARS,
     3               IPARS, NR, CVEC )
      ENDDO
C
      END
C*********************************************************************
C*********************************************************************
C subroutine Adapted Solution Vector ComPLete ************************
C            -       -        -      -  --    ************************
C Steve Gibbons Wed Oct 27 14:10:12 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C A solution vector whose end nodes are defined implicitly by the    C
C boundary conditions rather than given explicitly may have 1 or 2   C
C nodes 'missing' at each boundary.                                  C
C                                                                    C
C ASVCPL uses the finite difference coefficients in SVFDC to fill in C
C these points. It is not usually necessary to do this when further  C
C calculations are required as ASVDR will disregard all such values. C
C However, it is preferable when observing a solution that all       C
C elements are present.                                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     INARR     : Int. parameter array corresponding to ASV          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NR      See INDFUN for details        C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     MHIBC     : MHIBC( is ) describes the inner boundary           C
C                  condition for scheme IS.                          C
C                                                                    C
C  MHIBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( is ) = 7 --> insulating magnetic field.                    C
C                                                                    C
C     MHOBC     : MHOBC( is ) describes the outer boundary           C
C                  condition for scheme IS. (See above for key.)     C
C                                                                    C
C     NFDCM     : Leading order of the array SVFDC.                  C
C                 This must be atleast (2*NBN + 1)                   C
C     NDRVS     : Number of derivatives stored in SVFDC.             C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NBN       : Number of bounding nodes as supplied to SVFDCF.    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     ASV       : Solution vector. Dim ( * )                         C
C                 Length must be atleast NRI*NH                      C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ASVCPL( ASV, NR, NDCS, INARR, MHP, MHIBC, MHOBC,
     1                   NFDCM, NDRVS, NDRVM, NBN, SVFDC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NDCS, INARR( * ), MHP( * ), MHIBC( NDCS ),
     1        MHOBC( NDCS ), NFDCM, NDRVS, NDRVM, NBN
      DOUBLE PRECISION ASV( * ), SVFDC( NFDCM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IH, IS, IR, ISN, IEN, NH, IND, INDFUN, IHD, IBC
      DOUBLE PRECISION DERV( 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IHD = 0
      NH  = INARR( 3 )
C     .
C     . Loop around harmonics
C     .
      DO IH = 1, NH
        IS = MHP( IH )
        IF ( IS.LT.1 .OR. IS.GT.NDCS ) THEN
          PRINT *,' Subroutine ASVCPL.'
          PRINT *,' Harmonic ', IH
          PRINT *,' IS    = ', IS
          PRINT *,' NDCS  = ', NDCS
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
C Inner boundary
C
        IBC = MHIBC( IS )
        IF ( IBC.EQ.1 ) GOTO 50
C       .
        IF ( IBC.EQ.2 .OR. IBC.EQ.3 .OR. IBC.EQ.6
     1         .OR. IBC.EQ.7 ) THEN
           ISN = 1
           IEN = 1
           GOTO 51
        ENDIF
C       .
        IF ( IBC.EQ.4 .OR. IBC.EQ.5 ) THEN
           ISN = 1
           IEN = 2
           GOTO 51
        ENDIF
C       .
        PRINT *,' Subroutine ASVCPL'
        PRINT *,' MHIBC(',IS,') = ', IBC
        PRINT *,' Program aborted.'
        STOP
C       .
 51     CONTINUE
        DO IR = ISN, IEN
C         .
C         . Find index of element
C         .
          IND = INDFUN( IR, IH, INARR )
C         .
          CALL ASVDR ( ASV, IR, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                 NDRVM, DERV, INARR, SVFDC, NDCS )
C         .
          ASV( IND ) = DERV( 1 )
C         .
        ENDDO
C       .
 50     CONTINUE
C
C Outer boundary
C
        IBC = MHOBC( IS )
        IF ( IBC.EQ.1 ) GOTO 60
C       .
        IF ( IBC.EQ.2 .OR. IBC.EQ.3 .OR. IBC.EQ.6
     1         .OR. IBC.EQ.7 ) THEN
           ISN = NR
           IEN = NR
           GOTO 61
        ENDIF
C       .
        IF ( IBC.EQ.4 .OR. IBC.EQ.5 ) THEN
           ISN = NR-1
           IEN = NR
           GOTO 61
        ENDIF
C       .
        PRINT *,' Subroutine ASVCPL'
        PRINT *,' MHOBC(',IS,') = ', IBC
        PRINT *,' Program aborted.'
        STOP
C       .
 61     CONTINUE
        DO IR = ISN, IEN
C         .
C         . Find index of element
C         .
          IND = INDFUN( IR, IH, INARR )
C         .
          CALL ASVDR ( ASV, IR, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                 NDRVM, DERV, INARR, SVFDC, NDCS )
C         .
          ASV( IND ) = DERV( 1 )
C         .
        ENDDO
C       .
 60     CONTINUE
C
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Solution Vector File WriTe ******************************
C            -        -      -    -  -  ******************************
C Steve Gibbons Sat Nov 13 14:30:24 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Writes out a solution vector to a file.                            C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1, 3. INDFUN = ( IR - 1 )*NH + IH       C
C                   IFORMF = 2, 4. INDFUN = ( IH - 1 )*NR + IR       C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C     IFORM     : Specifies how the x values are stored on the file. C
C                 Current values are:-                               C
C                                                                    C
C                   IFORM = 1 --> (5(1PD16.7))                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Dim ( * ) but length atleast NR*NH.                C
C                  Solution vector defined by INARR.                 C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVFWT( INARR, LU, IFORM, SV, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), LU, IFORM
      CHARACTER *(*) FNAME
      DOUBLE PRECISION SV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IWR, ILEN, IFORMF, NR, NH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFORMF = INARR( 1 )
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C
      ILEN   = NR*NH
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine SVFWT.'
        PRINT *,' IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Open file for writing
C
      IWR = 3
      CALL FOPEN ( LU, FNAME, IWR )
C
C Write iformf, nr, nh, iform
C  
       WRITE ( LU, 40 ) IFORMF, NR, NH, IFORM
C
C OK, so write X values ...
C
      IF ( IFORM.EQ.1 ) WRITE ( LU, 41 ) ( SV( I ), I = 1, ILEN )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 40   FORMAT(I5,I5,I5,I5)
 41   FORMAT(5(1PD16.7))
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine Poloidal Velocity Vector ComPLete ***********************
C            -        -        -      -  --    ***********************
C Steve Gibbons Thu Nov  9 10:47:05 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Because of the additional boundary conditions on the poloidal      C
C velocity radial function, the values at grid nodes 2 and NR-1 are  C
C not solved for in the solution vector.                             C
C PVVCPL completes these two values in the solution vector.          C
C                                                                    C
C The solution vector must consist only of poloidal velocity         C
C harmonics, of which there are NH1, and the values of these         C
C radial functions, at NR grid nodes, must be given by the index     C
C                                                                    C
C   ind = ( ih - 1 )*nr + ir                                         C
C                                                                    C
C The arrays PVLC and PVRC have been prepared by a call to PVCCF.    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NH1       : Number of poloidal harmonic functions.             C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV1       : Solution vector for poloidal velocity.             C
C                  Dim (NR*NH1)                                      C
C     PVLC      : Dim (NBN). See PVCCF.                              C
C     PVRC      : Dim (NBN). See PVCCF.                              C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PVVCPL( NR, NH1, NBN, SV1, PVLC, PVRC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NR, NH1, NBN
      DOUBLE PRECISION SV1( NR*NH1 ), PVLC( NBN ), PVRC( NBN )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          I, IND, IH, IR, IND2
      DOUBLE PRECISION TEMP
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO IH = 1, NH1
C
        IR   = 2
        IND2 = ( IH - 1 )*NR + IR
        TEMP = 0.0d0
        DO I = 1, NBN
          IND = IND2 + I
          TEMP = TEMP + PVLC( I )*SV1( IND )
        ENDDO
        SV1( IND2 ) = TEMP
C
        IR   = NR - 1
        IND2 = ( IH - 1 )*NR + IR
        TEMP = 0.0d0
        DO I = 1, NBN
          IND = IND2 - I
          TEMP = TEMP + PVRC( I )*SV1( IND )
        ENDDO
        SV1( IND2 ) = TEMP
C
      ENDDO
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine VECtor OPeration ****************************************
C Steve Gibbons 22.4.97 Fills vector with a constant, multiplies a   C
C                       vector by a constant or adds a constant.     C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IOP	: Type of operation required.                        C
C                  IOP=0  -->  Each element of the vector = CONST    C
C                  IOP=1  -->  Each el. is multiplied by CONST       C
C                  IOP=2  -->  Each el. is added to CONST            C
C     N		: Length of the vector.                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VEC	: Vector - dimension ( N )                           C
C     CONST     : Double precision constant.                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VECOP ( VEC, CONST, N, IOP )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N, IOP
      DOUBLE PRECISION VEC( N ), CONST
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First do the case of making an element constant.
      IF ( IOP.EQ.0 ) THEN
         DO I = 1, N
            VEC ( I ) = CONST
         ENDDO
         RETURN
      ENDIF
C Now do multiplying a vector
      IF ( IOP.EQ.1 ) THEN
         DO I = 1, N
            VEC ( I ) = VEC( I )*CONST
         ENDDO
         RETURN
      ENDIF
C Now do adding a vector
      IF ( IOP.EQ.2 ) THEN
         DO I = 1, N
            VEC ( I ) = VEC( I ) + CONST
         ENDDO
         RETURN
      ENDIF
C____________________________________________________________________C

      PRINT *,' Subroutine VECOP. IOP must be 0, 1 or 2.'
      PRINT *,'Program aborted.'
      STOP
      END
C*********************************************************************


C*********************************************************************
C subroutine New Solution Vector Heat Source Term ********************
C            -   -        -      -    -      -    ********************
C Steve Gibbons Thu Nov  9 15:36:15 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C If \nabla^2 T_0( r ) = C (C is a constant), then T_0 has the       C
C general solution                                                   C
C                                                                    C
C            CB1 r^2                                                 C
C T_0( r ) = ------- - CB2 r^{-1}      (cb1 and cb2 are constants)   C
C               2                                                    C
C                                                                    C
C                                                                    C
C and so                                                             C
C                                                                    C
C              [                                     ]               C
C \nabla T_0 = [  CB1 r  +  CB2 r^{-2}     0     0   ]               C
C              [                       ,      ,      ]               C
C                                                                    C
C                                                                    C
C                                          [         CB2  ]          C
C and so v . \nabla T_0 = l(l+1)P(r) Y_l^m [ CB1 +  ----- ]          C
C                                          [         r^3  ]          C
C                                                                    C
C NSVHST adds this term to the appropriate temperature radial funcs. C
C                                                                    C
C The poloidal velocity terms are stored in the vector SV1.          C
C The radial function IH, at grid node IR is stored at the location  C
C    ind = ( ih - 1 )*nr + ir                                        C
C SV1 must contain ONLY poloidal velocity scalars and the l and m    C
C indices are stored in the arrays ML1 and MM1. There are NH1 pol.   C
C vel. harmonics.                                                    C
C                                                                    C
C The vector FT3 contains the forcing terms to the temperature       C
C equation and is stored in the same way. There are NH3 radial       C
C functions with l and m indices given by ML3 and MM3.               C
C                                                                    C
C The terms are added for nodes IR = 2 to IR = nr - 1.               C
C The solved vector SV1 contains only values for nodes 3 to NR - 2   C
C and so the routine PVVCPL must be called in advance to complete    C
C the poloidal velocity solution vector.                             C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                  See above for key. (corresponds to input vec.)    C
C     NH1       : Number of poloidal velocity radial functions.      C
C     ML1       : Array length ( * ) - atleast length NH1            C
C                  Sph. harm. degree, l, for poloidal velocity.      C
C     MM1       : Array length ( * ) - atleast length NH1            C
C                  Sph. harm. order, m, for cos m phi dep. (p.vel.)  C
C                 -Sph. harm. order, m, for sin m phi dep. (p.vel.)  C
C                                                                    C
C     NH3       : Number of temperature radial functions.            C
C     ML3       : Array length ( * ) - atleast length NH3            C
C                  Sph. harm. degree, l, for temperature.            C
C     MM3       : Array length ( * ) - atleast length NH3            C
C                  Sph. harm. order, m, for cos m phi dep. (temp.)   C
C                 -Sph. harm. order, m, for sin m phi dep. (temp.)   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV1       : Solution vector. Dim ( * ) ( poloidal velocity )   C
C                 Length must be atleast NR*NH1                      C
C     FT3       : Forcing term for heat eqn. Dim ( * )               C
C                 Length must be atleast NR*NH3                      C
C     FAC       : Multiplier of term to be added.                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C     CB1       : Constant - see above.                              C
C     CB2       : Constant - see above.                              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NSVHST( NR, NH1, ML1, MM1, NH3, ML3, MM3, SV1, FT3,
     1                   FAC, XARR, CB1, CB2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH1, ML1(NH1), MM1(NH1), NH3, ML3(NH3), MM3(NH3)
      DOUBLE PRECISION SV1( * ), FT3( * ), FAC, XARR( NR ), CB1, CB2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IR, L1, MMM1, IH1, IH3, IND1, IND3, IBEG1, IBEG3, INDL
      DOUBLE PRECISION DLOW, RAD, D0F, F1, R3, ZERO, DLL1
      PARAMETER ( ZERO = 0.0d0, DLOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     . early exit?
C
      IF ( FAC.EQ.ZERO ) RETURN
C     .
C     . Loop around poloidal velocity harmonics
C     .
      DO IH1 = 1, NH1
        L1   = ML1( IH1 )
        MMM1 = MM1( IH1 )
        INDL = L1*L1 + L1
        DLL1 = DBLE( INDL )
C
C so let's look for the corresponding temperature harmonic in FT3
C
        DO IH3 = 1, NH3
          IF (   ML3( IH3 ).EQ.L1 .AND. MM3( IH3 ).EQ.MMM1  ) THEN
C
C o.k. we've found the corresponding harmonic
C
            IBEG1 = (IH1 - 1)*NR
            IBEG3 = (IH3 - 1)*NR
            DO IR = 2, NR - 1
C
C Find locations of source and destination vectors ...
C
              IND1 = IBEG1 + IR
              IND3 = IBEG3 + IR
C
              RAD = XARR( IR )
              IF ( ABS( RAD ).LT.DLOW ) THEN
                PRINT *,' Subroutine NSVHST.'
                PRINT *,' Rad at node ',IR,' is ',RAD
                PRINT *,' Division by zero imminent.'
                PRINT *,' Program aborted.'
                STOP
              ENDIF
              R3  = RAD*RAD*RAD
              D0F = SV1( IND1 )
              F1  = D0F*DLL1
              FT3( IND3 ) = FT3( IND3 ) + FAC*( CB1 + CB2/R3 )*F1
C
            ENDDO
C           . End of loop ir = 2, nr-1
C
C Jump out of ih3 loop if we have found our harmonic
C
            GOTO 46
C
          ENDIF
C         . check ml3(ih3).eq.l1 and mm3(ih3).eq.mm1
        ENDDO
C       . End of loop ih3 = 1, nh3
C
 46   CONTINUE
      ENDDO
C     . End of loop ih1 = 1, nh1
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine New Solution Vector Buoyancy Term Add *******************
C            -   -        -      -        -    -   *******************
C Steve Gibbons Fri Nov 10 07:44:00 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Adds the buoyancy term from SV3, a solution vector containing the  C
C temperature function, to the toroidal part of the vorticity        C
C equation in FT1. CH multiplies the terms added.                    C
C                                                                    C
C The temperature terms are stored in the vector SV3.                C
C The radial function IH, at grid node IR is stored at the location  C
C    ind = ( ih - 1 )*nr + ir                                        C
C SV3 must contain ONLY temperature scalars and the l and m indices  C
C are stored in the arrays ML3 and MM3. There are NH3 temp. harm.s   C
C                                                                    C
C The vector FT1 contains the forcing terms to the toroidal part     C
C of the vorticity equation and is stored in the same way. There     C
C are NH1 functions with l and m indices given by ML1 and MM1.       C
C                                                                    C
C The terms are added for nodes IR = 3 to IR = nr - 2.               C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                  See above for key. (corresponds to input vec.)    C
C     NH3       : Number of poloidal velocity radial functions.      C
C     ML3       : Array length ( * ) - atleast length NH3            C
C                  Sph. harm. degree, l, for poloidal velocity.      C
C     MM3       : Array length ( * ) - atleast length NH3            C
C                  Sph. harm. order, m, for cos m phi dep. (temp.)   C
C                 -Sph. harm. order, m, for sin m phi dep. (temp.)   C
C                                                                    C
C     NH1       : Number of func.s in tor. vorticity eqn.            C
C     ML1       : Array length ( * ) - atleast length NH1            C
C                  Sph. harm. degree, l, for tor. vort. forcing t.   C
C     MM1       : Array length ( * ) - atleast length NH3            C
C                  Sph. harm. order, m, for cos m phi dep. (t. vort) C
C                 -Sph. harm. order, m, for sin m phi dep. (t. vort) C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV3       : Solution vector. Dim ( * ) ( temperature )         C
C                 Length must be atleast NR*NH3                      C
C     FT1       : Forcing term for toroidal part of vorticity eqn.   C
C                  Dim ( * ). Length must be atleast NR*NH1          C
C     CH        : Multiplier of term to be added.                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NSVBTA( NR, NH3, ML3, MM3, NH1, ML1, MM1, SV3, FT1,
     1                   CH )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH3, ML3(NH3), MM3(NH3), NH1, ML1(NH1), MM1(NH1)
      DOUBLE PRECISION SV3( * ), FT1( * ), CH
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IR, L3, MMM3, IH3, IH1, IND3, IND1, IBEG3, IBEG1
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     . early exit?
C
      IF ( CH.EQ.ZERO ) RETURN
C     .
C     . Loop around temperature harmonics
C     .
      DO IH3 = 1, NH3
        L3   = ML3( IH3 )
        MMM3 = MM3( IH3 )
C
C so let's look for the corresponding tor. vort. harmonic in FT1
C
        DO IH1 = 1, NH1
          IF (   ML1( IH1 ).EQ.L3 .AND. MM1( IH1 ).EQ.MMM3  ) THEN
C
C o.k. we've found the corresponding harmonic
C
            IBEG3 = (IH3 - 1)*NR
            IBEG1 = (IH1 - 1)*NR
            DO IR = 3, NR - 2
C
C Find locations of source and destination vectors ...
C
              IND3 = IBEG3 + IR
              IND1 = IBEG1 + IR
C
              FT1( IND1 ) = FT1( IND1 ) + CH*SV3( IND3 )
C
            ENDDO
C           . End of loop ir = 3, nr-2
C
C Jump out of ih1 loop if we have found our harmonic
C
            GOTO 46
C
          ENDIF
C         . check ml1(ih1).eq.l3 and mm1(ih1).eq.mm3
        ENDDO
C       . End of loop ih1 = 1, nh1
C
 46   CONTINUE
      ENDDO
C     . End of loop ih3 = 1, nh3
C
      RETURN
      END
C*********************************************************************


C*********************************************************************
C subroutine Radial Spectrally Decomposed Vector 2 xtra special C ****
C            -      -          -          -      -              - ****
C Steve Gibbons Mon Feb 12 10:14:54 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Prepares the arrays                                                C
C                                                                    C
C         FTFP( 3, NPH, NTHP ), INFP( 2, NPH ),                      C
C         FTFT( 2, NTH, NTHP ), INFT( 2, NTH )                       C
C                                                                    C
C for use by the vector transform routine RSDV2D.                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NTHP      : Number of theta points.                            C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     LH        : Maximum spherical harmonic degree, l.              C
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
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INFP      : Dim(2,NPH). Locations in FTF2/3 req. by RSDV2D     C
C                                                                    C
C     INFT      : Dim(2,NTH). Locations in FTF2/3 req. by RSDV2D     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     FTFP     : Dim (3,NPH,NTHP). Coeff.s required by RSDV2D        C
C     FTFT     : Dim (2,NTH,NTHP). Coeff.s required by RSDV2D        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RSDV2C( NTHP, M0, LH, NPH, MLP, MMP, NTH, MLT, MMT,
     1                   GAUX, PA, DPA, FTFP, FTFT, INFP, INFT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NTHP, M0, LH, NPH, MLP( NPH ),
     1                 MMP( NPH ), NTH, MLT( NTH ), MMT( NTH )
      INTEGER          INFP( 2, NPH ), INFT( 2, NTH )
      DOUBLE PRECISION GAUX( NTHP ),
     1                 FTFP( 3, NPH, NTHP ), FTFT( 2, NTH, NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER          ITHETA, IP, ICS, L, M, MPS, IHP, IHT,
     1                 INDCOS, INDSIN
      DOUBLE PRECISION ZERO, X, SINE, TERM1, DALF, DALFD,
     1                 DLFAC, DLSQR
C
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C   
C ............. start to loop around theta points ...................
      DO ITHETA = 1, NTHP
C
       X = GAUX( ITHETA )
C
       SINE = DSQRT( 1.0d0 - X*X )
C
C ............. start to loop around Harmonics ......................
C ............. First do poloidal harmonics (Q and S radial func.s) .
       DO IHP = 1, NPH
C        .
C        . Calculate factors including L.
C        . Protect against division by zero.
C        .
         L      = MLP( IHP )
         IF ( L.LT.1 ) THEN
           PRINT *,' Subroutine RSDV2C.'
           PRINT *,' MLP(',IHP,') = ', L
           PRINT *,' Program aborted.'
           STOP
         ENDIF
         DLFAC  = DBLE( L )
         DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
C        .
C        . Find wavenumber, m.
C        .
         M      = MMP( IHP )
         ICS    = 1
C        .
C        . Modify M and ICS if the harmonic is sin (m phi) dep.
C        .
         IF ( M.LT.0   ) ICS = 2
         IF ( ICS.EQ.2 ) M   = -M
C        .
C        . Store associated Legendre Functions
C        .
         IP     = L*(L+1)/2+M+1
         DALF   = PA( IP, ITHETA )
         DALFD  = DPA( IP, ITHETA )
C        .
C        . Make sure that harmonic is
C        . compatible with the symmetry imposition
C        . MPS (pseudo M) = M/M0
C        .
         MPS = M/M0
         IF ( MPS*M0.NE.M ) THEN
           PRINT *,' Subroutine RSDV2C.'
           PRINT *,' MMP(',IHP,') = ', MMP( IHP )
           PRINT *,' Program aborted.'
           STOP
         ENDIF
C        .
         IF ( ICS.EQ.1 ) THEN
C          .
C          . We are dealing with a cos ( m phi ) term
C          .
           IF ( M.EQ.0 ) THEN
             INDCOS         = 1
             INDSIN         = INDCOS + 1
             FTFP( 1, IHP, ITHETA ) = DALF
             FTFP( 2, IHP, ITHETA ) = DALFD/DLSQR
             FTFP( 3, IHP, ITHETA ) = ZERO
             INFP( 1, IHP ) = INDCOS
             INFP( 2, IHP ) = INDSIN
           ELSE
             INDCOS         = 1 + 2*M/M0
             INDSIN         = INDCOS + 1
             TERM1          = (-1.0d0)*DBLE( M )*DALF/SINE
             FTFP( 1, IHP, ITHETA ) = DALF
             FTFP( 2, IHP, ITHETA ) = DALFD/DLSQR
             FTFP( 3, IHP, ITHETA ) = TERM1/DLSQR
             INFP( 1, IHP ) = INDCOS
             INFP( 2, IHP ) = INDSIN
           ENDIF
C          .
         ELSE
C          .
C          . We are dealing with a sin ( m phi ) term
C          . m should _never_ be zero here!!
C          .
           INDCOS         = 1 + 2*M/M0
           INDSIN         = INDCOS + 1
           TERM1          = DBLE( M )*DALF/SINE
           FTFP( 1, IHP, ITHETA ) = DALF
           FTFP( 2, IHP, ITHETA ) = DALFD/DLSQR
           FTFP( 3, IHP, ITHETA ) = TERM1/DLSQR
           INFP( 1, IHP ) = INDSIN
           INFP( 2, IHP ) = INDCOS
C          .
         ENDIF
C        .
 70      CONTINUE
       ENDDO
C      .
C      . End loop ihp = 1, nhp
C
C............. Now do toroidal harmonics (T radial func.s) .
       DO IHT = 1, NTH
C        .
C        . Calculate factors including L.
C        . Protect against division by zero.
C        .
         L      = MLT( IHT )
         IF ( L.LT.1 ) THEN
           PRINT *,' Subroutine RSDV2C.'
           PRINT *,' MLT(',IHT,') = ', L
           PRINT *,' Program aborted.'
           STOP
         ENDIF
         DLFAC  = DBLE( L )
         DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
C        .
C        . Find wavenumber, m.
C        .
         M      = MMT( IHT )
         ICS    = 1
C        .
C        . Modify M and ICS if the harmonic is sin (m phi) dep.
C        .
         IF ( M.LT.0   ) ICS = 2
         IF ( ICS.EQ.2 ) M   = -M
C        .
C        . Store associated Legendre Functions
C        .
         IP     = L*(L+1)/2+M+1
         DALF   = PA( IP, ITHETA )
         DALFD  = DPA( IP, ITHETA )
C        .
C        . Make sure that harmonic is
C        . compatible with the symmetry imposition
C        . MPS (pseudo M) = M/M0
C        .
         MPS = M/M0
         IF ( MPS*M0.NE.M ) THEN
           PRINT *,' Subroutine RSDV2C.'
           PRINT *,' MMT(',IHT,') = ', MMT( IHT )
           PRINT *,' Program aborted.'
           STOP
         ENDIF
C        .
         IF ( ICS.EQ.1 ) THEN
C          .
C          . We are dealing with a cos ( m phi ) term
C          .
           IF ( M.EQ.0 ) THEN
             INDCOS         = 1
             INDSIN         = INDCOS + 1
             FTFT( 1, IHT, ITHETA ) = ZERO
             FTFT( 2, IHT, ITHETA ) = DALFD/DLSQR
             INFT( 1, IHT ) = INDSIN
             INFT( 2, IHT ) = INDCOS
           ELSE
             INDCOS         = 1 + 2*M/M0
             INDSIN         = INDCOS + 1
             TERM1          = DBLE( M )*DALF/SINE
             FTFT( 1, IHT, ITHETA ) = TERM1/DLSQR
             FTFT( 2, IHT, ITHETA ) = DALFD/DLSQR
             INFT( 1, IHT ) = INDSIN
             INFT( 2, IHT ) = INDCOS
           ENDIF
C          .
         ELSE
C          .
C          . We are dealing with a sin ( m phi ) term
C          . m should _never_ be zero here!!
C          .
           INDCOS         = 1 + 2*M/M0
           INDSIN         = INDCOS + 1
           TERM1          = (-1.0d0)*DBLE( M )*DALF/SINE
           FTFT( 1, IHT, ITHETA ) = TERM1/DLSQR
           FTFT( 2, IHT, ITHETA ) = DALFD/DLSQR
           INFT( 1, IHT ) = INDCOS
           INFT( 2, IHT ) = INDSIN
C          .
         ENDIF
C        .
 71      CONTINUE
       ENDDO
C        .
C        . End loop iht = 1, nht
C        .
C ............. ended looping around Harmonics ......................
C
      ENDDO
C ............. ended looping around theta points ...................
C
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine Scalar Function 2 Vector Gradient xtra special vector C *
C            -      -        - -      -                            - *
C Steve Gibbons Tue Feb 13 09:28:18 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Prepares the 2 arrays                                              C
C                                                                    C
C        VGFA( 3, NH3, NTHP )                                        C
C        IVGFA( 2, NH3 )                                             C
C                                                                    C
C  for use by the gradient transform routine SF2VGD.                 C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NH3       : Number of scalar functions                         C
C     ML3       : Array dim ( NH3 ). Sph. harm degree, L.            C
C     MM3       : Array dim ( NH3 ). Sph. harm order, M, or -M.      C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     NTHP      : Number of theta points.                            C
C                                                                    C
C     IVGFA     : Dim ( 2, NH3 ). Locations for SF2VGD.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     VGFA      : Dim ( 3, NH3, NTHP ). Coeff.s for SF2VGD           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SF2VGC( LH, NH3, ML3, MM3, M0, NTHP,
     1                   GAUX, PA, DPA, VGFA, IVGFA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          LH, NH3, ML3( NH3 ), MM3( NH3 ), M0, NTHP,
     1                 IVGFA( 2, NH3 )
      DOUBLE PRECISION GAUX( NTHP ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2, NTHP ),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2, NTHP )
      DOUBLE PRECISION VGFA( 3, NH3, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ITHETA, IP, ICS, L, M, MPS, IH3,
     1                 INDCOS, INDSIN
      DOUBLE PRECISION X, SINE, DALF, DALFD, TERM1
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C ............. start to loop around theta points ...................
      DO ITHETA = 1, NTHP
        X      = GAUX( ITHETA )
        SINE   = DSQRT( 1.0d0 - X*X )
C
C ............. start to loop around Harmonics ......................
C
        DO IH3 = 1, NH3 
C
           L      = ML3( IH3 )
C
C First do monopole term
C (We do not fill in any coefficients for this term)
C
           IF ( L.EQ.0 ) GOTO 50
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
           IF ( MPS*M0.NE.M ) THEN
             PRINT *,' Subroutine SF2VGC.'
             PRINT *,' MM3(',IH3,') = ', MM3( IH3 )
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
               INDSIN         = INDCOS + 1
               VGFA( 1, IH3, ITHETA ) = DALF
               VGFA( 2, IH3, ITHETA ) = DALFD
               VGFA( 3, IH3, ITHETA ) = 0.0d0
               IVGFA( 1, IH3 ) = INDCOS
               IVGFA( 2, IH3 ) = INDSIN
             ELSE
               INDCOS         = 1 + 2*M/M0
               INDSIN         = INDCOS + 1
               TERM1          = (-1.0d0)*DBLE( M )*DALF/SINE
               VGFA( 1, IH3, ITHETA ) = DALF
               VGFA( 2, IH3, ITHETA ) = DALFD
               VGFA( 3, IH3, ITHETA ) = TERM1
               IVGFA( 1, IH3 ) = INDCOS
               IVGFA( 2, IH3 ) = INDSIN
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS         = 1 + 2*M/M0
             INDSIN         = INDCOS + 1
             TERM1          = DBLE( M )*DALF/SINE
             VGFA( 1, IH3, ITHETA ) = DALF
             VGFA( 2, IH3, ITHETA ) = DALFD
             VGFA( 3, IH3, ITHETA ) = TERM1
             IVGFA( 1, IH3 ) = INDSIN
             IVGFA( 2, IH3 ) = INDCOS
C            .
           ENDIF
C          .
 50     CONTINUE
        ENDDO
C ............. ended looping around Harmonics ......................
      ENDDO
C ............. ended looping around theta points ...................
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Scalar Function 2 Spectral Decomposition add xsv C ******
C            -      -        - -        -                     - ******
C Steve Gibbons Thu Feb 22 13:45:12 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Pre-calculates the array SF2SA( NH3, NTHP ) which is used by       C
C scalar real --> spectral transform routine  SF2SDD.                C
C It also produces the integer array ISF2S( NH3 ) for SF2SDD.        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NH3       : Number of scalar functions                         C
C     ML3       : Array dim ( NH3 ). Sph. harm degree, L.            C
C     MM3       : Array dim ( NH3 ). Sph. harm order, M, or -M.      C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     NTHP      : Number of theta points.                            C
C     NPHP      : Number of phi points.                              C
C                                                                    C
C     ISF2S     : Dim ( NH3 ). Array for use by SF2SDD.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C                                                                    C
C     GAUW      : Gauss weights as calculated by GAUWTS. ( NTHP )    C
C                                                                    C
C     FAC       : Coefficient of function to be added.               C
C                                                                    C
C     SF2SA     : Dim ( NH3, NTHP ). Array for use by SF2SDD.        C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SF2SDC( LH, NH3, ML3, MM3, M0, NTHP, NPHP, PA, GAUW,
     1                   FAC, SF2SA, ISF2S )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          LH, NH3, ML3( NH3 ), MM3( NH3 ), M0, NTHP,
     1                 NPHP, ISF2S( NH3 )
      DOUBLE PRECISION PA( ( LH + 1 )*( LH + 2 )/2, NTHP ),
     1                 GAUW( NTHP ), FAC, SF2SA( NH3, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ITHETA, IH3, ICS, L, M, IP, MPS
      DOUBLE PRECISION FACL, DALF, WEIGHT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C ............... begin looping around theta points ..................
      DO ITHETA = 1, NTHP
        WEIGHT = GAUW( ITHETA )
C       .
C       . Loop around harmonics
C       .
        DO IH3 = 1, NH3
C
           L      = ML3( IH3 )
C
C First do monopole term
C
           IF ( L.EQ.0 ) THEN
             SF2SA( IH3, ITHETA ) = 0.5d0*FAC*WEIGHT/DBLE( NPHP )
             GOTO 70
           ENDIF
C
           FACL   = 0.25d0*DBLE( 2*L + 1 )
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
C          .
C          . Make sure that harmonic is
C          . compatible with the symmetry imposition
C          . MPS (pseudo M) = M/M0
C          .
           MPS = M/M0
           IF ( MPS*M0.NE.M ) THEN
             PRINT *,' Subroutine SF2SDC.'
             PRINT *,' MM3(',IH3,') = ', MM3( IH3 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           ISF2S( IH3 )         = ICS + 2*M/M0
C          .
           SF2SA( IH3, ITHETA ) = FAC*WEIGHT*DALF*FACL
C          .
 70     CONTINUE
        ENDDO
C       .      Ended loop ih3 = 1, nh3
C       .
      ENDDO
C ............. ended looping around theta points ...................
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Xtra Special Vector function 2 Spectrally Decomp. vec A *
C            -    -       -                 -          -           - *
C Steve Gibbons Tue Feb 13 14:02:50 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Creates the four arrays                                            C
C                                                                    C
C       PFA( 3, NPH, NTHP ), JPFA( 2, NPH ),                         C
C       TFA( 2, NTH, NTHP ), JTFA( 2, NTH )                          C
C                                                                    C
C required by the transform routine XSVSDD.                          C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NTHP      : Number of theta points.                            C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NPH       : Number of poloidal harmonics [number of radial     C
C                  functions Q(r) and S(r)]                          C
C     MLP       : Array dim ( NPH ). Sph. harm degree, L.            C
C     MMP       : Array dim ( NPH ). Sph. harm order, M, or -M.      C
C     NTH       : Number of toroidal harmonics [number of radial     C
C                  functions T(r)]                                   C
C     MLT       : Array dim ( NTH ). Sph. harm degree, L.            C
C     MMT       : Array dim ( NTH ). Sph. harm order, M, or -M.      C
C                                                                    C
C     JPFA      : Dim(2,NPH). Locations for use by XSVSDD            C
C     JTFA      : Dim(2,NTH). Locations for use by XSVSDD            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUX      : Dim (NTHP). Gauss weights from the GAUWTS routine. C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                  Both PA and DPA are formed by the routine SCHNLA. C
C                                                                    C
C     PFA       : Dim ( 3, NPH, NTHP ). Coeffs for XSVSDD.           C
C     TFA       : Dim ( 2, NTH, NTHP ). Coeffs for XSVSDD.           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XSVSDC( NTHP, M0, LH, NPH, MLP, MMP, NTH, MLT, MMT,
     1                   GAUX, GAUW, PA, DPA, PFA, TFA, JPFA, JTFA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHP, M0, LH, NPH, MLP( NPH ), MMP( NPH ), NTH,
     1        MLT( NTH ), MMT( NTH ), JPFA( 2, NPH ), JTFA( 2, NTH )
      DOUBLE PRECISION GAUX( NTHP ), GAUW( NTHP ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP ),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP )
      DOUBLE PRECISION PFA( 3, NPH, NTHP ), TFA( 2, NTH, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHETA, MPS, ICS, IHP, IHT, IP, L, M, INDCOS, INDSIN
      DOUBLE PRECISION X, SINE, TERM, WEIGHT, W1, W2,
     1                 DALF, DALFD, DLFAC, DLSQR
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C ....................................................................
C ...... now start to loop around theta points .......................
      DO ITHETA = 1, NTHP
C
         X = GAUX( ITHETA )
         SINE   = DSQRT( 1.0d0 - X*X )
         WEIGHT = GAUW( ITHETA )
C        .
C ...................................................................
C .                  Now let's loop around the Harmonics..          .
C ...................................................................
C ............. First do poloidal harmonics (Q and S radial func.s) .
         DO IHP = 1, NPH
C          .
C          . Calculate factors including L.
C          . Protect against division by zero.
C          .
           L      = MLP( IHP )
           IF ( L.LT.1 ) THEN
             PRINT *,' Subroutine XSVSDC.'
             PRINT *,' MLP(',IHP,') = ', L
             PRINT *,' Program aborted.'
             STOP
           ENDIF
           DLFAC  = DBLE( L )
           DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
           W1     = 0.25d0*WEIGHT*(DLFAC+DLFAC+1.0d0)
           W2     = W1/DLSQR
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
             PRINT *,' Subroutine XSVSDC.'
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
               INDCOS  = 1
               INDSIN  = INDCOS + 1
               PFA( 1, IHP, ITHETA ) = W1*DALF
               PFA( 2, IHP, ITHETA ) = W2*DALFD
               PFA( 3, IHP, ITHETA ) = 0.0d0
               JPFA( 1, IHP ) = INDCOS
               JPFA( 2, IHP ) = INDSIN
             ELSE
               INDCOS  = 1 + 2*M/M0
               INDSIN  = INDCOS + 1
               TERM    = (-1.0d0)*W2*DBLE( M )*DALF/SINE
               PFA( 1, IHP, ITHETA ) = W1*DALF
               PFA( 2, IHP, ITHETA ) = W2*DALFD
               PFA( 3, IHP, ITHETA ) = TERM
               JPFA( 1, IHP ) = INDCOS
               JPFA( 2, IHP ) = INDSIN
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS  = 1 + 2*M/M0
             INDSIN  = INDCOS + 1
             TERM    = W2*DBLE( M )*DALF/SINE
             PFA( 1, IHP, ITHETA ) = W1*DALF
             PFA( 2, IHP, ITHETA ) = W2*DALFD
             PFA( 3, IHP, ITHETA ) = TERM
             JPFA( 1, IHP ) = INDSIN
             JPFA( 2, IHP ) = INDCOS
C            .
           ENDIF
C          .
         ENDDO
C        .
C        . End loop ihp = 1, nhp
C        .
C............. Now do toroidal harmonics (T radial func.s) .
         DO IHT = 1, NTH
C          .
C          . Calculate factors including L.
C          . Protect against division by zero.
C          .
           L      = MLT( IHT )
           IF ( L.LT.1 ) THEN
             PRINT *,' Subroutine XSVSDC.'
             PRINT *,' MLT(',IHT,') = ', L
             PRINT *,' Program aborted.'
             STOP
           ENDIF
           DLFAC  = DBLE( L )
           DLSQR  = DSQRT( DLFAC*DLFAC + DLFAC )
           W1     = 0.25d0*WEIGHT*(DLFAC+DLFAC+1.0d0)
           W2     = W1/DLSQR
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
             PRINT *,' Subroutine XSVSDC.'
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
               INDSIN         = INDCOS + 1
               TFA( 1, IHT, ITHETA ) = 0.0d0
               TFA( 2, IHT, ITHETA ) = W2*DALFD
               JTFA( 1, IHT ) = INDSIN
               JTFA( 2, IHT ) = INDCOS
             ELSE
               INDCOS  = 1 + 2*M/M0
               INDSIN  = INDCOS + 1
               TERM    = W2*DBLE( M )*DALF/SINE
               TFA( 1, IHT, ITHETA ) = TERM
               TFA( 2, IHT, ITHETA ) = W2*DALFD
               JTFA( 1, IHT ) = INDSIN
               JTFA( 2, IHT ) = INDCOS
             ENDIF
C            .
           ELSE
C            .
C            . We are dealing with a sin ( m phi ) term
C            . m should _never_ be zero here!!
C            .
             INDCOS  = 1 + 2*M/M0
             INDSIN  = INDCOS + 1
             TERM    = (-1.0d0)*W2*DBLE( M )*DALF/SINE
             TFA( 1, IHT, ITHETA ) = TERM
             TFA( 2, IHT, ITHETA ) = W2*DALFD
             JTFA( 1, IHT ) = INDCOS
             JTFA( 2, IHT ) = INDSIN
C            .
           ENDIF
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
C
      RETURN
      END
C*********************************************************************
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
C*********************************************************************
C subroutine Scalar Function 2 Spectral Decomposition add xsv ********
C            -      -        - -        -                 -   ********
C Steve Gibbons Thu Nov 30 15:24:44 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C The function XSV( NCMX, NPHP, NTHP, NR ) defines, in real space,   C
C the values of a function. SF( ICMP, iphi, ithe, ir ) is the value  C
C at  RAD = xarr( ir )                                               C
C                                                                    C
C  THETA = ACOS[ GAUX( ithe ) ]                                      C
C                                 and                                C
C  PHI   = (iphi-1)*DELTAP                                           C
C                                 with                               C
C  deltap = 2*pi/(NPHP*M0).                                          C
C                                                                    C
C                                                                    C
C SV3 is a solution vector containing a spectral decomposition       C
C of a scalar function, psi, at NR radial grid nodes such that       C
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
C  SV3[ ( ih3 - 1 )*NR + IR ]   for IFORMF = 4    or                 C
C  SV3[ ( IR  - 1 )*NH3 + ih3 ]   for IFORMF = 3                     C
C                                                                    C
C If the array SF denotes values of a function at grid node IR       C
C then SF2SDD will add FAC*SF to the appropriate elements in SV3.    C
C                                                                    C
C This routine requires the arrays SF2SA( NH3, NTHP ) and            C
C ISF2S( NH3 ) to be pre-calculated by a call to SF2SDC              C
C (this array must incorporate the scaling FAC).                     C
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
C     ICM       : Component of XSV which stores scalar function.     C
C                                                                    C
C     IFORMF    : Defines the arrangement of the solution vec.s      C
C                   (see above: either set to 3 or 4)                C
C                                                                    C
C     ISF2S     : Dim ( NH3 ). Array supplied by SF2SDC.             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XSV       : Xtra Special Function: is an array containing a    C
C                  function over a set of theta, phi and r points.   C
C                    Dimensions are                                  C
C                      ( NCMX, NPHP, NTHP, NR )                      C
C                                                                    C
C     SV3       : Dim ( NR*NH3 ) Radial functions psi_{ip3}(r)       C
C                                                                    C
C     FTF       : Work array: dim ( 2*NPHP )                         C
C                                                                    C
C     SF2SA     : Dim ( NH3, NTHP ). Array prepared by SF2SDC.       C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SF2SDD( ILNR, IRNR, NR, NH3, ML3, NTHP,
     1                   NPHP, NCMX, ICM, SF, SV3, FTF, SF2SA,
     2                   IFORMF, ISF2S )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          ILNR, IRNR, NR, NH3, ML3( NH3 ),
     1                 NTHP, NPHP, NCMX, ICM, IFORMF, ISF2S( NH3 )
      DOUBLE PRECISION SF( NCMX, NPHP, NTHP, NR ),
     1                 FTF( 2*NPHP ), SV3( * ), SF2SA( NH3, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ISIGN, ITHETA, LENV, IPHI, INDLOC, IH3,
     1                 L, INDEX, IR, INDINC
      DOUBLE PRECISION ZERO, ZCOEF, FAC2
      PARAMETER        ( ZERO = 0.0d0, ISIGN = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
        PRINT *,' Subroutine SF2SDD'
        PRINT *,' IFORMF = ', IFORMF
        PRINT *,' Program aborted,'
        STOP
      ENDIF
C
C Check value of ICM
C
      IF ( ICM.LT.1 .OR. ICM.GT.NCMX ) THEN
        PRINT *,' Subroutine SF2SDD.'
        PRINT *,' ICM = ', ICM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      LENV = 2*NPHP
C ............... begin looping around radial grid nodes .............
      DO IR = ILNR, IRNR
C
C ............... begin looping around theta points ..................
       DO ITHETA = 1, NTHP
        CALL DVECZ( FTF, LENV )
        ZCOEF  = ZERO
        INDLOC = -1
        DO IPHI = 1, NPHP
          FAC2   = SF( ICM, IPHI, ITHETA, IR )
          INDLOC = INDLOC + 2
          FTF( INDLOC ) = FAC2
          ZCOEF = ZCOEF + FAC2
        ENDDO
C       .
C       . Do Fast Fourier Transform
C       .
        CALL FFTRLV( FTF, NPHP, ISIGN )
C       .
C       . ftf now contains spectral coefficients
C       . Loop around harmonics
C       .
        IF ( IFORMF.EQ.3 ) THEN
          INDLOC = IR*NH3 - NH3
          INDINC = 1
        ENDIF
C       .
        IF ( IFORMF.EQ.4 ) THEN
          INDLOC = IR - NR
          INDINC = NR
        ENDIF
C       .
        DO IH3 = 1, NH3
           INDLOC = INDLOC + INDINC
C
           L      = ML3( IH3 )
C
C First do monopole term
C
           IF ( L.EQ.0 ) THEN
             FAC2 = SF2SA( IH3, ITHETA )*ZCOEF
             SV3( INDLOC ) = SV3( INDLOC ) + FAC2
             GOTO 70
           ENDIF
C          .
           INDEX = ISF2S( IH3 )
C          .
           FAC2 = SF2SA( IH3, ITHETA )*FTF( INDEX )
           SV3( INDLOC ) = SV3( INDLOC ) + FAC2
C          .
 70     CONTINUE
        ENDDO
C       .      Ended loop ih3 = 1, nh3
C       .
       ENDDO
C ............. ended looping around theta points ...................
C     .
      ENDDO
C ............. ended looping around radial grid nodes ..............
C     .
      RETURN
      END
C*********************************************************************
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
C*********************************************************************
C subroutine Non-Magnetic Convection Xtra Special Evaluation *********
C            -   -        -          -    -       -          *********
C Steve Gibbons Thu Nov 30 16:24:54 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C We have received an array XSV( NCMX, NPHP, NTHP, NR ) into which   C
C we have placed the following:-                                     C
C                                                                    C
C  ICM             Function                                          C
C  ---             --------                                          C
C                                                                    C
C    1             v_r                                               C
C    2             v_{theta}                                         C
C    3             v_{phi}                                           C
C    4             (Grad Theta)_r                                    C
C    5             (Grad Theta)_{theta}                              C
C    6             (Grad Theta)_{phi}                                C
C    7             (curl v)_r                                        C
C    8             (curl v)_{theta}                                  C
C    9             (curl v)_{phi}                                    C
C                                                                    C
C  Into ICM  = 10, we put v.Grad( theta ) and we put the r, theta    C
C and phi components of   (-CG)*( k x V ) + CF*( V x curl V )        C
C                    respectively into ICM  = 11, 12 and 13.         C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHP      : Number of theta points.                            C
C     NPHP      : Number of phi points.                              C
C     NR        : Number of radial grid nodes.                       C
C     ILNR      : Lowest radial grid node.                           C
C     IRNR      : Highest radial grid node.                          C
C     NCMX      : Maximum number of components stored in XSV         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XSV       : eXtra Special Vector. An array of dimensions       C
C                ( NCMX, NPHP, NTHP, NR)                             C
C                 At the radial grid node, IR, and theta point       C
C                 ITHE and phi point IPHI, the components listed     C
C                 above are stored in                                C
C                   XSV( ICM, IPHI, ITHE, IR )                       C
C                 where ICM corresponds to the number above.         C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     CG        : Coefficient of Coriolis term.                      C
C     CF        : Coefficient of ( V x curl V ) term.                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NMCXSE( NTHP, NPHP, NR, ILNR, IRNR, NCMX, XSV, GAUX,
     1                   CG, CF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NTHP, NPHP, NR, ILNR, IRNR, NCMX
      DOUBLE PRECISION XSV( NCMX, NPHP, NTHP, NR ), GAUX( NTHP ),
     1                 CG, CF
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      INTEGER          ITHE, IPHI, IR
      DOUBLE PRECISION VRAD, VTHE, VPHI, GRAD, GTHE, GPHI,
     1                 CRAD, CTHE, CPHI, SCALF, OUTPHI, OUTRAD,
     2                 OUTTHE, COSTH, SINTH, COSTHG, SINTHG, COSTHM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of NCMX
C
      IF ( NCMX.LT.13 ) THEN
        PRINT *,' Subroutine NMCXSE.'
        PRINT *,' NCMX = ', NCMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Loop around radial grid nodes
C
      DO IR = ILNR, IRNR
C      ............. loop around theta and phi points ................
       DO ITHE = 1, NTHP
         COSTH  = GAUX( ITHE )
         SINTH  = DSQRT( 1.0d0 - COSTH*COSTH )
         COSTHG = CG*COSTH
         COSTHM = COSTHG*(-1.0d0)
         SINTHG = CG*SINTH
         DO IPHI = 1, NPHP
C          .
C          . Fill all the scalar variables
C          . First, the velocity components
C          .
           VRAD   = XSV( 1, IPHI, ITHE, IR )
           VTHE   = XSV( 2, IPHI, ITHE, IR )
           VPHI   = XSV( 3, IPHI, ITHE, IR )
C          .
C          . The Grad( theta ) components
C          .
           GRAD   = XSV( 4, IPHI, ITHE, IR )
           GTHE   = XSV( 5, IPHI, ITHE, IR )
           GPHI   = XSV( 6, IPHI, ITHE, IR )
C          .
C          . Now curl v components
C          .
           CRAD   = XSV( 7, IPHI, ITHE, IR )*CF
           CTHE   = XSV( 8, IPHI, ITHE, IR )*CF
           CPHI   = XSV( 9, IPHI, ITHE, IR )*CF
C          .
           SCALF  = VRAD*GRAD + VTHE*GTHE + VPHI*GPHI
           XSV( 10, IPHI, ITHE, IR ) = SCALF
C          .
C          . Evaluate Coriolis terms
C          .
c          OUTRAD = CG*SINTH*VPHI
c          OUTTHE = CG*COSTH*VPHI
c          OUTPHI = (-1.0d0)*CG*COSTH*VTHE - CG*SINTH*VRAD
           OUTRAD = SINTHG*VPHI
           OUTTHE = COSTHG*VPHI
           OUTPHI = COSTHM*VTHE - SINTHG*VRAD
C          .
C          . Finally evaluate forcing terms
C          .
           XSV( 11, IPHI, ITHE, IR ) = OUTRAD +
     1           VTHE*CPHI - VPHI*CTHE
C
           XSV( 12, IPHI, ITHE, IR ) = OUTTHE +
     1           VPHI*CRAD - VRAD*CPHI
C
           XSV( 13, IPHI, ITHE, IR ) = OUTPHI +
     1           VRAD*CTHE - VTHE*CRAD
C          .
         ENDDO
       ENDDO
C ............. ended looping around theta, phi points ...............
      ENDDO
C
C End loop around radial grid nodes
C
      RETURN
      END
C*********************************************************************
C********************************************************************
C subroutine File CLOSE *********************************************
C            -    ----- *********************************************
C Steve Gibbons 14.4.97                                             C
C  ( note that this is essentially the routine of Dan Gordon        C
C___________________________________________________________________C
C Closes file with integer logical unit LU, filename FNAME.         C
C LABEL contains any other information regarding the nature of the  C
C the file.							    C
C___________________________________________________________________C
C Input Variables :-						    C
C ===============						    C
C  Integer							    C
C  -------							    C
C     LU	: Number of file.				    C
C  Character							    C
C  ---------							    C
C     FNAME	: Name of file. Undefined length		    C
C     LABEL	: Any further information. Undefined length         C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE FCLOSE ( LU, FNAME, LABEL )
      IMPLICIT NONE
C___________________________________________________________________C
C Variable Declarations - Parameters ...............................C
      INTEGER LU
      CHARACTER *(*) FNAME
      CHARACTER *(*) LABEL
C___________________________________________________________________C
C START OF PROGRAM *************************************************C
C___________________________________________________________________C
      IF ( LU.EQ.0 ) THEN
         PRINT *,' Subroutine FCLOSE '
         PRINT *,' I bet you ve forgotten to set LU ...??'
         PRINT *,' Think again and come back when you have'
         PRINT *,' remembered that LU must be a non zero integer!'
         PRINT *,' See you later. Bye for now!!'
         STOP
      ENDIF
C----------------------------
      CLOSE (UNIT=LU, STATUS='KEEP', ERR=989 )
      RETURN
C
 989  PRINT *,' Error.  Failed to close ', LABEL, ' file ', FNAME
      STOP
      END
C********************************************************************
C*********************************************************************
C subroutine Boundary Independent Harmonic File ReaD *****************
C            -        -           -        -    -  - *****************
C Steve Gibbons Wed Dec  8 09:58:50 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Reads in the integer indices of spherical harmonic sets BUT NOT    C
C the appropriate boundary conditions from a file.                   C
C This is a simpler routine than HMFRD and avoids the need to set    C
C up additional arrays if knowledge of boundary conditions are not   C
C necessary.                                                         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Number of vector spherical harmonics. (Output)     C
C     NHMAX     : Maximum number of vector spherical harmonics.      C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NHMAX          C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NHMAX          C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NHMAX          C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE BIHFRD( NH, NHMAX, MHT, MHL, MHM, LU, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, NHMAX, MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ), LU
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH, IWR, IIBF, IOBF
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
C Read in number of spherical harmonics
C
      READ ( LU, * ) NH
C
C Check value of NH
C
      IF ( NH.GT.NHMAX ) THEN
        PRINT *, ' Subroutine BIHFRD.'
        PRINT *, ' From file, NH = ', NH
        PRINT *, ' NHMAX = ', NHMAX
        PRINT *, ' Program aborted.'
        STOP
      ENDIF
C
C     Now loop around the harmonics and read in their
C     properties
C
      DO IH = 1, NH
        READ ( LU, * ) MHT( IH ), MHL( IH ), MHM( IH ), IIBF, IOBF
      ENDDO
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine Imposed Temperature Vector Auxiliary Arrays Form ********
C            -       -           -      -         -      -    ********
C Steve Gibbons Fri Dec 15 10:22:21 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C There are NH3 temperature harmonics in the time-stepping solution  C
C vector. They are indexed by the integer arrays ML3 and MM3.        C
C Each are represented at NR radial grid points and the element for  C
C node ir and harmonic ih is stored in SV3( ( IH - 1 )*NR + IR ).    C
C The radial spacings are given by the array XARR.                   C
C                                                                    C
C Now the "foreign" temperature is stored in SVC with indexing       C
C by INARRC, MHTC, MHLC, MHMC and with grid nodes stored in the      C
C array XARRC. Now XARRC and XARR may be different but the first and C
C last points must really be the same! So, we have a variable called C
C DTOL which checks how close the inner and outer values are         C
C between the two arrays. If ( DABS( RI - RIC ).GT.DTOL   ) or       C
C ( DABS( RO - ROC ).GT.DTOL   ) then the program will abort.        C
C                                                                    C
C The routine uses the subroutine SVRINT to interpolate the 0th (F), C
C 1st (DF) and 2nd (DDF) derivatives of SVC and put into the         C
C appropriate positions of the arrays DV3C, SV3C and DSV3C           C
C respectively                                                       C
C                                                                    C
C SCAL*CD*DL( L , R , F , DF , DDF ),                                C
C SCAL*F                                                             C
C SCAL*DF.                                                           C
C                                                                    C
C DV3C, SV3C and DSV3C are blanked upon entry.                       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH3       : Number of temperature harmonics in main soln.      C
C     NR        : Number of radial grid nodes in main soln.          C
C     ML3       : Dim ( NH3 ). Sph. harm. degrees of temp. harms.    C
C     MM3       : Dim ( NH3 ). Sph. harm. orders of temp. harms.     C
C                                                                    C
C     INARRC    : Indexing array (dim 3) for foreign temperature.    C
C                   INARRC( 1 ) = iformfc                            C
C                   INARRC( 2 ) = NRC                                C
C                   INARRC( 3 ) = NHC                                C
C     MHTC      : Type array (dim NHC) for foreign temperature.      C
C     MHLC      : l array (dim NHC) for foreign temperature.         C
C     MHMC      : m array (dim NHC) for foreign temperature.         C
C                                                                    C
C     NNDS      : Number of nodes for interpolation. (atleast 3 )    C
C     IWORK     : Dimension ( NNDS ). Work array.                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim (NR). Radial grid nodes of main solution.      C
C     XARRC     : Dim (NRC). Radial grid nodes of foreign temp.      C
C     SVC       : Dim (NRC*NHC). foreign temp. solution vector.      C
C                                                                    C
C     CD        : Coefficient of diffusion for temperature.          C
C     SCAL      : Multiplier for inhomogeneous part of the temp.     C
C                                                                    C
C     DTOL      : Tolerance allowed for difference of inner and      C
C                 outer boundaries between the two spacings arrays.  C
C                                                                    C
C     DV3C      : Diffusion forcing term from inhomog. temp (NR*NH3) C
C     SV3C      : Zero^th derivative interpolation      Dim (NR*NH3) C
C     DSV3C     : First derivative interpolation      Dim (NR*NH3)   C
C                                                                    C
C     WORK1     : Dimension ( NNDS ). Work array.                    C
C     WORK2     : Dimension ( NNDS ). Work array.                    C
C     WORKM     : Dimension ( NNDS, NNDS ). Work array.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITVAAF( NH3, NR, ML3, MM3, INARRC, MHTC, MHLC, MHMC,
     1                  NNDS, IWORK, XARR, XARRC, SVC, CD, SCAL,
     2                  DTOL, DV3C, SV3C, DSV3C, WORK1, WORK2, WORKM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NH3, NR, ML3( NH3 ), MM3( NH3 ), INARRC( 3 ),
     1                 MHTC( * ), MHLC( * ), MHMC( * ), NNDS,
     2                 IWORK( NNDS )
      DOUBLE PRECISION XARR( NR ), XARRC( * ), SVC( * ), CD, SCAL,
     1                 DTOL, DV3C( NR*NH3 ), SV3C( NR*NH3 ),
     2                 DSV3C( NR*NH3 ), WORK1( NNDS ), WORK2( NNDS ),
     3                 WORKM( NNDS, NNDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IR, IH, IHC, IND, NRC, NHC, IOP, IHARM, L, N
      DOUBLE PRECISION RAD, RI, RIC, RO, ROC, ZERO, F, DF, DDF, DL
      EXTERNAL         DL
      PARAMETER      ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Zero arrays.
C
      IOP     = 0
      N       = NR*NH3
      CALL VECOP(  DV3C, ZERO, N, IOP )
      CALL VECOP(  SV3C, ZERO, N, IOP )
      CALL VECOP( DSV3C, ZERO, N, IOP )
C
      NRC     = INARRC( 2 )
      NHC     = INARRC( 3 )
C
      RI      = XARR(   1 )
      RO      = XARR(  NR )
C
      RIC     = XARRC(   1 )
      ROC     = XARRC( NRC )
C
      IF ( DABS( RI-RIC ).GT.DTOL  .OR.
     1     DABS( RO-ROC ).GT.DTOL      ) THEN
        PRINT *,' Subroutine ITVAAF.'
        PRINT *,' RI = ', RI,' RIC = ', RIC
        PRINT *,' RO = ', RO,' ROC = ', ROC
        PRINT *,' Incompatible radial functions.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO IHC = 1, NHC
        IF ( MHTC( IHC ).NE.3 ) GOTO 50
C       .
C       . OK so IHC is a temperature harmonic -
C       . so let's try to find the corresponding harmonic.
C       .
        L     = MHLC( IHC )
        IHARM = 0
        DO IH = 1, NH3
          IF (     ML3( IH ).EQ.MHLC( IHC )   .AND.
     1             MM3( IH ).EQ.MHMC( IHC )        ) IHARM = IH
        ENDDO
        IF ( IHARM.EQ.0 ) THEN
          PRINT *,' Subroutine ITVAAF.'
          PRINT *,' Solution vector contains no harmonic, ih,'
          PRINT *,' with ML3( ih ) = ', MHLC( IHC )
          PRINT *,' and MM3( ih ) = ', MHMC( IHC )
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
        DO IR = 2, NR - 1
          IND    = ( IHARM - 1 )*NR + IR
          RAD    = XARR( IR )
          CALL SVRINT( RAD, SVC, XARRC, INARRC, IHC, NNDS, WORK1,
     1                   IWORK, WORK2, WORKM )
C         .
          F      = WORK1( 1 )
          DF     = WORK1( 2 )
          DDF    = WORK1( 3 )
C         .
          DV3C ( IND ) = SCAL*CD*DL( L , RAD , F , DF , DDF )
          SV3C ( IND ) = SCAL*F
          DSV3C( IND ) = SCAL*DF
C         .
        ENDDO
C       .
 50   CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine File NAME giveR *****************************************
C            -    ----     - *****************************************
C Steve Gibbons 14.4.97 (Adapted from Dan Gordon's course.)          C
C____________________________________________________________________C
C Asks the user for a file name FNAME. Pretty simple really ...      C
C____________________________________________________________________C
C Input Variable :-						     C
C ==============   						     C
C  Character							     C
C  ---------							     C
C     LABEL	: Message arbitrary length  			     C
C____________________________________________________________________C
C Output Variable :-						     C
C ===============   						     C
C  Character							     C
C  ---------							     C
C     FNAME	: Filename arbitrary length			     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FNAMER ( FNAME, LABEL )
      IMPLICIT NONE
      CHARACTER *(*) FNAME
      CHARACTER *(*) LABEL

      PRINT *,' Please enter a filename for FNAME.'
      PRINT *, LABEL
      READ (5, 200) FNAME
 200  FORMAT (A)
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Adapted Solution Vector DeRivative **********************
C            -       -        -      - -        **********************
C Steve Gibbons Mon Oct 25 15:17:59 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C V is a solution vector with NH (= inarr(3) ) harmonic radial       C
C functions and NR (= inarr(2) ) grid nodes in each function.        C
C The position of the j^{th} node of the i^{th} harmonic is given by C
C INDFUN( j, i, INARR) and the radial value is given by              C
C  XARR( j ) - an array which is passed to the routine FDCMBD        C
C in order to calculate the finite difference coefficients, FDCM.    C
C XARR itself is not referenced by ASVDR.                            C
C ASVDR returns the radial derivatives 0, 1, ..., IHD of radial      C
C function IH evaluated at node IR.                                  C
C                                                                    C
C IFORMF = INARR(1) should be either 3 or 4 since this is the        C
C arbitrarily spaced mesh version of the code.                       C
C                                                                    C
C NBN is the maximum number of nodes on                              C
C either side which maybe used in central differences.               C
C For instance, if to calculate the derivative of                    C
C f at r_j, you may use the values of f at r = r_{j-2}, r_{j-1},     C
C r_j, r_{j+1} and r_{j+2} then NBN = 2. The value of IHD is checked C
C only for being positive and no greater than NDRVS (the             C
C number of the highest derivative for which coefficients            C
C are stored by the array SVFDC), as SVFDC must be calculated in     C
C advance by a call to SVFDCF which checks NDRVS against             C
C the physical restrictions imposed by the value of NBN.             C
C                                                                    C
C NDRVM restricts the size of NDRVS and is a defining parameter      C
C of the array SVFDCF.                                               C
C                                                                    C
C NBN must be as supplied to SVFDCF.                                 C
C                                                                    C
C IS must be supplied to indicate the finite difference scheme       C
C being employed.                                                    C
C                                                                    C
C ASVDR will use whichever coefficients SVFDC( k, j, nd + 1, K ),    C
C that SVFDCF formed with K = IS.                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Number of radial grid node.                        C
C     IS        : Number of finite difference scheme used.           C
C     IH        : Number of radial function (harmonic).              C
C     NBN       : Number of bounding nodes. See above.               C
C     IHD       : Highest derivative requested.                      C
C     NFDCM     : Leading dimension of SVFDC. At least (2*NBN+1)     C
C     NR        : Number of radial grid nodes in each function.      C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NRR. Must be consistent with NR.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V         : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C     DERV      : Derivatives. Dim ( * ) but length atleast IHD      C
C                  DERV( i ) is returned containing the value of     C
C                  the i^[th} derivative.                            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ASVDR ( V, IR, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                   NDRVM, DERV, INARR, SVFDC, NDCS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, IS, IH, NBN, IHD, NFDCM, NR, NDRVS, NDRVM,
     1        INARR( * ), NDCS
      DOUBLE PRECISION V( * ), DERV( * ), 
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION COEF
      INTEGER INODE, ILN, IRN, INDFUN, ID, IND, NRR, IFORMF, IK, ID1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameters ...
C
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
C
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine ASVDR.'
         PRINT *,' INARR( 2 ) = ', NRR
         PRINT *,' NR = ', NR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
         PRINT *,' Subroutine ASVDR.'
         PRINT *,' INARR( 1 ) = ', IFORMF
         PRINT *,' This is an irregular grid routine.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IR.LT.1 .OR. IR.GT.NR ) THEN
        PRINT *,' Subroutine ASVDR.'
        PRINT *,' IR   = ', IR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IHD.LT.0 .OR. IHD.GT.NDRVS .OR. NDRVS.GT.NDRVM ) THEN
        PRINT *,' Subroutine ASVDR.'
        PRINT *,' IHD   = ', IHD
        PRINT *,' NDRVS = ', NDRVS
        PRINT *,' NDRVM = ', NDRVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IS.LT.1 .OR. IS.GT.NDCS ) THEN
        PRINT *,' Subroutine ASVDR.'
        PRINT *,' IS    = ', IS
        PRINT *,' NDCS  = ', NDCS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Calculate furthest left and furthest right mode
C     . to be used to form derivative.
C     .
      ILN = MAX(  1, IR - NBN )
      IRN = MIN( NR, IR + NBN )
C
      DO ID = 0, IHD
        ID1 = ID + 1
        DERV( ID1 ) = 0.0d0
        DO INODE = ILN, IRN
          IK = INODE - IR + NBN + 1
          COEF = SVFDC( IK, IR, ID1, IS )
          IND = INDFUN( INODE, IH, INARR )
          DERV( ID1 ) = DERV( ID1 ) + COEF*V( IND )
        ENDDO
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine Linear Dependence of Grid Node Matrix Form **************
C            -      -             -    -    -      -    **************
C Steve Gibbons Sat Oct 23 15:01:52 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Let f be a function of x and let f_j denote the value of f at      C
C the grid node j (x value is x_j given by XARR( j ) ... )           C
C                                                                    C
C If f has to satisfy a particular boundary condition then           C
C all the f_j may not be linearly independent.                       C
C                                                                    C
C For a group of n ( = NNDS ) grid nodes, from s = k+1 to k+n        C
C for some integer, k,                                               C
C                                                                    C
C  f_{k+i} = \sum_{s=1}^n DMAT( i, s ) f_{k+s},                      C
C                                                                    C
C and the routine LDGNMF returns the matrix DMAT.                    C
C                                                                    C
C The number of linearly dependent nodes to the left and right are   C
C given by NALF and NARF respectively.                               C
C                                                                    C
C The boundary conditions at the inner and outer boundaries are      C
C specified by the integers IIBC and IOBC which may take the         C
C following values                                                   C
C                                                                    C
C    IIBC              Inner Boundary Condition                      C
C    ====              ========================                      C
C                                                                    C
C      1       No condition to be applied                            C
C      2       Function must vanish                                  C
C      3       First derivative must vanish                          C
C      4       Function AND first derivative must vanish             C
C      5       Function AND second derivative must vanish            C
C      6       rdf/dr - f(r) = 0                                     C
C      7       r df/dr - l f(r) = 0                                  C
C                                                                    C
C    IOBC              Outer Boundary Condition                      C
C    ====              ========================                      C
C                                                                    C
C      1       No condition to be applied                            C
C      2       Function must vanish                                  C
C      3       First derivative must vanish                          C
C      4       Function AND first derivative must vanish             C
C      5       Function AND second derivative must vanish            C
C      6       rdf/dr - f(r) = 0                                     C
C      7       r df/dr + (l+1) f(r) = 0                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Total number of radial grid nodes.                 C
C     NNDS      : Number of nodes needed to take derivative.         C
C     NALF      : Number of nodes to left which are a linear         C
C                 combination of the other nodes.                    C
C     NARF      : Number of nodes to right which are a linear        C
C                 combination of the other nodes.                    C
C     L         : Spherical harmonic degree, l.                      C
C     IIBC      : Inner boundary flag - see above.                   C
C     IOBC      : Outer boundary flag - see above.                   C
C     NCFM      : Leading dimension of array DMAT etc.               C
C     IPCM      : Dimension ( NCFM ). Working array.                 C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension ( NR ).                         C
C                  XARR( i ) contains the value of r at the          C
C                   i^{th} grid node.                                C
C                                                                    C
C     DMAT      : Dimension ( NCFM, NCFM ). See above.               C
C     WMAT      : Dimension ( NCFM, NCFM ). Working array.           C
C     WORK1     : Dimension ( NCFM ). Working array.                 C
C     WORK2     : Dimension ( NCFM ). Working array.                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LDGNMF( NR, NNDS, NALF, NARF, L, IIBC, IOBC, NCFM,
     1                   XARR, DMAT, WMAT, WORK1, WORK2, IPCM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NNDS, NALF, NARF, L, IIBC, IOBC, NCFM,
     1        IPCM( NCFM )
      DOUBLE PRECISION XARR( NR ), DMAT( NCFM, NCFM),
     1                 WMAT( NCFM, NCFM), WORK1( NCFM ),
     2                 WORK2( NCFM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER I, IFN, ILN, NDS2, ID1, IHND
      DOUBLE PRECISION ZERO, TOL, X0, FAC, QUOT
      PARAMETER ( ZERO = 0.0d0, TOL = 1.0d-9 )
C ifn is first linearly independent node
C iln is last linearly independent node
C ihnd is the highest number derivative which
C will be needed to be calculated during this
C routine (by GFDCFD)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First check on values of NALF and NARF
C
      I = NALF*NARF
      IF ( I.NE.0 .OR. (NALF.LT.0) .OR. (NALF.GT.2) .OR.
     1                 (NARF.LT.0) .OR. (NARF.GT.2)       ) THEN
        PRINT *,' Subroutine LDGNMF.'
        PRINT *,' NALF = ', NALF
        PRINT *,' NARF = ', NARF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Ok - the number of affected nodes is valid
C     . so let's proceed. First, zero the matrix DMAT
C     .
      I = 0
      CALL MATOP( DMAT, ZERO, NCFM, NCFM, I )
C     .
C     . Now add the diagonal elements for the
C     . linearly independent variables
C     .
      IF ( NALF.EQ.0 ) IFN = 1
      IF ( NALF.EQ.1 ) IFN = 2
      IF ( NALF.EQ.2 ) IFN = 3
C     .
      IF ( NARF.EQ.0 ) ILN = NNDS
      IF ( NARF.EQ.1 ) ILN = NNDS - 1
      IF ( NARF.EQ.2 ) ILN = NNDS - 2
C     .
      DO I = IFN, ILN
        DMAT( I, I ) = 1.0d0
      ENDDO
C     .
C     . We can now return if the identity matrix is required
C     .
      IF ( NALF.EQ.0 .AND. NARF.EQ.0 ) RETURN
C     .
C     .
C     . OK - now decide whether we are doing inner or
C     . outer boundary
C     .
      IF ( NARF.EQ.0 ) THEN
C       .
C       . We are considering the inner boundary
C       . First, return if IIBC = 1 or 2
C       .
        IF ( IIBC.EQ.1 .OR. IIBC.EQ.2 ) RETURN
C       .
C       . So we need to calculate deriv.s at
C       . the inner boundary, XARR( 1 ).
C       . NDS2 is the number of nodes we can use at
C       . the inner boundary.
C       .
        IF ( IIBC.EQ.5 ) THEN
          IHND = 2
        ELSE
          IHND = 1
        ENDIF
C       .
        NDS2 = 1 + NNDS - NALF
C       .
        IF ( (NDS2-1).LT.IHND ) THEN
           PRINT *,' Subroutine LDGNMF '
           PRINT *,' NDS2                = ', NDS2
           PRINT *,' Required derivative = ', IHND
           PRINT *,' Program aborted.'
           STOP
        ENDIF
C       .
        DO I = 1, NDS2
          WORK1( I ) = XARR( I )
        ENDDO
        X0 = XARR( 1 )
        CALL GFDCFD ( X0, WORK1, NDS2, WMAT, NCFM, IPCM, WORK2 )
C       .
C       . WMAT( m + 1, i ) now contains the coefficient
C       . with which you multiply f_i to get the m^{th}
C       . derivative of f at x_0.
C       .
C       . Consider the case IIBC = 3
C       . We require the first derivative = 0
C       .
        IF ( IIBC.EQ.3 ) THEN
C          .
C          . Just quickly check that NALF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NALF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IIBC = ',IIBC,' and NALF = ',NALF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . Check for imminent division by zero
C          .
           IF ( ABS( WMAT( 2, 1 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT( 2, 1 ) = ',WMAT( 2, 1 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 2, NDS2
             DMAT( 1, I ) = (-1.0d0)*WMAT( 2, I )/WMAT( 2, 1 )
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . Consider the case IIBC = 4(5)
C       . We need both the function and the first
C       . (second) derivative to vanish. We do not need to
C       . to change anything to make the function
C       . vanish so just need do the first (second) deriv. cond.
C       .
        IF ( IIBC.EQ.4 .OR. IIBC.EQ.5 ) THEN
C          .
C          . Check for imminent division by zero
C          .
           IF ( IIBC.EQ.4 ) ID1 = 2
           IF ( IIBC.EQ.5 ) ID1 = 3
C          .
           IF ( ABS( WMAT( ID1, 2 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT(',ID1,', 2 ) = ',WMAT( ID1, 2 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 3, NDS2
             DMAT( NALF, NALF - 2 + I ) =
     1                    (-1.0d0)*WMAT( ID1, I )/WMAT( ID1, 2 )
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . We may treat the cases IIBC.EQ.6 and IIBC.EQ.7
C       . together as they are both conditions of the form
C       . rdf/dr + FAC f(r) = 0
C       .
        IF ( IIBC.EQ.6 .OR. IIBC.EQ.7 ) THEN
C          .
C          . First we make an early escape if r_{inner bnd} = 0
C          . this is equivlent to setting the function to zero
C          . which is what the current status of the matrix is
C          .
           IF ( ABS( X0 ).LT.TOL ) RETURN
C          .
C          . OK - so it's non-trivial!
C          . First, let's allocate the correct value of FAC
C          .
           IF ( IIBC.EQ.6 ) FAC = -1.0d0
           IF ( IIBC.EQ.7 ) FAC = -1.0d0*DBLE( L )
C          .
C          . Just quickly check that NALF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NALF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IIBC = ',IIBC,' and NALF = ',NALF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . Check for imminent division by zero
C          . We can safely divide by X0 now since the
C          . X0 = 0.0 case has been discounted above.
C          .
           QUOT = FAC/X0 + WMAT( 2, 1 )
           IF ( ABS( QUOT ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' QUOT = ', QUOT
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 2, NDS2
             DMAT( 1, I ) = (-1.0d0)*WMAT( 2, I )/QUOT
           ENDDO
           RETURN
C          .
        ENDIF
C       .
      ELSE
C       .
C       . We are considering the outer boundary
C       . First, return if IOBC = 1 or 2
C       .
        IF ( IOBC.EQ.1 .OR. IOBC.EQ.2 ) RETURN
C       .
C       . So we need to calculate deriv.s at
C       . the outer boundary, XARR( NR ).
C       . NDS2 is the number of nodes we can use at
C       . the outer boundary.
C       .
        IF ( IOBC.EQ.5 ) THEN
          IHND = 2
        ELSE
          IHND = 1
        ENDIF
C       .
        NDS2 = 1 + NNDS - NARF
C       .
        IF ( (NDS2-1).LT.IHND ) THEN
           PRINT *,' Subroutine LDGNMF '
           PRINT *,' NDS2                = ', NDS2
           PRINT *,' Required derivative = ', IHND
           PRINT *,' Program aborted.'
           STOP
        ENDIF
C       .
        DO I = 1, NDS2
          WORK1( I ) = XARR( NR - NDS2 + I )
        ENDDO
        X0 = XARR( NR )
        CALL GFDCFD ( X0, WORK1, NDS2, WMAT, NCFM, IPCM, WORK2 )
C       .
C       . WMAT( m + 1, i ) now contains the coefficient
C       . with which you multiply f_{nr-nds2+i} to get the m^{th}
C       . derivative of f at x_0.
C       .
C       . Consider the case IOBC = 3
C       . We require the first derivative = 0
C       .
        IF ( IOBC.EQ.3 ) THEN
C          .
C          . Just quickly check that NARF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NARF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IOBC = ',IOBC,' and NARF = ',NARF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . Check for imminent division by zero
C          .
           IF ( ABS( WMAT( 2, NDS2 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT( 2,',NDS2,') = ',WMAT( 2, NDS2 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          .
           DO I = 1, NDS2-1
             DMAT( NNDS, NNDS - NDS2 + I ) =
     1             (-1.0d0)*WMAT( 2, I )/WMAT( 2, NDS2 )
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . Consider the case IOBC = 4(5)
C       . We need both the function and the first
C       . (second) derivative to vanish. We do not need to
C       . to change anything to make the function
C       . vanish so just need do the first (second) deriv. cond.
C       .
        IF ( IOBC.EQ.4 .OR. IOBC.EQ.5 ) THEN
C          .
C          . Check for imminent division by zero
C          .
           IF ( IOBC.EQ.4 ) ID1 = 2
           IF ( IOBC.EQ.5 ) ID1 = 3
C          .
           IF ( ABS( WMAT( 2, NDS2-1 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT(',ID1,',',NDS2-1,') = ',WMAT( ID1,NDS2-1)
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 1, NDS2-2
             DMAT( NNDS + 1 - NARF, NNDS - NDS2 - NARF + 2 + I ) =
     1                    (-1.0d0)*WMAT( ID1, I )/WMAT( ID1, NDS2-1)
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . We may treat the cases IOBC.EQ.6 and IOBC.EQ.7
C       . together as they are both conditions of the form
C       . rdf/dr + FAC f(r) = 0
C       .
        IF ( IOBC.EQ.6 .OR. IOBC.EQ.7 ) THEN
C          .
C          . Just quickly check that NARF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NARF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IOBC = ',IOBC,' and NARF = ',NARF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . First, let's allocate the correct value of FAC
C          .
           IF ( IOBC.EQ.6 ) FAC = -1.0d0
           IF ( IOBC.EQ.7 ) FAC = DBLE( L + 1 )
C          .
C          . Check for imminent division by zero
C          . (Once again, we safely divide by X0)
C          .
           QUOT = FAC/X0 + WMAT( 2, NDS2 )
           IF ( ABS( QUOT ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' QUOT = ', QUOT
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          .
           DO I = 1, NDS2-1
             DMAT( NNDS, NNDS - NDS2 + I ) =
     1             (-1.0d0)*WMAT( 2, I )/QUOT
           ENDDO
           RETURN
C          .
        ENDIF
C       .
      ENDIF
C     .
      PRINT *,' Subroutine LDGNMF '
      PRINT *,' Problem with inputs.'
      PRINT *,' IIBC = ', IIBC
      PRINT *,' IOBC = ', IOBC
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
C*********************************************************************
C subroutine General Finite Difference Coefficient Find **************
C            -       -      -          -           -    **************
C Steve Gibbons Mon Sep 20 16:57:54 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Given a value of X and the values of x_i at NNDS distinct points,  C
C (X need not necessarily be one of the x_i) then the array COEFM    C
C is returned with the finite difference coefficients such that      C
C the ND^{th} derivative of a function f, evaluated at x = X,        C
C is given by                                                        C
C                                                                    C
C f^{ ND }( X ) = \sum_{j = 1, NNDS} COEFM( ND + 1, j )*f( x_j )     C
C                                                                    C
C This is a general version of the routine FDCINV which is valid     C
C only for equally spaced grid nodes.                                C
C                                                                    C
C Coefficients for up to the (NNDS-1)^{th} derivative are given      C
C although care must be taken to ensure the highest derivatives      C
C are sufficiently accurate.                                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Abscissa at which derivatives are to be            C
C                  evaluated.                                        C
C     XARR      : Array of dimension ( NNDS ).                       C
C                  XARR( i ) contains the value of x/r at the        C
C                   i^{th} grid node.                                C
C                                                                    C
C     COEFM     : Dimension ( NCFM, NCFM).                           C
C     WORK      : Workspace array for LAPACK inversion routine.      C
C                 Dimension ( NCFM )                                 C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NNDS      : Number of grid nodes.                              C
C     NCFM      : Leading order of coefficient matrix.               C
C                                                                    C
C     IPCM      : Work array for LAPACK routines to perform          C
C                 pivotting in the matrix inversion.                 C
C                 Dimension ( NCFM )                                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE GFDCFD ( X, XARR, NNDS, COEFM, NCFM, IPCM, WORK)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NNDS, NCFM, IPCM( NCFM )
      DOUBLE PRECISION X, COEFM( NCFM, NCFM ), WORK( NCFM ),
     1                 XARR( NNDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER I, J, NDER, INODE, INFO, ICOL, IROW
      DOUBLE PRECISION DZERO, LOW, FAC
      PARAMETER ( DZERO = 0.0d0, LOW = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check bounds of integer parameters
C
      IF ( NNDS.GT.NCFM ) THEN
         PRINT *,' Subroutine GFDCFD: '
         PRINT *,' NNDS = ', NNDS,'. NCFM = ', NCFM
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C Calculate the distances h_i = ( x_i - X )
C These h_i must be distinct and this will be
C checked for - otherwise matrix is singular.
C We can store the h_i in WORK as this will not
C be needed until the inversion, by which time
C it will not be needed by us.
C
      DO I = 1, NNDS
        WORK( I ) = XARR( I ) - X
      ENDDO
C
C Now check for the uniqueness of the points ...
C
      DO I = 1, NNDS - 1
        DO J = I + 1, NNDS
          IF ( ABS( WORK( I ) - WORK( J ) ).LT.LOW ) THEN
            PRINT *,' Subroutine GFDCFD.'
            PRINT *,' X values ',I,' and ',J,' are'
            PRINT *,' identical.'
            PRINT *,' Program aborted.'
            STOP
          ENDIF
        ENDDO
      ENDDO
C
C____________________________________________________________________C
C Parameters are ok so let's zero COEFM
C
      I = 0
      CALL MATOP( COEFM, DZERO, NCFM, NCFM, I )
C
C (nder+1) is the number of the matrix column being filled in.
C inode is the number of the matrix row being filled in.
C
      DO NDER = 0, NNDS - 1
        ICOL = NDER + 1
        DO INODE = 1, NNDS
          IROW = INODE
          IF ( NDER.EQ.0 ) THEN
            COEFM( IROW, ICOL )  = 1.0d0
          ELSE
            FAC = WORK( INODE )/DBLE( NDER )
            COEFM( IROW, ICOL ) = COEFM( IROW, ICOL-1 )*FAC
          ENDIF
        ENDDO
      ENDDO
C
C Ok - this matrix is now ready for inversion -
C For this we use the LAPACK routines DGETRF and DGETRI
C First perform LU decomposition
C
      CALL DGETRF( NNDS, NNDS, COEFM, NCFM, IPCM, INFO )
C
C     . Check that LU decomposition has gone without
C     . problem.
C     .
C
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine GFDCFD.'
         PRINT *,' The LAPACK subroutine DGETRF has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in LU decomposition of COEFM matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Now compute the inverse with the LAPACK routine
C     . DGETRI.
C     .
      CALL DGETRI( NNDS, COEFM, NCFM, IPCM, WORK, NCFM, INFO )
C     .
C     . Check that inversion has gone without problem.
C     .
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine GFDCFD.'
         PRINT *,' The LAPACK subroutine DGETRI has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in inversion of COEFM matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine MATrix OPeration ****************************************
C Steve Gibbons 23.4.97 Does operation on a two-dimensional array.   C
C                                                Can set equal to a  C
C                       constant; multiply by a constant or have a   C
C                       constant added to it.                        C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IOP	: Type of operation to be done.                      C
C                    WHOLE MATRIX OPERATIONS                         C
C                  IOP=0  -->  Each element of the matrix = CONST    C
C                  IOP=1  -->  Each el. is multiplied by CONST       C
C                  IOP=2  -->  Each el. is added to CONST            C
C     NDIM1     : First dimension of the matrix.	             C
C     NDIM2     : Second dimension of the matrix.	             C
C  Double Precision                                                  C
C  ----------------                                                  C
C     MAT	: Matrix with dimension ( NDIM1, NDIM2 )             C
C     CONST     : Double precision constant.                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MATOP ( MAT, CONST, NDIM1, NDIM2, IOP)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM1, NDIM2, IOP
      DOUBLE PRECISION MAT ( NDIM1, NDIM2 ), CONST
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, J
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First do case IOP=0 
      IF ( IOP.EQ.0 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               MAT ( I, J) = CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
C Now do case IOP=1
      IF ( IOP.EQ.1 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               MAT ( I, J) = MAT ( I, J)*CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
C Now do case IOP=2
      IF ( IOP.EQ.2 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               MAT ( I, J) = MAT ( I, J) + CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
      PRINT *,' Subroutine MATOP. IOP must be 0,1 or 2.'
      STOP
      END
C*********************************************************************
C*********************************************************************
C subroutine Adapted Matrix Linear Interaction Contribution Add ******
C            -       -      -      -           -            -   ******
C Steve Gibbons Tue Oct 26 15:58:14 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C  Adds to the double precision matrix, A, the terms which a single  C
C harmonic IHC contributes to a single harmonic IHR in a single term C
C of a coupled o.d.e. when the grid in radius is not necessarily     C
C uniform and the j^[th} abscissa has a radial value of XARR( j )    C
C                                                                    C
C A Linear contribution indicates one which is merely the function   C
C of radius and a finite number of integer and double precision      C
C constants.                                                         C
C                                                                    C
C See below for the form of the subroutine which gives this          C
C function.                                                          C
C                                                                    C
C Let the radial function of the IHR harmonic be denoted f_{IHR)(r)  C
C Let the radial function of the IHC harmonic be denoted f_{IHC)(r)  C
C                                                                    C
C If the contribution for a given radial node, $r_j$, is of the form C
C                                                                    C
C f_{IHR)(r_j) = \sum_{nd=0,IHD} f_{IHC}^{nd}(r_j) c_{IHC}^{nd}(r_j) C
C                                                                    C
C  (where f^{nd} denotes the (nd)^{th} derivative of f with respect  C
C  to r, evaluated at r_j)                                           C
C                                                                    C
C then the coefficient c_{IHC}^{nd}(r_j) must be supplied by the     C
C subroutine SUB1 (declared EXTERNAL in the calling (sub)program )   C
C which returns c_{IHC}^{nd}(r_j) in the array element CVEC(nd + 1). C
C                                                                    C
C The subroutine SUB1 *MUST* have the calling sequence               C
C                                                                    C
C SUB1( CVEC, RAD, IPARS, DPARS, IHD )                               C
C                                                                    C
C where  RAD is the double precision value of the radius,            C
C        IPARS is an integer array to provide SUB1 with parameters.  C
C          IPARS is not referred to by AMLICA other than to pass     C
C          this information.                                         C
C        DPARS is like IPARS but contains double precision elements. C
C        IHD is the highest derivative which will be refered to by   C
C        AMLICA. So all CVEC( i ) must be zero with 1.le.i.le.ihd+1  C
C        if it is not assigned another value.                        C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C                                                                    C
C     IMF       : Matrix format flag.                                C
C                                                                    C
C          imf = 1; Matrix is in LAPACK banded format                C
C                   ie element a_{i,j} is stored in                  C
C                   A( kle + ku + 1 + i - j , j )                    C
C                                                                    C
C          imf = 2; Matrix is banded but with element a_{i,j}        C
C                   stored in A( kl + 1 + j - i , i ).               C
C                                                                    C
C          imf = 3; Matrix is square - ie a_{i,j} is stored          C
C                   in A( i, j ).                                    C
C                                                                    C
C     IHC       : Number of the input (column) harmonic              C
C     IHR       : Number of the output (row) harmonic                C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1. INDFUN = ( IR - 1 )*NH + IH          C
C                   IFORMF = 2. INDFUN = ( IH - 1 )*NR + IR          C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     IHD       : Highest derivative involved.                       C
C     NBN       : Number of diagonal elements in radius.             C
C     ILNR      : Left-most node at which to start adding rows to    C
C                   matrix.                                          C
C     IRNR      : Right-most node at which to start adding rows to   C
C                   matrix.                                          C
C                                                                    C
C     NFDCM     : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     IS        : Number of required finite difference scheme.       C
C                  (This points to the fourth dim. element in        C
C                   the coefficient array SVFDC.)                    C
C                                                                    C
C     NDRVS     : Highest derivative stored in SVFDC.                C
C                 (Must be atleast 1).                               C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     IPARS     : Information array for SUB1. Dimension ( * )        C
C                                                                    C
C  Subroutines                                                       C
C  -----------                                                       C
C     SUB1      : Determines what multiplies each derivative in      C
C                  the matrix. Must have calling sequence ...        C
C                                                                    C
C     CALL SUB1 ( CVEC, RAD, IPARS, DPARS, IHD )                     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     AMAT      : Matrix. Dimensions ( N1, N2 )                      C
C              Will generally be banded due to the nature of the     C
C              numerical scheme. KL, KU and KLE parameterise this.   C
C                                                                    C
C     FAC       : Multiplication factor.                             C
C                                                                    C
C     XARR      : Array of dimension ( NR )                          C
C                 XARR( j ) = element x_j                            C
C                                                                    C
C     WORK      : Working array of dimension ( NDRVS + 1 )           C
C                                                                    C
C     DPARS     : Information array for SUB1. Dimension ( * )        C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AMLICA( N1, N2, KL, KU, KLE, IMF, IHC, IHR, INARR,
     1                   IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2                   NDRVS, NDRVM, IPARS, SUB1, AMAT, FAC, XARR,
     3                   WORK, DPARS, SVFDC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, IMF, IHC, IHR, INARR( * ), IHD,
     1        NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS, NDRVS, NDRVM,
     2        IPARS( * )
      EXTERNAL SUB1
      DOUBLE PRECISION AMAT( N1, N2 ), FAC, XARR( NR ),
     1                 WORK( NDRVS + 1 ), DPARS( * ),
     2                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IRAD, INODE, NLN, NRN, INDFUN, INDR, INDC,
     1        IROW, ICOL, NH, IFORMF, NLCS, NRCS, I,
     2        NDER, ICOROW, ICOCOL, NRR
      DOUBLE PRECISION RAD, LOW
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First we check for an easy exit - 
C if FAC = 0.0d0 there is no point in doing this  ...
C
C     .
      IF ( ABS( FAC ).LT.LOW ) RETURN
C     .
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C     .
C     . Check validity of NRR. Because it is an array
C     . dimension, NR has to be passed in the parameter
C     . list and so we must ensure we are consistent
C     .
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine AMLICA.'
         PRINT *,' INARR( 2 ) = ', NRR
         PRINT *,' NR = ', NR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Check that we really are dealing with
C     . a non-uniform grid.
C     .
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
        PRINT *,' Subroutine AMLICA. IFORMF = ',IFORMF
        PRINT *,' This is not the uniform grid code.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check size of IHD
C     .
      IF ( IHD.GT.NDRVS ) THEN
        PRINT *,' Subroutine AMLICA. '
        PRINT *,' IHD   = ', IHD
        PRINT *,' NDRVS = ', NDRVS
        STOP
      ENDIF
C     .
C     . Check N2 for case of IFORMF = 3 and 4
C     .
      IF ( IFORMF.EQ.3 .OR. IFORMF.EQ.4 ) THEN
        IF ( N2.NE.NH*NR ) THEN
          PRINT *,' Subroutine AMLICA. N2 = ', N2
          PRINT *,' NH = ', NH,' NR = ', NR
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDIF
C     .
C     . Check N1 for IMF.EQ.1 and IMF.EQ.3
C     .
      IF ( ( IMF.EQ.1 .AND. N1.NE.(KLE+KL+KU+1) ) .OR.
     1     ( IMF.EQ.3 .AND. N1.NE.N2           )      ) THEN
         PRINT *,' Subroutine AMLICA. IMF = ',IMF
         PRINT *,' N1 = ', N1,' N2 = ',N2
         PRINT *,' KL = ', KL,' KU = ',KU
         PRINT *,' KLE = ', KLE
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Now we will check that we have sufficient
C     . band width to be able to calculate the
C     . IHD^{th} derivative of a function.
C     . IHD must be less than or equal to the
C     . minimum value of (NLN + NRN)
C     . 
C     . To do this, we calculate NLCS and NRCS
C     . which are respectively the minimum number
C     . of points you will ever find on the left
C     . and the right of your node, IRAD.
C     .
      NLCS = ILNR - 1
      NRCS = NR   - IRNR
C
C Just check that we have sufficiently many grid nodes
C to be able to enforce all equations.
C
      IF ( NLCS.LT.0 ) THEN
        PRINT *,' Subroutine AMLICA.'
        PRINT *,' ILNR = ', ILNR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NRCS.LT.0 ) THEN
        PRINT *,' Subroutine AMLICA.'
        PRINT *,' IRNR = ', IRNR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      I = MIN( NLCS, NRCS) + NBN
      IF ( IHD.GT.I ) THEN
        PRINT *,' Subroutine AMLICA.'
        PRINT *,' You want a derivative of order ', IHD
        PRINT *,' However, there is a node at which you '
        PRINT *,' only have ',I,' points for a derivative.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Ok - we now loop around the radial nodes for the
C     . rows ... 
C     .
      DO IRAD = ILNR, IRNR
C       .
C       . Calculate the corresponding index and radius
C       .
        INDR = INDFUN( IRAD, IHR, INARR )
        RAD = XARR( IRAD )
C       .
C       . Calculate NLN and NRN
C       .
        NLN = MIN( NBN, IRAD - 1 )
        NRN = MIN( NBN, NR - IRAD )
C       .
        CALL SUB1( WORK, RAD, IPARS, DPARS, IHD )
C       .
C       . So now loop around the nodes from (IRAD - NLN)
C       . to (IRAD + NRN)
C       .
        DO INODE = IRAD - NLN, IRAD + NRN
C          .
C          . This is the row of the finite
C          . difference coefficients matrix
C          . which will contain the correct values.
C          .
           ICOCOL = INODE - IRAD + NBN + 1
C          .
C          . find the corresponding matrix column
C          .
           INDC = INDFUN( INODE, IHC, INARR )
C          .
C          . Find the actual location in the AMAT array
C          .
           CALL MATIND (INDR,INDC,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
C          .
           DO NDER = 0, IHD
             ICOROW = NDER + 1
             AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + 
     1         FAC*WORK( ICOROW )*SVFDC( ICOCOL, IRAD, ICOROW, IS )
           ENDDO
C          .
        ENDDO
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Adapted Matrix Surplus Diagonal Element Addition ********
C            -       -      -       -        -       -        ********
C Steve Gibbons Wed Oct 27 09:12:14 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C By using finite difference schemes where homogeneous boundary      C
C conditions are implicit in the coefficients used at points close   C
C to the boundary, either the end point or the two end points at     C
C each boundary become obsolete and must be arbitrarily filled       C
C to prevent singularity of the matrix.                              C
C                                                                    C
C AMSDEA loops around each harmonic, checks to see which scheme is   C
C used for this harmonic, looks to see what condition is implied     C
C at the boundary and then adds 0, 1 or 2 diagonal values of         C
C DIAGEL. If we are merely solving a linear system, DIAGEL is        C
C arbitrary. If we solving an eigensystem, DIAGEL will become an     C
C eigenvalue and so should be as highly negative as is necessary     C
C to prevent confusion with the 'real' eigenvalues.                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C                                                                    C
C     IMF       : Matrix format flag.                                C
C                                                                    C
C          imf = 1; Matrix is in LAPACK banded format                C
C                   ie element a_{i,j} is stored in                  C
C                   A( kle + ku + 1 + i - j , j )                    C
C                                                                    C
C    ALL other values of IMF are disqualified for this routine.      C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1, 3. INDFUN = ( IR - 1 )*NH + IH       C
C                   IFORMF = 2, 4. INDFUN = ( IH - 1 )*NR + IR       C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHPI( ih ) = is, which is the 4th index of        C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     IBCARR    : Dimension ( NDCS ). This must be the array MHIBC   C
C                 submitted to SVFDCF in the case CHBND(1:1) = 'I'   C
C                 and the array MHOBC submitted to SVFDCF in the     C
C                 case CHBND(1:1) = 'O'                              C
C                                                                    C
C     IBCARR( IS ) contains the integer IBC which corresponds        C
C    to the following boundary conditions :-                         C
C                                                                    C
C   IBC       Boundary condition                   No. of nodes      C
C   ---       ------------------                   ------------      C
C                                                                    C
C    1           None                                    0           C
C    2           f = 0                                   1           C
C    3           df/dr = 0                               1           C
C    4           f = 0   and     df/dr = 0               2           C
C    5           f = 0   and   d^2f/dr^2 = 0             2           C
C    6           f/r - df/dr = 0                         1           C
C    7           Insulating magnetic field               1           C
C                                                                    C
C                                                                    C
C                                                                    C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHBND     : Boundary flag (*). Either 'Inner' or 'Outer'       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 )                             C
C     DIAGEL    : D.p. const. to be added to diagonals.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AMSDEA( A, N1, N2, KL, KU, KLE, IMF, INARR, 
     1                   MHP, IBCARR, CHBND, DIAGEL, NDCS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, IMF, INARR( * ), MHP( * ), NDCS,
     1        IBCARR( NDCS )
      DOUBLE PRECISION A( N1, N2 ), DIAGEL
      CHARACTER *(*)   CHBND
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NR, NH, IRAD, INDFUN, ISN, IEN, IH, IS, IOF,
     1        NRC, IRC, IBC
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR = INARR( 2 )
      NH = INARR( 3 )
C     .
C     . Check on IMF (1 is only acceptable value)
C     .
      IF ( IMF.NE.1 ) THEN
        PRINT *,' Subroutine AMSDEA'
        PRINT *,' IMF = ', IMF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check on size of N2
C     .
      IF ( N2.NE.NH*NR ) THEN
        PRINT *,' Subroutine AMSDEA'
        PRINT *,' N2 = ', N2
        PRINT *,' NR = ', NR
        PRINT *,' NH = ', NH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Set inner outer flag
C     .
      IF ( CHBND(1:1).EQ.'I' .OR. CHBND(1:1).EQ.'i' ) THEN
        IOF = 1
        GOTO 50
      ENDIF
C     .
      IF ( CHBND(1:1).EQ.'O' .OR. CHBND(1:1).EQ.'o' ) THEN
        IOF = 2
        GOTO 50
      ENDIF
C     .
      PRINT *,' Subroutine AMSDEA'
      PRINT *,' CHBND = ', CHBND
      PRINT *,' Program aborted.'
      STOP
C     .
 50   CONTINUE
C     .
C     . Begin loop around harmonics
C     .
      DO IH = 1, NH
        IS  = MHP( IH )
        IF ( IS.LT.1 .OR. IS.GT.NDCS ) THEN
          PRINT *,' Subroutine AMSDEA'
          PRINT *,' IS   = ', IS
          PRINT *,' NDCS = ', NDCS
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        IBC = IBCARR( IS )
        IF ( IBC.EQ.1 ) GOTO 51
C ............................ inner boundary case
C
        IF ( IOF.EQ.1 ) THEN
C         .
          IF ( IBC.EQ.2 .OR. IBC.EQ.3 .OR. IBC.EQ.6
     1           .OR. IBC.EQ.7 ) THEN
             ISN = 1
             IEN = 1
             GOTO 61
          ENDIF
C         .
          IF ( IBC.EQ.4 .OR. IBC.EQ.5 ) THEN
             ISN = 1
             IEN = 2
             GOTO 61
          ENDIF
C         .
          PRINT *,' Subroutine AMSDEA'
          PRINT *,' IBCARR(',IS,') = ', IBC
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C ............................ outer boundary case
C
        IF ( IOF.EQ.2 ) THEN
C         .
          IF ( IBC.EQ.2 .OR. IBC.EQ.3 .OR. IBC.EQ.6
     1           .OR. IBC.EQ.7 ) THEN
             ISN = NR
             IEN = NR
             GOTO 61
          ENDIF
C         .
          IF ( IBC.EQ.4 .OR. IBC.EQ.5 ) THEN
             ISN = NR - 1
             IEN = NR
             GOTO 61
          ENDIF
C         .
          PRINT *,' Subroutine AMSDEA'
          PRINT *,' IBCARR(',IS,') = ', IBC
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C ............................ 
 61     CONTINUE
        DO IRAD = ISN, IEN
C           .
C           . Zero the row of the matrix
C           .
            IRC = 1
            NRC = INDFUN( IRAD, IH, INARR )
            CALL BMRCOP( KL, KU, KLE, N2, IRC, NRC, A,
     1                   ZERO, ZERO )
C           .
C           . Enter value of the diagonal element
C           .
            A( KLE + KU + 1, NRC ) = DIAGEL
C           .
        ENDDO
C ............................ 
 51   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Adapted Matrix DL Terms *********************************
C            -       -      -- -     *********************************
C Steve Gibbons Tue Oct 26 15:24:17 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Supplies coefficients to the routine AMLICA for terms with the     C
C function D_l to the power 0, 1 or 2. Must be declared EXTERNAL     C
C in the (sub)program which calls AMLICA and must replave SUB1 in    C
C the calling sequence.                                              C
C                                                                    C
C Set IPARS( 1 ) to IFLAG                                            C
C                                                                    C
C Now if IFLAG = 1, D_l^0 is added to the matrix - i.e. the identity C
C                   IHD can be zero.                                 C
C                                                                    C
C        IFLAG = 2, D_l is added to the matrix.                      C
C                   This returns                                     C
C                                                                    C
C                     CVEC( 1 ) = (-1.0d0)*DBLE( L*L+L )/(RAD*RAD)   C
C                     CVEC( 2 ) =   2.0d0/RAD                        C
C                     CVEC( 3 ) =   1.0d0                            C
C                                                                    C
C                   IHD must be atleast 2.                           C
C                                                                    C
C        IFLAG = 3, D_l^2 is added to the matrix                     C
C                                                                    C
C                     CVEC( 1 ) = (L+2)(L+1)L(L-1)/RAD**4            C
C                     CVEC( 2 ) =  0.0d0                             C
C                     CVEC( 3 ) = -2L(L+1)/RAD**2                    C
C                     CVEC( 4 ) = 4/RAD                              C
C                     CVEC( 5 ) =  1.0d0                             C
C                                                                    C
C  The integer L is passed in as IPARS( 2 ).                         C
C                                                                    C
C  All other parameters are looked after by AMLICA and most are      C
C  dummy parameters.                                                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AMDLT( CVEC, RAD, IPARS, DPARS, IHD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IPARS( * ), IHD
      DOUBLE PRECISION CVEC( * ), RAD, DPARS( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IHDMIN( 3 ), IFLAG, L, ND
      DOUBLE PRECISION TOL
      PARAMETER ( TOL = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check on value of RAD
C     .
      IF ( ABS( RAD ).LT.TOL ) THEN
         PRINT *,' Subroutine AMDLT.'
         PRINT *,' RAD = ', RAD
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C
C Just apply the minimum values of IHD to the IHDMIN array.
C
C     .
      IHDMIN( 1 ) = 0
      IHDMIN( 2 ) = 2
      IHDMIN( 3 ) = 4
C     .
      IFLAG = IPARS( 1 )
      L     = IPARS( 2 )
C     .
C     . Check for valid value of IFLAG
C     .
      IF ( IFLAG.LT.1 .OR. IFLAG.GT.3 ) THEN
         PRINT *,' Subroutine AMDLT.'
         PRINT *,' IFLAG = ', IFLAG
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . ok - so our chosen boundary condition is valid
C     . now check that IHD is large enough
C     .
      IF (     IHD.NE.IHDMIN( IFLAG )    ) THEN
         PRINT *,' Subroutine AMDLT.'
         PRINT *,' IHD = ',IHD,' and must be atleast'
         PRINT *,  IHDMIN( IFLAG ),' for option ',IFLAG
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Zero all coefficients up to IHD
C     .
      DO ND = 1, IHD + 1
        CVEC( ND ) = 0.0d0
      ENDDO
C     .
C     . Option 1:
C     . D_l^0
C     .
      IF ( IFLAG.EQ.1 ) THEN
         CVEC( 1 ) = 1.0d0
         RETURN
      ENDIF
C     .
C     . Option 2:
C     . D_l
C     .
      IF ( IFLAG.EQ.2 ) THEN
         CVEC( 1 ) = (-1.0d0)*DBLE( L*L+L )/(RAD*RAD)
         CVEC( 2 ) = 2.0d0/RAD
         CVEC( 3 ) = 1.0d0
         RETURN
      ENDIF
C     .
C     . Option 3:
C     . D_l^2
C     .
      IF ( IFLAG.EQ.3 ) THEN
         CVEC( 1 ) = DBLE( (L+2)*(L+1)*L*(L-1) )/RAD**4
         CVEC( 3 ) = (-2.0d0)*DBLE( L*L + L )/RAD**2
         CVEC( 4 ) = 4.0d0/RAD
         CVEC( 5 ) = 1.0d0
         RETURN
      ENDIF
C     .
      END
C*********************************************************************

C*********************************************************************
C subroutine Vector Operations Banded Matrix Auxiliary Routine *******
C            -      -          -      -      -         -       *******
C Steve Gibbons Mon Nov 27 15:20:35 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C We wish to form banded matrices with which to perform many         C
C different operations. This subroutine may be called as an          C
C EXTERNAL routine from subroutines such as AMLICA and NMLICA        C
C in order to form these matrices.                                   C
C                                                                    C
C DPARS is not referred to.                                          C
C                                                                    C
C RAD is the radius.                                                 C
C                                                                    C
C IPARS( 1 ) contains the integer IOPT which specifies which         C
C matrix we are building. This may take any of the values listed     C
C below.                                                             C
C                                                                    C
C IPARS( 2 ) contains the spherical harmonic degree, l.              C
C                                                                    C
C IHD is dependent upon IOPT (see below)                             C
C                                                                    C
C Possible values for IOPT are                                       C
C ----------------------------                                       C
C                                                                    C
C  1: multiply p(r) to get Q(r)                                      C
C                                                                    C
C        now Q( r ) = L(L+1)/RAD p( r )                              C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = L(L+1)/RAD                    C
C                 *******  **********************                    C
C                                                                    C
C  2: multiply p(r) to get S(r)                                      C
C                                                                    C
C        now S( r ) = DSQRT( L(L+1) )( p/RAD + dp/dr )               C
C                                                                    C
C        and so   IHD = 1, CVEC( 1 ) = DSQRT( L(L+1) )/RAD           C
C                 *******  CVEC( 2 ) = DSQRT( L(L+1) )               C
C                          ***************************               C
C                                                                    C
C  3: multiply tau(r) to get T(r)                                    C
C                                                                    C
C        now T( r ) = -DSQRT( L(L+1) ) tau( r )                      C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = -DSQRT( L(L+1) )              C
C                 *******  ****************************              C
C                                                                    C
C In the following q, s and t refer to the vectors                   C
C                                                                    C
C  4: curl of scaloidal vector Q(r) q                                C
C                                                                    C
C       curl[ Q(r) q ] = -DSQRT( L(L+1) ) Q( r )/RAD  t              C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = -DSQRT( L(L+1) )/RAD          C
C                 *******  ********************************          C
C                                                                    C
C  5: curl of spheroidal vector S(r) s                               C
C                                                                    C
C       curl[ S(r) s ] = t[ dS/dr + S(r)/RAD ]                       C
C                                                                    C
C        and so   IHD = 1, CVEC( 1 ) = 1/RAD                         C
C                 *******  CVEC( 2 ) = 1.0d0                         C
C                          *****************                         C
C                                                                    C
C  6: curl of toroidal vector T(r) t: (scaloidal component)          C
C                                                                    C
C       curl[ T(r) t ] = - q[ T(r) DSQRT( L(L+1) )/RAD ]             C
C                        - s[ dT/dr + T(r)/RAD ]                     C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = -DSQRT( L(L+1) )/RAD          C
C                 *******  ********************************          C
C                                                                    C
C  7: curl of toroidal vector T(r) t: (spheroidal component)         C
C                                                                    C
C       curl[ T(r) t ] = - q[ T(r) DSQRT( L(L+1) )/RAD ]             C
C                        - s[ dT/dr + T(r)/RAD ]                     C
C                                                                    C
C        and so   IHD = 1, CVEC( 1 ) = - 1.0d0 / RAD                 C
C                 *******  CVEC( 2 ) = - 1.0d0                       C
C                          *************************                 C
C                                                                    C
C  8: Mulitply Q(r) to get p(r)                                      C
C                                                                    C
C       p( r ) = RAD/( L*L + L ) Q( r )                              C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = RAD/( L*L + L )               C
C                 *******  ***************************               C
C                                                                    C
C  9: Mulitply T(r) to get tau(r)                                    C
C                                                                    C
C       tau( r ) = - T( r ) / SQRT( L*L + L )                        C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = - 1 / SQRT( L*L + L )         C
C                 *******  *********************************         C
C                                                                    C
C 10: Calculate a pure first derivative.                             C
C                                                                    C
C        and so   IHD = 1, CVEC( 1 ) = 0                             C
C                 *******  CVEC( 2 ) = 1.0d0                         C
C                          *****************                         C
C                                                                    C
C 11: Curl of scaloidal vector Q(r) q to give tau( r ) radial        C
C      function                                                      C
C                                                                    C
C        curl[ Q(r) q ] = T_c( r ) t                                 C
C                                                                    C
C     where T_c( r ) = -DSQRT( L(L+1) ) Q( r )/RAD                   C
C                                                                    C
C     and so tau( r ) = Q( r )/RAD                                   C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = 1.0d0 / RAD                   C
C                 *******  ***********************                   C
C                                                                    C
C 12: Curl of spheroidal vector S(r) s to give tau( r ) radial       C
C      function                                                      C
C                                                                    C
C        curl[ S(r) s ] = T_c( r ) t                                 C
C                                                                    C
C     where T_c( r ) = dS/dr + S(r)/RAD                              C
C                                                                    C
C       and so tau( r ) =  -( dS/dr + S(r)/RAD )/DSQRT( L(L+1) )     C
C                                                                    C
C        and so   IHD = 1, CVEC( 1 ) = -1.0d0/(RAD*DSQRT( L(L+1) ))  C
C                 *******  CVEC( 2 ) = -1.0d0/DSQRT( L(L+1) )        C
C                          ****************************************  C
C                                                                    C
C 13: Curl of toroidal vector T( r ) t to give p( r )                C
C                                                                    C
C        curl[ T(r) t ] = Q_c( r ) q   +   S_c( r ) s                C
C                                                                    C
C         where  Q_c( r ) = -T(r) DSQRT( L(L+1) )/RAD   and          C
C                S_c( r ) = -[ dT/dr + T(r)/RAD ]                    C
C                                                                    C
C   hence p( r ) = -T(r)/DSQRT( L(L+1) )                             C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = -1.0d0/DSQRT( L(L+1) )        C
C                 *******  **********************************        C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VOBMAR( CVEC, RAD, IPARS, DPARS, IHD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IPARS( * ), IHD
      DOUBLE PRECISION CVEC( * ), RAD, DPARS( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IHDMIN( 13 ), IOPT, L, ND
      DOUBLE PRECISION TOL, DLFAC, DLL1, DSLL1
      PARAMETER ( TOL = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check on value of RAD
C     .
      IF ( DABS( RAD ).LT.TOL ) THEN
        PRINT *,' Subroutine VOBMAR.'
        PRINT *,' RAD = ', RAD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IOPT = IPARS( 1 )
      L    = IPARS( 2 )
C     .
      IF ( IOPT.LT.1 .OR. IOPT.GT.13 ) THEN
        PRINT *,' Subroutine VOBMAR.'
        PRINT *,' IOPT = ', IOPT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      DLFAC = DBLE( L )
      DLL1  = DLFAC*DLFAC + DLFAC
      DSLL1 = DSQRT( DLL1 )
C     .
      IHDMIN(  1 ) = 0
      IHDMIN(  2 ) = 1
      IHDMIN(  3 ) = 0
      IHDMIN(  4 ) = 0
      IHDMIN(  5 ) = 1
      IHDMIN(  6 ) = 0
      IHDMIN(  7 ) = 1
      IHDMIN(  8 ) = 0
      IHDMIN(  9 ) = 0
      IHDMIN( 10 ) = 1
      IHDMIN( 11 ) = 0
      IHDMIN( 12 ) = 1
      IHDMIN( 13 ) = 0
C     .
      IF ( IHD.NE.IHDMIN( IOPT ) ) THEN
        PRINT *,' Subroutine VOBMAR.'
        PRINT *,' IHD = ', IOPT
        PRINT *,' Should be ',IHDMIN( IOPT ),' for option.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      DO ND = 1, IHD + 1
        CVEC( ND ) = 0.0d0
      ENDDO
C     .
      IF ( IOPT.EQ.1 ) THEN
        CVEC( 1 ) = DLL1/RAD
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.2 ) THEN
        CVEC( 1 ) = DSLL1/RAD
        CVEC( 2 ) = DSLL1
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.3 ) THEN
        CVEC( 1 ) = (-1.0d0)*DSLL1
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.4 ) THEN
        CVEC( 1 ) = (-1.0d0)*DSLL1/RAD
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.5 ) THEN
        CVEC( 1 ) = 1.0d0/RAD
        CVEC( 2 ) = 1.0d0
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.6 ) THEN
        CVEC( 1 ) = (-1.0d0)*DSLL1/RAD
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.7 ) THEN
        CVEC( 1 ) = (-1.0d0)/RAD
        CVEC( 2 ) = (-1.0d0)
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.8 ) THEN
        CVEC( 1 ) = RAD/DLL1
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.9 ) THEN
        CVEC( 1 ) = (-1.0d0)/DSLL1
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.10 ) THEN
        CVEC( 1 ) = 0.0d0
        CVEC( 2 ) = 1.0d0
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.11 ) THEN
        CVEC( 1 ) = 1.0d0/RAD
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.12 ) THEN
        CVEC( 1 ) = -1.0d0/(DSLL1*RAD)
        CVEC( 2 ) = -1.0d0/DSLL1
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.13 ) THEN
        CVEC( 1 ) = -1.0d0/DSLL1
        RETURN
      ENDIF
C     .
      END
C*********************************************************************
C*********************************************************************
C subroutine Non-adapted Matrix Linear Interaction Contribution Add **
C            -           -      -      -           -            -   **
C Steve Gibbons Mon Nov 13 13:20:34 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C  Adds to the double precision matrix, A, the terms which a single  C
C harmonic IHC contributes to a single harmonic IHR in a single term C
C of a coupled o.d.e. when the grid in radius is not necessarily     C
C uniform and the j^[th} abscissa has a radial value of XARR( j )    C
C                                                                    C
C Let the radial function of the IHR harmonic be denoted f_{IHR)(r)  C
C Let the radial function of the IHC harmonic be denoted f_{IHC)(r)  C
C                                                                    C
C If the contribution for a given radial node, $r_j$, is of the form C
C                                                                    C
C f_{IHR)(r_j) = \sum_{nd=0,IHD} f_{IHC}^{nd}(r_j) c_{IHC}^{nd}(r_j) C
C                                                                    C
C  (where f^{nd} denotes the (nd)^{th} derivative of f with respect  C
C  to r, evaluated at r_j)                                           C
C                                                                    C
C then the coefficient c_{IHC}^{nd}(r_j) must be supplied by the     C
C subroutine SUB1 (declared EXTERNAL in the calling (sub)program )   C
C which returns c_{IHC}^{nd}(r_j) in the array element CVEC(nd + 1). C
C                                                                    C
C The subroutine SUB1 *MUST* have the calling sequence               C
C                                                                    C
C SUB1( CVEC, RAD, IPARS, DPARS, IHD )                               C
C                                                                    C
C where  RAD is the double precision value of the radius,            C
C        IPARS is an integer array to provide SUB1 with parameters.  C
C          IPARS is not referred to by NMLICA other than to pass     C
C          this information.                                         C
C        DPARS is like IPARS but contains double precision elements. C
C        IHD is the highest derivative which will be refered to by   C
C        NMLICA. So all CVEC( i ) must be zero with 1.le.i.le.ihd+1  C
C        if it is not assigned another value.                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C                                                                    C
C     IMF       : Matrix format flag.                                C
C                                                                    C
C          imf = 1; Matrix is in LAPACK banded format                C
C                   ie element a_{i,j} is stored in                  C
C                   A( kle + ku + 1 + i - j , j )                    C
C                                                                    C
C          imf = 2; Matrix is banded but with element a_{i,j}        C
C                   stored in A( kl + 1 + j - i , i ).               C
C                                                                    C
C          imf = 3; Matrix is square - ie a_{i,j} is stored          C
C                   in A( i, j ).                                    C
C                                                                    C
C     IHC       : Number of the input (column) harmonic              C
C     IHR       : Number of the output (row) harmonic                C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1. INDFUN = ( IR - 1 )*NH + IH          C
C                   IFORMF = 2. INDFUN = ( IH - 1 )*NR + IR          C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     IHD       : Highest derivative involved.                       C
C     NBN       : Number of diagonal elements in radius.             C
C     ILNR      : Left-most node at which to start adding rows to    C
C                   matrix.                                          C
C     IRNR      : Right-most node at which to start adding rows to   C
C                   matrix.                                          C
C     ILNC      : Left-most node which may be used to form deriv.s   C
C     IRNC      : Right-most node which may be used to form deriv.s  C
C                                                                    C
C  ILNC must never be greater than ILNR, and similarly,              C
C  IRNC must never be less than IRNR.                                C
C                                                                    C
C  ILNC and IRNC must be the same as in the call to the routine      C
C  fdcmbd when calculating the array FDCM. In matrix problems,       C
C  ILNC and IRNC should almost always be set to 1 and NR resp.       C
C  Other options are generally written for the case where a          C
C  function of r may not be defined at one extremity.                C
C                                                                    C
C     NFDCM     : Leading dimension of FDCM.                         C
C     NDRVS     : Highest derivative stored in FDCM.                 C
C                                                                    C
C     IPARS     : Information array for SUB1. Dimension ( * )        C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C  Subroutines                                                       C
C  -----------                                                       C
C                                                                    C
C     SUB1      : Determines what multiplies each derivative in      C
C                  the matrix. Must have calling sequence ...        C
C                                                                    C
C     CALL SUB1 ( CVEC, RAD, IPARS, DPARS, IHD )                     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     AMAT      : Matrix. Dimensions ( N1, N2 )                      C
C              Will generally be banded due to the nature of the     C
C              numerical scheme. KL, KU and KLE parameterise this.   C
C                                                                    C
C     FAC       : Multiplication factor.                             C
C                                                                    C
C     XARR      : Array of dimension ( NR )                          C
C                 XARR( j ) = element x_j                            C
C                                                                    C
C     FDCM      : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVS ).                   C
C                   Array is generated by the routine fdcmbd         C
C                 See documentation for FDCMBD for details.          C
C                                                                    C
C     DPARS     : Information array for SUB1. Dimension ( * )        C
C                                                                    C
C     WORK      : Working array of dimension ( NDRVS + 1 )           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NMLICA( N1, N2, KL, KU, KLE, IMF, IHC, IHR, INARR,
     1                   IHD, NBN, ILNR, IRNR, ILNC, IRNC, NFDCM,
     2                   NDRVS, SUB1, AMAT, FAC, XARR, FDCM, DPARS,
     3                   IPARS, NR, WORK )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, IMF, IHC, IHR, INARR( * ), IHD,
     1        NBN, ILNR, IRNR, ILNC, IRNC, NFDCM, NDRVS, IPARS( * ),
     2        NR
      EXTERNAL SUB1
      DOUBLE PRECISION AMAT( N1, N2 ), FAC, XARR( NR ),
     1                 FDCM( NFDCM, NR, * ), DPARS( * ),
     2                 WORK( NDRVS + 1 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IRAD, INODE, NLN, NRN, INDFUN, INDR, INDC,
     1        IROW, ICOL, NH, IFORMF, NLCS, NRCS, I,
     2        NDER, ICOROW, ICOCOL, NRR
      DOUBLE PRECISION RAD, ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First we check for an easy exit - 
C if FAC = 0.0d0 there is no point in doing this  ...
C
C     .
      IF ( FAC.EQ.ZERO ) RETURN
C     .
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C     .
C     . Check validity of NRR. Because it is an array
C     . dimension, NR has to be passed in the parameter
C     . list and so we must ensure we are consistent
C     .
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine NMLICA.'
         PRINT *,' INARR( 2 ) = ', NRR
         PRINT *,' NR = ', NR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Check size of IHD
C     .
      IF ( IHD.GT.NDRVS ) THEN
        PRINT *,' Subroutine NMLICA. '
        PRINT *,' IHD   = ', IHD
        PRINT *,' NDRVS = ', NDRVS
        STOP
      ENDIF
C     .
C     . Check N2 for case of IFORMF = 1,2,3 and 4
C     .
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 .OR. IFORMF.EQ.3 .OR.
     1                     IFORMF.EQ.4            ) THEN
        IF ( N2.NE.NH*NR ) THEN
          PRINT *,' Subroutine NMLICA. N2 = ', N2
          PRINT *,' NH = ', NH,' NR = ', NR
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDIF
C     .
C     . Check N1 for IMF.EQ.1 and IMF.EQ.3
C     .
      IF ( ( IMF.EQ.1 .AND. N1.NE.(KLE+KL+KU+1) ) .OR.
     1     ( IMF.EQ.3 .AND. N1.NE.N2           )      ) THEN
         PRINT *,' Subroutine NMLICA. IMF = ',IMF
         PRINT *,' N1 = ', N1,' N2 = ',N2
         PRINT *,' KL = ', KL,' KU = ',KU
         PRINT *,' KLE = ', KLE
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Now we will check that we have sufficient
C     . band width to be able to calculate the
C     . IHD^{th} derivative of a function.
C     . IHD must be less than or equal to the
C     . minimum value of (NLN + NRN)
C     . 
C     . To do this, we calculate NLCS and NRCS
C     . which are respectively the minimum number
C     . of points you will ever find on the left
C     . and the right of your node, IRAD.
C     .
      NLCS = ILNR - ILNC
      NRCS = IRNC - IRNR
C
C Just check that we have sufficiently many grid nodes
C to be able to enforce all equations.
C
      IF ( NLCS.LT.0 ) THEN
        PRINT *,' Subroutine NMLICA.'
        PRINT *,' ILNR = ', ILNR
        PRINT *,' ILNC = ', ILNC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NRCS.LT.0 ) THEN
        PRINT *,' Subroutine NMLICA.'
        PRINT *,' IRNR = ', IRNR
        PRINT *,' IRNC = ', IRNC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      I = MIN( NLCS, NRCS) + NBN
      IF ( IHD.GT.I ) THEN
        PRINT *,' Subroutine NMLICA.'
        PRINT *,' You want a derivative of order ', IHD
        PRINT *,' However, there is a node at which you '
        PRINT *,' only have ',I,' points for a derivative.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Ok - we now loop around the radial nodes for the
C     . rows ... first set our logical flag to .FALSE.
C     . i.e. calculate differencing coefficients each
C     . time unless told otherwise
C     .
      DO IRAD = ILNR, IRNR
C       .
C       . Calculate the corresponding index and radius
C       .
        INDR = INDFUN( IRAD, IHR, INARR )
        RAD  = XARR( IRAD )
C       .
C       . Calculate NLN and NRN
C       .
        NLN = MIN( NBN, IRAD - ILNC )
        NRN = MIN( NBN, IRNC - IRAD )
C       .
        CALL SUB1( WORK, RAD, IPARS, DPARS, IHD )
C       .
C       . So now loop around the nodes from (IRAD - NLN)
C       . to (IRAD + NRN)
C       .
        DO INODE = IRAD - NLN, IRAD + NRN
C          .
C          . This is the row of the finite
C          . difference coefficients matrix
C          . which will contain the correct values.
C          .
           ICOCOL = INODE - IRAD + NBN + 1
C          .
C          . find the corresponding matrix column
C          .
           INDC = INDFUN( INODE, IHC, INARR )
C          .
C          . Find the actual location in the AMAT array
C          .
           CALL MATIND (INDR,INDC,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
C          .
C          . Add the zero derivative part ...
C          .
           IF ( INODE.EQ.IRAD ) THEN
             AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) +
     1           FAC*WORK( 1 )
           ENDIF
C          .
           DO NDER = 1, IHD
             ICOROW = NDER + 1
             AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + 
     1           FAC*WORK( ICOROW )*FDCM( ICOCOL, IRAD, NDER )
           ENDDO
C          .
        ENDDO
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Double precision VECtor Zero ****************************
C Steve Gibbons Wed Feb 14 09:01:11 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Sets to zero a double precision vector, VEC, of length N.          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N		: Length of the vector.                              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VEC	: Vector - dimension ( N )                           C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DVECZ( VEC, N )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          N
      DOUBLE PRECISION VEC( N )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          I
      DOUBLE PRECISION DZERO
      PARAMETER      ( DZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO I = 1, N
        VEC( I ) = DZERO
      ENDDO
      RETURN
      END
C*********************************************************************


C*********************************************************************
C subroutine Fast Fourier Transform ReaL Version *********************
C            -    -       -         -  - -       *********************
C Steve Gibbons 22.4.97 (Based on FOUR1 from Numerical Recipes       C
C      I've made slight alterations for double precision, error      C
C      checking at the start and rescaling the coefficients after    C
C      the forward transform.)                                       C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NN	: Number of data to be transformed.		     C
C                  NN MUST be a power of two. This is checked for by C
C                  the subroutine POWTWO.                            C
C                                                                    C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       C
C                                                                    C
C    The real function f( PHI ) must be periodic, satisfying         C
C    f( PHI ) = f( PHI + 2*PI )                                      C
C                                                                    C
C    f( PHI ) has the formulation                                    C
C                                                                    C
C    f( PHI ) = sum_{m=0}^M [ c_m cos(m PHI) + s_m sin( m PHI) ]     C
C                                                                    C
C    FFTRLV will transform between the coefficients {c_m, s_m}       C
C    and discretised values of f at NN equally spaced points PHI_i   C
C    with                                                            C
C                                                                    C
C    PHI_i = (i-1)*DELTA_phi          and                            C
C                                                                    C
C    DELTA_phi = 2*PI/dble( NN )                                     C
C                                                                    C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       C
C                                                                    C
C     ISIGN	:                                                    C
C                                                                    C
C    If ISIGN = 1, on entry DATA( 2*I - 1 ) must contain the value   C
C    of f at the point PHI = PHI_i as defined above.                 C
C                                                                    C
C    On exit, DATA( 2*m + 1 ) will contain the coeff. c_m and        C
C             DATA( 2*m + 2 ) will contain the coeff. s_m and        C
C                                                                    C
C    for m = 0, M.                                                   C
C      For accuracy, M **MUST** be less than NN/2                    C
C    (This is the Nyquist criterion - see Numerical Recipes ).       C
C                                                                    C
C    If ISIGN = -1, DATA must be zero except for the values          C
C    c_m contained in DATA( 2*m + 1 )  and                           C
C    s_m contained in DATA( 2*m + 2 ) on input.                      C
C                                                                    C
C    On output, f( PHI_i ) will be contained in DATA( 2*I - 1 ).     C
C                                                                    C
C    Note that this routine SCALES the coefficients after            C
C    transforming, which the Numerical Recipes version doesn't.      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     DATA	: Array of dimension (2*NN). See above ...           C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FFTRLV ( DATA, NN, ISIGN )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NN, ISIGN
      DOUBLE PRECISION DATA( 2 * NN )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TEMPR,TEMPI,THETA,WPR,WPI,WR,WI,WTEMP,FAC
      INTEGER I,J,N,M,MMAX,ISTEP
c     LOGICAL POWT
      DOUBLE PRECISION PI, HTHETA
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C 
C  SJG Mon Feb  5 18:18:33 WET 2001
C  I have commented out the power of two checking
C  for speed: reinstate if necessary:
C 
C First check that the number NN supplied is infact a power of 2.
c     CALL POWTWO ( NN, POWT)
c     IF ( .NOT. POWT) THEN
c        PRINT *,' Subroutine FFTRLV. NN is not a power of 2.'
c        PRINT *,'Program aborted.'
c        STOP
c     ENDIF
C Now check that ISIGN = 1 or -1.
      IF ( ISIGN*ISIGN .NE. 1) THEN
         PRINT *,' Subroutine FFTRLV. ISIGN must be 1 or -1.'
         PRINT *,'Program aborted.'
         STOP
      ENDIF
C____________________________________________________________________C
C Here we will zero the imaginary part of the function
C just incase it contains non-zero data which may harm the
C outcome.
C
      IF ( ISIGN.EQ.1 ) THEN
        DO I = 1, NN
          DATA( 2*I ) = 0.0d0
        ENDDO
      ENDIF
C
C____________________________________________________________________C
C
      N = 2 * NN
      J = 1
      DO I = 1, N, 2
         IF ( J.GT.I ) THEN
            TEMPR = DATA( J )
            TEMPI = DATA( J+1 )
            DATA( J ) = DATA( I )
            DATA( J+1 ) = DATA( I+1 )
            DATA( I ) = TEMPR
            DATA( I+1 ) = TEMPI
         ENDIF
C ... now calculate inverse bit map for next value of I.
         M = N/2
 500     IF ( M.GE.2 .AND. J.GT.M ) THEN
            J = J-M
            M = M/2
            GOTO 500
         ENDIF
         J = J + M
      ENDDO
C ............. now do the real part of transform.
      MMAX = 2
 502  IF ( N.GT.MMAX ) THEN
         ISTEP = 2*MMAX
         HTHETA = PI/DBLE(ISIGN*MMAX)
c        THETA = 2.0d0*PI/DBLE(ISIGN*MMAX)
C set THETA temporarily
         THETA = DSIN(HTHETA)
         WPR = -2.0d0*THETA*THETA
c        WPR = -2.0d0*DSIN(HTHETA)**2
         THETA = 2.0d0*HTHETA
         WPI = DSIN(THETA)
         WR = 1.0d0
         WI = 0.0d0
         DO M = 1, MMAX, 2
            DO I = M, N, ISTEP
               J = I + MMAX
               TEMPR = WR*DATA( J ) - WI*DATA( J + 1 )
               TEMPI = WR*DATA( J + 1 ) + WI*DATA( J )
               DATA( J ) = DATA( I ) - TEMPR
               DATA( J + 1 ) = DATA( I + 1 ) - TEMPI
               DATA( I ) = DATA( I ) + TEMPR
               DATA( I + 1 ) = DATA( I + 1 ) + TEMPI
            ENDDO
            WTEMP = WR
            WR = WR*WPR - WI*WPI + WR
            WI = WI*WPR + WTEMP*WPI + WI
         ENDDO
         MMAX = ISTEP
         GOTO 502
      ENDIF
C
      IF ( ISIGN.EQ.1 ) THEN
        FAC = 2.0d0/DBLE( NN )
        DO I = 1, NN
          DATA( I ) = FAC*DATA( I )
        ENDDO
        DO I = NN+1, 2*NN
          DATA( I ) = 0.0d0
        ENDDO
      ENDIF
C
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine Solution Vector Radial INTerpolate **********************
C            -        -      -      ---         **********************
C Steve Gibbons Sun Oct 31 14:12:34 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Evaluates the zero^{th} to the (NNDS-1)^{th} derivative of a       C
C solution vector (harmonic IH) at the arbitrary radial point, R.    C
C                                                                    C
C  The double precision array element VALS( i ) is returned          C
C  with f^{i-1}( r ) where f is the function of the IH^{th}          C
C  harmonic in SV.                                                   C
C                                                                    C
C Note that ALL points may be referred to and so if solution         C
C vector contains points which are stored implicitly, these must     C
C be filled in by a call to ASVCPL.                                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NR. Number of radial grid nodes.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C     IH        : Number of radial function (harmonic).              C
C     NNDS      : Number of nodes to be used in interpolation.       C
C                 Must be atleast 2.                                 C
C     IWORK     : Dimension ( NNDS ). Work array.                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C     XARR      : Dim ( * ) but length atleast NR. Location of       C
C                  radial grid nodes.                                C
C                                                                    C
C     VALS      : Returned containing the {i-1}^{th} derivative      C
C                  at r=R in VALS( i )                               C
C                                                                    C
C     WORK      : Dimension ( NNDS ). Work array.                    C
C     COEFM     : Dimension ( NNDS, NNDS ). Work array.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVRINT( R, SV, XARR, INARR, IH, NNDS, VALS, IWORK,
     1                   WORK, COEFM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), NNDS, IWORK( NNDS ), IH
      DOUBLE PRECISION R, SV( * ), XARR( * ), VALS( NNDS ),
     1                 WORK( NNDS ), COEFM( NNDS, NNDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NR, IR, IRAD, IFN, ILN, IND, I, INC, INDFUN, NH
      DOUBLE PRECISION RI, RO, ZERO, ONE, RAD, DTOL
      CHARACTER *(1) TRANS
      PARAMETER ( ZERO = 0.0d0, ONE = 1.0d0, TRANS = 'N',
     1            DTOL = 2.0d-6, INC = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameters ...
C
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C
      IF ( IH.LT.1 .OR. IH.GT.NH ) THEN
        PRINT *,' Subroutine SVRINT.'
        PRINT *,' IH = ',IH
        PRINT *,' NH = ',NH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RI  = XARR(  1 )
      RO  = XARR( NR )
      RAD = R
C
      IF ( RAD.LT.RI .OR. RAD.GT.RO ) THEN
C       .
C       . Ammendment. Here it is possible
C       . that RAD is only less than RI or
C       . greater than RO because of numerical roundoff.
C       . For example:
C             RAD  =   0.6666666666000000    
C             RI   =   0.6666666700000000 
C       . To rectify this, we will try and set a 
C       . new condition on RAD
C       .
        IF ( DABS( RAD - RI ).LT.DTOL ) THEN
          RAD = RI
          GOTO 51
        ENDIF
C       .
        IF ( DABS( RAD - RO ).LT.DTOL ) THEN
          RAD = RO
          GOTO 51
        ENDIF
C       .
        PRINT *,' Subroutine SVRINT.'
        PRINT *,' RAD  = ',RAD
        PRINT *,' RI   = ',RI
        PRINT *,' RO   = ',RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
 51   CONTINUE
C
      IF ( NNDS.LT.2 ) THEN
        PRINT *,' Subroutine SVRINT.'
        PRINT *,' NNDS = ', NNDS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IRAD = NR
      DO IR = 1, NR - 1
        IF ( XARR( IR ).LE.RAD .AND. XARR( IR + 1 ).GT.RAD ) THEN
          IRAD = IR
          GOTO 50
        ENDIF
      ENDDO
C
 50   CONTINUE
      IFN = IRAD + 1 - (NNDS+1)/2
      ILN = IRAD + NNDS/2
C
      DO IR = 1, NNDS
        IWORK( IR ) = IFN + IR - 1
        IF ( IWORK( IR ).LT.1 ) IWORK( IR ) = ILN - IWORK( IR ) + 1
        IF ( IWORK( IR ).GT.NR ) IWORK( IR ) = NR + IFN - IWORK( IR )
        I = IWORK( IR )
        VALS( IR ) = XARR( I )
      ENDDO
C
C Calculate finite difference coefficients
C
      CALL GFDCFD( RAD, VALS, NNDS, COEFM, NNDS, IWORK, WORK )
C
      DO IR = 1, NNDS
        IWORK( IR ) = IFN + IR - 1
        IF ( IWORK( IR ).LT.1 ) IWORK( IR ) = ILN - IWORK( IR ) + 1
        IF ( IWORK( IR ).GT.NR ) IWORK( IR ) = NR + IFN - IWORK( IR )
        I = IWORK( IR )
        IND = INDFUN( I, IH, INARR )
        WORK( IR ) = SV( IND )
      ENDDO
C
C WORK now contains the function values at the nodes
C so we can now multiply COEFM by WORK to get VALS
C
      CALL DGEMV ( TRANS, NNDS, NNDS, ONE, COEFM, NNDS, WORK, INC,
     1             ZERO, VALS, INC )
C  
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine MATrix INDices ******************************************
C            ---    ---     ******************************************
C Steve Gibbons                                                      C
C____________________________________________________________________C
C                                                                    C
C  Tues. Nov 3rd. 1998                                               C
C  -------------------                                               C
C In order that subroutines which write matrices can be adapted to   C
C write in different storage formats, this routine takes IR and IC   C
C (the ACTUAL coordinates of the 'physical' matrix) and calculates   C
C IROW and ICOL (the row and column number in which the number is    C
C to be placed). Hence, the same code will be able to write to a     C
C matrix of any format.                                              C
C  The code also checks that the bounds of the matrix have not       C
C been exceeded.                                                     C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Row from theoretical matrix.                       C
C     IC        : Column from theoretical matrix.                    C
C     IMF       : Matrix format flag.                                C
C                                                                    C
C          imf = 1; Matrix is in LAPACK banded format                C
C                   ie element a_{i,j} is stored in                  C
C                   A( kle + ku + 1 + i - j , j )                    C
C                                                                    C
C          imf = 2; Matrix is banded but with element a_{i,j}        C
C                   stored in A( kl + 1 + j - i , i ). The matrices  C
C                   in Dave Gubbins' dynamo codes are in this        C
C                   format.                                          C
C                                                                    C
C          imf = 3; Matrix is square - ie a_{i,j} is stored          C
C                   in A( i, j ).                                    C
C                                                                    C
C     KL        : Number of lower diagonals in band matrix.          C
C     KU        : Number of upper diagonals in band matrix.          C
C     KLE       : Number of additional lower diagonals               C
C                  in band matrix. This should be set to zero for    C
C                  most applications but should be equal to KL if    C
C                  the matrix is to be LU decomposed by a LAPACK or  C
C                  NAG subroutine.                                   C
C                                                                    C
C  ( Note - KL, KU and KLE are not referenced if IMF is set to 3.)   C
C                                                                    C
C     N1        : Leading dimension of matrix A.                     C
C     N2        : Second dimension of matrix A.                      C
C                                                                    C
C       ( Note - i.e. A is dimensioned A( N1, N2 ) ... )             C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Integer                                                           C
C  -------                                                           C
C     IROW      : 'i' coordinate in the final matrix.                C
C     ICOL      : 'j' coordinate in the final matrix.                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MATIND (IR,IC,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR,IC,IMF,KL,KU,KLE,N1,N2,IROW,ICOL
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NJ
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C  First check for errors on the definition of IMF
C
      IF ( IMF.NE.1 .AND. IMF.NE.2 .AND. IMF.NE.3 ) THEN
         PRINT *,' Subroutine MATIND '
         PRINT *,' IMF = ', IMF
         PRINT *,' Illegal value. Program aborted.'
         STOP
      ENDIF
C
C  Check that N1 is equal to (KLE + KL + KU + 1) for
C  the banded matrices.
C
      IF ( IMF.EQ.1 .OR. IMF.EQ.2 ) THEN
         NJ = KLE + KL + KU + 1
         IF ( N1.NE.NJ ) THEN
           PRINT *,' Subroutine MATIND '
           PRINT *,' KL  = ', KL 
           PRINT *,' KU  = ', KU 
           PRINT *,' KLE = ', KLE
           PRINT *,' N1  = ', N1 
           PRINT *,' Illegal value. Program aborted.'
           STOP
         ENDIF
      ENDIF
C
C  Do case IMF = 1. So this is a LAPACK format band matrix.
C
      IF ( IMF.EQ.1 ) THEN
         ICOL = IC
         IROW = KLE + KU + 1 + IR - IC
      ENDIF
C
C  Do case IMF = 2. So this is a DG format band matrix.
C
      IF ( IMF.EQ.2 ) THEN
         ICOL = IR
         IROW = KL + 1 + IC - IR
      ENDIF
C
C  Do case IMF = 3. So this is a square matrix.
C
      IF ( IMF.EQ.3 ) THEN
         ICOL = IC
         IROW = IR
      ENDIF
C
C  Final error check on IROW and ICOL
C
      IF ( IROW.LT.1 .OR. IROW.GT.N1 ) THEN
         PRINT *,' Subroutine MATIND '
         PRINT *,' IMF  = ',IMF
         PRINT *,' IROW = ', IROW
         PRINT *,' 1st. dimension = ', N1
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( ICOL.LT.1 .OR. ICOL.GT.N2 ) THEN
         PRINT *,' Subroutine MATIND '
         PRINT *,' IMF  = ',IMF
         PRINT *,' ICOL = ', ICOL
         PRINT *,' 2nd. dimension = ', N2
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Banded Matrix Row or Column OPeration *******************
C            -      -      -      -      --        *******************
C Steve Gibbons 3.12.97                                              C
C____________________________________________________________________C
C                                                                    C
C If a double precision matrix ABAND is stored in LAPACK format      C
C i.e. with the element a(i,j) being stored in                       C
C  ABAND( KLE + KU + 1 + i - j, j ),                                 C
C then BMRCOP will take either a row or a column of A and            C
C multiply it by ALPHA and add BETA.                                 C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C     NDIM	: Dimension of matrix.                               C
C     IRC	: = 1 to act on ROW                                  C
C                 = 2 to act on COLUMN                               C
C     NRC	: Number of row or column.                           C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     ABAND     : Banded matrix. Dimensions ( NDIM1, NDIM2 ) with    C
C                  NDIM1 = KL + KU + KLE + 1                         C
C                  NDIM2 = Dimension of problem (NDIM ).             C
C     ALPHA	: First multiplies current matrix entry by ALPHA     C
C     BETA 	: Then adds on BETA.                                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE BMRCOP( KL, KU, KLE, NDIM, IRC, NRC, ABAND, 
     1                   ALPHA, BETA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KL, KU, KLE, NDIM, IRC, NRC
      DOUBLE PRECISION ABAND( KLE + KU + KL + 1, NDIM ) , ALPHA, BETA
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IC, IR, IROW, IND1, IND2
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
C     .                    First deal with the case acting on ROWS ...
      IF ( IRC.EQ.1 ) THEN
         IR = NRC
         IND1 = MAX( 1, IR - KL )
         IND2 = MIN( IR + KU, NDIM )
         DO IC = IND1, IND2
            IROW = KLE + KU + 1 + IR - IC
            ABAND( IROW, IC ) = ABAND( IROW, IC )*ALPHA + BETA
         ENDDO
         RETURN
      ENDIF
c
C     .                    First deal with the case acting on COLUMNS.
      IF ( IRC.EQ.2 ) THEN
         IC = NRC
         IND1 = MAX( 1, IC - KU )
         IND2 = MIN( IC + KL, NDIM )
         DO IR = IND1, IND2
            IROW = KLE + KU + 1 + IR - IC
            ABAND( IROW, IC ) = ABAND( IROW, IC )*ALPHA + BETA
         ENDDO
         RETURN
      ENDIF
C
      PRINT *,' Subroutine BMRCOP '
      PRINT *,' IRC must be set to either 1 or 2.'
      PRINT *,' Program aborted.'
      STOP
C
      END
C*********************************************************************
C*********************************************************************
C subroutine POWer of TWO checker ************************************
C            ---      ---         ************************************
C Steve Gibbons 11.4.97        					     C
C____________________________________________________________________C
C Checks an integer INPUT and returns a .TRUE. logical variable POT  C
C if and only if INPUT is a power of 2.				     C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     INPUT	: Integer to be tested				     C
C____________________________________________________________________C
C Output :-                                                          C
C ======                                                             C
C  Logical 							     C
C  ------- 							     C
C     POT	: Power Of Two ? 				     C
C____________________________________________________________________C
C*********************************************************************
      SUBROUTINE POWTWO(INPUT,POT)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INPUT
      LOGICAL POT
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITWO,N,IREM
      PARAMETER (ITWO=2)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C First of all put N = INPUT so that INPUT won't be altered
      N=INPUT
C Then check that N is atleast equal to 2 - otherwise
C it is obviously not a power of 2.
      IF ( N.LT.2 ) THEN
         POT=.FALSE.
         RETURN
      ENDIF
C Begin our iteration ....
 500  CONTINUE
      IREM = MOD ( N , ITWO )
      IF ( IREM.NE.0 ) THEN
         POT=.FALSE.
         RETURN
      ENDIF
      IF ( N.EQ.2 ) THEN
         POT=.TRUE.
         RETURN
      ELSE
         N=N/2
         GOTO 500
      ENDIF
      END
C*********************************************************************
C*********************************************************************
C function SQuare Root of L*(L+1) ************************************
C          --     -       -  - -  ************************************
C Steve Gibbons 25.4.97                                              C
C____________________________________________________________________C
C
      FUNCTION SQRLL1 ( L )
      IMPLICIT NONE
      DOUBLE PRECISION SQRLL1,Q
      INTEGER L
   
      IF ( L.LT.0 ) THEN
         PRINT *,' Function SQRLL1. L less than 0. Program aborted.'
         STOP
      ENDIF
      Q = DBLE( L )
      SQRLL1 = DSQRT( Q*Q + Q )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function Element of Matrix MULTiplication evaluate *****************
C          -          -      ----                    *****************
C Steve Gibbons Fri Oct 22 17:23:46 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Let the matrices B and C be respectively 'm by k' and 'k by n'     C
C double precision matrices. Then A = B C is an 'm by n' matrix      C
C with a_{ij} = \sum_{l=1,k} b_{il} c_{lj}                           C
C                                                                    C
C EMMULT returns the value of A( I, J)                               C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     I         : Row containing desired element of A.               C
C     J         : Column containing desired element of A.            C
C     LDB       : Leading dimension of matrix B.                     C
C     LDC       : Leading dimension of matrix C.                     C
C     K         : Number of columns of matrix B and                  C
C                 Number of rows of matrix C.                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     B         : First matrix. Dim ( LDB, * )                       C
C     C         : First matrix. Dim ( LDC, * )                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION EMMULT( I, J, LDB, LDC, K, B, C )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER I, J, LDB, LDC, K
      DOUBLE PRECISION EMMULT, B( LDB, * ), C( LDB, * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First check the validity of the integer inputs.
C K must not exceed LDC
C I must not exceed LDB
C     .
      IF ( K.GT.LDC .OR. I.GT.LDB ) THEN
        PRINT *,' Function EMMULT.'
        PRINT *,' K    = ', K,' I    = ', I
        PRINT *,' LDB  = ', LDB,' LDC  = ', LDC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      EMMULT = 0.0d0
C     . 
      DO L = 1, K
        EMMULT = EMMULT + B( I, L )*C( L, J )
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

C*********************************************************************
C function PMM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Gives the Schmidt Normalised Legendre Function P_m^m ( X )         C
C from eqn. 175 in my notes ie                                       C
C                                                                    C
C                   ( 2m - 1)!! (1- XX)^(m/2) * SQRT (2.0d0 )        C
C   P_m^m( X ) =  ---------------------------------------------      C
C                        SQRT ( (2m)! )                              C
C                                                                    C
C       for m non-zero and                                           C
C                                                                    C
C   P_0^0( X ) = 1.0d0                                               C
C                                                                    C
C N.B. The double factorial sign means the product of all ODD        C
C integers between 1 and ( 2m - 1 ).                                 C
C Best to use the form in eq 184 to calculate PMM.                   C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SINE      : Sin (theta)                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMM ( M, S )
      IMPLICIT NONE
      DOUBLE PRECISION PMM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION S
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
      DOUBLE PRECISION TWOI
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM. M < 0 error. Stopped.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         PMM = 1.0d0
         RETURN
      ENDIF
C ................................. so M is greater than 0
      PMM = DSQRT ( 2.0d0 )
      DO I = 1, M
         TWOI = 2.0d0*DBLE(I)
         PMM = PMM * S * DSQRT ( (TWOI - 1.0d0 )/TWOI )
      ENDDO
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function PMM1 ******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Evaluates the Schmidt Normalised Legendre Function P_(m+1)^m (X)   C
C according to equation 179 in my notes ; i.e.                       C
C                                                                    C
C    P_(m+1)^m (X) = SQRT( 2m+1 ).X.P^m_m(X)                         C
C                                                                    C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     PMM0      : P_m^m(X) as evaluated by Function PMM.             C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMM1 ( M, X, PMM0 )
      IMPLICIT NONE
      DOUBLE PRECISION PMM1
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION X, PMM0
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = DBLE( M )
      PMM1 = X*PMM0*DSQRT( 2.0d0*RM+1.0d0 )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function PLM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Calculates the Schmidt Normalised Legendre Function P_l^m (x)      C
C given P_(l-1)^m and P_(l-2)^m according to equation 183 in my notesC
C i.e.                                                               C
C   P_l^m( X ) = { A * P_(l-1)^m - B * P_(l-2)^m }/C                 C
C                                                                    C
C where A = (2*l - 1)*X ,                                            C
C                                                                    C
C B = SQRT( (L+M-1)*(L-M-1) ) and C = SQRT( (L+M)*(L-M) )            C
C                                                                    C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C     L         : Well it's L isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     PLMIN1    : P_(l-1)^m ( X )                                    C
C     PLMIN2    : P_(l-2)^m ( X )                                    C
C____________________________________________________________________C
C Local Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     RL        : DBLE  ( L )                                        C
C     RM        : DBLE  ( M )                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PLM ( L, M, X, PLMIN1, PLMIN2 )
      IMPLICIT NONE
      DOUBLE PRECISION PLM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M,L
      DOUBLE PRECISION X, PLMIN1, PLMIN2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM,RL
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function PLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function PLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' PLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = DBLE( M )
      RL = DBLE( L )
      PLM = ( 2.0d0*RL - 1.0d0 )*X*PLMIN1
      PLM = PLM - PLMIN2*DSQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      PLM = PLM/DSQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function DPMM ******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Calculates the derivative of P_m^m(theta) according to equation    C
C 185 in my notes i.e.                                               C
C                                                                    C
C d P_m^m(theta)/ d(theta) = sqrt(2.0)*M*cos(theta)/sin(theta) * A   C
C  with A = product_(i=1)^m ( SQRT( (2i-1)/2i ) * sin(theta) )       C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SINE      : Sin (theta)                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPMM ( M , C , S )
      IMPLICIT NONE
      DOUBLE PRECISION DPMM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION C, S
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
      DOUBLE PRECISION TWOI
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function DPMM. M < 0 error. Stopped.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         DPMM = 0.0d0
         RETURN
      ENDIF
C ................................. so M is greater than 0
      DPMM = DSQRT ( 2.0d0 )*DBLE(M)*C/S
      DO I = 1, M
         TWOI = 2.0d0*DBLE(I)
         DPMM = DPMM * S * DSQRT ( (TWOI - 1.0d0 )/TWOI )
      ENDDO
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function DPMM1 *****************************************************
C Steve Gibbons 18.4.97                                              C
C____________________________________________________________________C
C Calculates d(P_(m+1)^m/d(theta) by equation 185 in my notes        C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     S         : Sin(theta)                                         C
C     PMM0      : P_m^m(X) as evaluated by Function PMM.             C
C     DPMM0     : d(P_m^m(X))/d(theta) as evaluated by Function DPMM C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPMM1 ( M, X, S, PMM0, DPMM0)
      IMPLICIT NONE
      DOUBLE PRECISION DPMM1
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION X, S, PMM0, DPMM0
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function DPMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = DBLE( M )
      DPMM1 = DSQRT( 2.0d0*RM+1.0d0 )*( X*DPMM0 - S*PMM0 )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function DPLM *******************************************************
C Steve Gibbons 18.4.97                                              C
C____________________________________________________________________C
C                                                                    C
C Calculates general P_l^m derivative from eqn 187 in my notes       C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C     L         : Well it's L isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     S         : Sin(theta)                                         C
C     PLMIN1    : P_(l-1)^m ( X )                                    C
C     DPLMN1   : P_(l-1)^m ( X ) derivative                         C
C     DPLMN2   : P_(l-2)^m ( X ) derivative                         C
C____________________________________________________________________C
C Local Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     RL        : DBLE  ( L )                                        C
C     RM        : DBLE  ( M )                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPLM ( L, M, X, S, PLMIN1, DPLMN1, DPLMN2 )
      IMPLICIT NONE
      DOUBLE PRECISION DPLM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M,L
      DOUBLE PRECISION X, S, PLMIN1, DPLMN1, DPLMN2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM,RL
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function DPLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function DPLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' DPLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = DBLE( M )
      RL = DBLE( L )
      DPLM = ( 2.0d0*RL - 1.0d0 )*(X*DPLMN1-S*PLMIN1)
      DPLM = DPLM - DPLMN2*DSQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      DPLM = DPLM/DSQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function DL ********************************************************
C          -- ********************************************************
C Steve Gibbons 8.4.97						     C
C____________________________________________________________________C
C Spherical Harmonic Differential Operator Dl			     C
C 								     C
C Dl f = f'' + 2/r f' - l(l+1)/(rr).f				     C
C____________________________________________________________________C
C Input variables :-						     C
C ===============						     C
C  Double Precision						     C
C  ----------------						     C
C     R		: Radius					     C
C     F		: Value of function at r(i)			     C
C     DF	: Value of first derivative at r(i)		     C
C     DDF	: Value of second derivative at r(i)		     C
C  Integer						             C
C  -------						             C
C     L		: Spherical harmonic degree			     C
C____________________________________________________________________C
C Output :-						    	     C
C ======						  	     C
C  Double Precision						     C
C  ----------------						     C
C     DL	: as above					     C
C____________________________________________________________________C
C*********************************************************************
      FUNCTION DL(L,R,F,DF,DDF)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION DL,R,F,DF,DDF
      INTEGER L
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION FACT,TOL
      PARAMETER (TOL=1.0d-10)
C Note that both of these working variables are confirmed to
C be explicitly set sjg 15.3.97 .
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Avoid division by zero ....
      IF (ABS(R).LT.TOL) THEN
         WRITE (6,987)
         WRITE (6,988)
         WRITE (6,989)
 987     FORMAT ('Function DL has received too small a ')
 988     FORMAT ('value of R. This would result in a ')
 989     FORMAT ('division by zero error. Program aborted. ')
         STOP
      ENDIF
C
      FACT=DBLE(L*L+L)
      DL=DDF+2.0d0*DF/R-FACT*F/(R*R)
      RETURN
      END
C*********************************************************************
C*********************************************************************
C integer function INDex FUNction ************************************
C                  ---   ---      ************************************
C Steve Gibbons Thu Sep 16 11:13:38 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Returns the location in the solution vector of the value of the    C
C IR^{th} grid node of the IH^{th} harmonic.                         C
C                                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Number of radial grid node.                        C
C     IH        : Number of harmonic.                                C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1/3. INDFUN = ( IR - 1 )*NH + IH        C
C                   IFORMF = 2/4. INDFUN = ( IH - 1 )*NR + IR        C
C                                                                    C
C  Here, NH is the TOTAL number of harmonics in the solution         C
C  vector - or atleast the part of which is visible to that part     C
C  of the program. NR is the number of radial grid nodes             C
C  corresponding to that harmonic (which for cases IFORMF = 1 and    C
C  IFORMF = 2 is identical for all harmonics - more complicated      C
C  options (for example magnetic fields which have to be resolved    C
C  beyond the region of fluid flow) may be added later and this      C
C  routine should be flexible to all possibilities with extra        C
C  constraints being added in other elements of INARR.               C
C                                                                    C
C  IFORMF = 1/3 is the option likely to be used for the purpose of   C
C  solution as it allows the banded formation of a matrix.           C
C                                                                    C
C  IFORMF = 2/4 is the option likely to be used for the purpose of   C
C  displaying solutions as it stores adjacent nodes for each         C
C  harmonic together.                                                C
C                                                                    C
C                 INARR( 2 ) = NR. See above.                        C
C                 INARR( 3 ) = NH. See above.                        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION INDFUN ( IR, IH, INARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INDFUN, IR, IH, INARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFORMF, NR, NH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFORMF = INARR( 1 )
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C
      IF ( IR.LT.1 .OR. IR.GT.NR ) THEN
        PRINT *,' Function INDFUN. IR = ', IR
        PRINT *,' NR = ', NR,'. Program aborted.'
        STOP
      ENDIF
C
      IF ( IH.LT.1 .OR. IH.GT.NH ) THEN
        PRINT *,' Function INDFUN. IH = ', IH
        PRINT *,' NH = ', NH,'. Program aborted.'
        STOP
      ENDIF
C
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.3 ) THEN
        INDFUN = ( IR - 1 )*NH + IH
        RETURN
      ENDIF
C
      IF ( IFORMF.EQ.2 .OR. IFORMF.EQ.4 ) THEN
        INDFUN = ( IH - 1 )*NR + IR
        RETURN
      ENDIF
C
      PRINT *,' Function INDFUN. IFORMF = ', IFORMF
      PRINT *,' Not current option. Program aborted.'
      STOP
      END
C*********************************************************************
      DOUBLE PRECISION FUNCTION DNRM2 ( N, X, INCX )
*     .. Scalar Arguments ..
      INTEGER                           INCX, N
*     .. Array Arguments ..
      DOUBLE PRECISION                  X( * )
*     ..
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*
*     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      INTEGER               IX
      DOUBLE PRECISION      ABSXI, NORM, SCALE, SSQ
*     .. Intrinsic Functions ..
      INTRINSIC             ABS, SQRT
*     ..
*     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE IF( N.EQ.1 )THEN
         NORM  = ABS( X( 1 ) )
      ELSE
         SCALE = ZERO
         SSQ   = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( X( IX ).NE.ZERO )THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI )THEN
                  SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SSQ   = SSQ   +     ( ABSXI/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF
*
      DNRM2 = NORM
      RETURN
*
*     End of DNRM2.
*
      END
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
      SUBROUTINE DGBMV ( TRANS, M, N, KL, KU, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, KL, KU, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGBMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  KL     - INTEGER.
*           On entry, KL specifies the number of sub-diagonals of the
*           matrix A. KL must satisfy  0 .le. KL.
*           Unchanged on exit.
*
*  KU     - INTEGER.
*           On entry, KU specifies the number of super-diagonals of the
*           matrix A. KU must satisfy  0 .le. KU.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading ( kl + ku + 1 ) by n part of the
*           array A must contain the matrix of coefficients, supplied
*           column by column, with the leading diagonal of the matrix in
*           row ( ku + 1 ) of the array, the first super-diagonal
*           starting at position 2 in row ku, the first sub-diagonal
*           starting at position 1 in row ( ku + 2 ), and so on.
*           Elements in the array A that do not correspond to elements
*           in the band matrix (such as the top left ku by ku triangle)
*           are not referenced.
*           The following program segment will transfer a band matrix
*           from conventional full matrix storage to band storage:
*
*                 DO 20, J = 1, N
*                    K = KU + 1 - J
*                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
*                       A( K + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( kl + ku + 1 ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry, the incremented array Y must contain the
*           vector y. On exit, Y is overwritten by the updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, K, KUP1, KX, KY,
     $                   LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( KL.LT.0 )THEN
         INFO = 4
      ELSE IF( KU.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.( KL + KU + 1 ) )THEN
         INFO = 8
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 10
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGBMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the band part of A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      KUP1 = KU + 1
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  K    = KUP1 - J
                  DO 50, I = MAX( 1, J - KU ), MIN( M, J + KL )
                     Y( I ) = Y( I ) + TEMP*A( K + I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  K    = KUP1 - J
                  DO 70, I = MAX( 1, J - KU ), MIN( M, J + KL )
                     Y( IY ) = Y( IY ) + TEMP*A( K + I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
               IF( J.GT.KU )
     $            KY = KY + INCY
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               K    = KUP1 - J
               DO 90, I = MAX( 1, J - KU ), MIN( M, J + KL )
                  TEMP = TEMP + A( K + I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               K    = KUP1 - J
               DO 110, I = MAX( 1, J - KU ), MIN( M, J + KL )
                  TEMP = TEMP + A( K + I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
               IF( J.GT.KU )
     $            KX = KX + INCX
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGBMV .
*
      END
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DGER  .
*
      END
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      SUBROUTINE DTBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, K, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTBSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular band matrix, with ( k + 1 )
*  diagonals.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows:
*
*              TRANS = 'N' or 'n'   A*x = b.
*
*              TRANS = 'T' or 't'   A'*x = b.
*
*              TRANS = 'C' or 'c'   A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with UPLO = 'U' or 'u', K specifies the number of
*           super-diagonals of the matrix A.
*           On entry with UPLO = 'L' or 'l', K specifies the number of
*           sub-diagonals of the matrix A.
*           K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           The following program segment will transfer an upper
*           triangular band matrix from conventional full matrix storage
*           to band storage:
*
*                 DO 20, J = 1, N
*                    M = K + 1 - J
*                    DO 10, I = MAX( 1, J - K ), J
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the
*           array A is not referenced.
*           The following program segment will transfer a lower
*           triangular band matrix from conventional full matrix storage
*           to band storage:
*
*                 DO 20, J = 1, N
*                    M = 1 - J
*                    DO 10, I = J, MIN( N, J + K )
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Note that when DIAG = 'U' or 'u' the elements of the array A
*           corresponding to the diagonal elements of the matrix are not
*           referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KPLUS1, KX, L
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( K.LT.0 )THEN
         INFO = 5
      ELSE IF( LDA.LT.( K + 1 ) )THEN
         INFO = 7
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTBSV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed by sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := inv( A )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     L = KPLUS1 - J
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( KPLUS1, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, MAX( 1, J - K ), -1
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 40, J = N, 1, -1
                  KX = KX - INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = KPLUS1 - J
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( KPLUS1, J )
                     TEMP = X( JX )
                     DO 30, I = J - 1, MAX( 1, J - K ), -1
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      - INCX
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     L = 1 - J
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( 1, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, MIN( N, J + K )
                        X( I ) = X( I ) - TEMP*A( L + I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  KX = KX + INCX
                  IF( X( JX ).NE.ZERO )THEN
                     IX = KX
                     L  = 1  - J
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( 1, J )
                     TEMP = X( JX )
                     DO 70, I = J + 1, MIN( N, J + K )
                        X( IX ) = X( IX ) - TEMP*A( L + I, J )
                        IX      = IX      + INCX
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := inv( A')*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            KPLUS1 = K + 1
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  L    = KPLUS1 - J
                  DO 90, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( KPLUS1, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  L    = KPLUS1  - J
                  DO 110, I = MAX( 1, J - K ), J - 1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( KPLUS1, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
                  IF( J.GT.K )
     $               KX = KX + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  L    = 1      - J
                  DO 130, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( 1, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  L    = 1       - J
                  DO 150, I = MIN( N, J + K ), J + 1, -1
                     TEMP = TEMP - A( L + I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( 1, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
                  IF( ( N - J ).GE.K )
     $               KX = KX - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTBSV .
*
      END
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*inv( A' )*B.
*
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*inv( A' ).
*
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END
      SUBROUTINE DTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*A*B.
*
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'.
*
            IF( UPPER )THEN
               DO 110, J = 1, N
                  DO 100, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
            ELSE
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 120, K = I + 1, M
                        TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*A.
*
            IF( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     END IF
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
  220          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'.
*
            IF( UPPER )THEN
               DO 260, K = 1, N
                  DO 240, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  END IF
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     END IF
  280             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMM .
*
      END
      SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,
*
*  where x is an n element vector and  A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := A'*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := A*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := A'*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMV .
*
      END
      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
     $                   INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGBTRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general band matrix A using the LU factorization computed
*  by DGBTRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*          Details of the LU factorization of the band matrix A, as
*          computed by DGBTRF.  U is stored as an upper triangular band
*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
*          the multipliers used during the factorization are stored in
*          rows KL+KU+2 to 2*KL+KU+1.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= N, row i of the matrix was
*          interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSWAP, DTBSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      KD = KU + KL + 1
      LNOTI = KL.GT.0
*
      IF( NOTRAN ) THEN
*
*        Solve  A*X = B.
*
*        Solve L*X = B, overwriting B with X.
*
*        L is represented as a product of permutations and unit lower
*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
*        where each transformation L(i) is a rank-one modification of
*        the identity matrix.
*
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J )
     $            CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL DGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ),
     $                    LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
*
         DO 20 I = 1, NRHS
*
*           Solve U*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU,
     $                  AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
*
      ELSE
*
*        Solve A'*X = B.
*
         DO 30 I = 1, NRHS
*
*           Solve U'*X = B, overwriting B with X.
*
            CALL DTBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB,
     $                  LDAB, B( 1, I ), 1 )
   30    CONTINUE
*
*        Solve L'*X = B, overwriting B with X.
*
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL DGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ),
     $                     LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L.NE.J )
     $            CALL DSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of DGBTRS
*
      END
      SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGBTRF computes an LU factorization of a real m-by-n band matrix A
*  using partial pivoting with row interchanges.
*
*  This is the blocked version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  M = N = 6, KL = 2, KU = 1:
*
*  On entry:                       On exit:
*
*      *    *    *    +    +    +       *    *    *   u14  u25  u36
*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*
*  Array elements marked * are not used by the routine; elements marked
*  + need not be set on entry, but are required by the routine to store
*  elements of U because of fill-in resulting from the row interchanges.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,
     $                   JU, K2, KM, KV, NB, NW
      DOUBLE PRECISION   TEMP
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   WORK13( LDWORK, NBMAX ),
     $                   WORK31( LDWORK, NBMAX )
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX, ILAENV
      EXTERNAL           IDAMAX, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGBTF2, DGEMM, DGER, DLASWP, DSCAL,
     $                   DSWAP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     KV is the number of superdiagonals in the factor U, allowing for
*     fill-in
*
      KV = KU + KL
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment
*
      NB = ILAENV( 1, 'DGBTRF', ' ', M, N, KL, KU )
*
*     The block size must not exceed the limit set by the size of the
*     local arrays WORK13 and WORK31.
*
      NB = MIN( NB, NBMAX )
*
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
*
*        Use unblocked code
*
         CALL DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
*
*        Use blocked code
*
*        Zero the superdiagonal elements of the work array WORK13
*
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
*
*        Zero the subdiagonal elements of the work array WORK31
*
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
*
*        Gaussian elimination with partial pivoting
*
*        Set fill-in elements in columns KU+2 to KV to zero
*
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
*
*        JU is the index of the last column affected by the current
*        stage of the factorization
*
         JU = 1
*
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
*
*           The active part of the matrix is partitioned
*
*              A11   A12   A13
*              A21   A22   A23
*              A31   A32   A33
*
*           Here A11, A21 and A31 denote the current block of JB columns
*           which is about to be factorized. The number of rows in the
*           partitioning are JB, I2, I3 respectively, and the numbers
*           of columns are JB, J2, J3. The superdiagonal elements of A13
*           and the subdiagonal elements of A31 lie outside the band.
*
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
*
*           J2 and J3 are computed after JU has been updated.
*
*           Factorize the current block of JB columns
*
            DO 80 JJ = J, J + JB - 1
*
*              Set fill-in elements in column JJ+KV to zero
*
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
*
*              Find pivot and test for singularity. KM is the number of
*              subdiagonal elements in the current column.
*
               KM = MIN( KL, M-JJ )
               JP = IDAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
*
*                    Apply interchange to columns J to J+JB-1
*
                     IF( JP+JJ-1.LT.J+KL ) THEN
*
                        CALL DSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1,
     $                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
*
*                       The interchange affects columns J to JJ-1 of A31
*                       which are stored in the work array WORK31
*
                        CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL DSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,
     $                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
*
*                 Compute multipliers
*
                  CALL DSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ),
     $                        1 )
*
*                 Update trailing submatrix within the band and within
*                 the current block. JM is the index of the last column
*                 which needs to be updated.
*
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )
     $               CALL DGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,
     $                          AB( KV, JJ+1 ), LDAB-1,
     $                          AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
*
*                 If pivot is zero, set INFO to the index of the pivot
*                 unless a zero pivot has already been found.
*
                  IF( INFO.EQ.0 )
     $               INFO = JJ
               END IF
*
*              Copy current column of A31 into the work array WORK31
*
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )
     $            CALL DCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,
     $                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
*
*              Apply the row interchanges to the other blocks.
*
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
*
*              Use DLASWP to apply the row interchanges to A12, A22, and
*              A32.
*
               CALL DLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,
     $                      IPIV( J ), 1 )
*
*              Adjust the pivot indices.
*
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
*
*              Apply the row interchanges to A13, A23, and A33
*              columnwise.
*
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
*
*              Update the relevant part of the trailing submatrix
*
               IF( J2.GT.0 ) THEN
*
*                 Update A12
*
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,
     $                        AB( KV+1-JB, J+JB ), LDAB-1 )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A22
*
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J2,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Update A32
*
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J2,
     $                           JB, -ONE, WORK31, LDWORK,
     $                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,
     $                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
*
               IF( J3.GT.0 ) THEN
*
*                 Copy the lower triangle of A13 into the work array
*                 WORK13
*
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
*
*                 Update A13 in the work array
*
                  CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit',
     $                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,
     $                        WORK13, LDWORK )
*
                  IF( I2.GT.0 ) THEN
*
*                    Update A23
*
                     CALL DGEMM( 'No transpose', 'No transpose', I2, J3,
     $                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,
     $                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ),
     $                           LDAB-1 )
                  END IF
*
                  IF( I3.GT.0 ) THEN
*
*                    Update A33
*
                     CALL DGEMM( 'No transpose', 'No transpose', I3, J3,
     $                           JB, -ONE, WORK31, LDWORK, WORK13,
     $                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
*
*                 Copy the lower triangle of A13 back into place
*
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
*
*              Adjust the pivot indices.
*
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
*
*           Partially undo the interchanges in the current block to
*           restore the upper triangular form of A31 and copy the upper
*           triangle of A31 back into place
*
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
*
*                 Apply interchange to columns J to JJ-1
*
                  IF( JP+JJ-1.LT.J+KL ) THEN
*
*                    The interchange does not affect A31
*
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
*
*                    The interchange does affect A31
*
                     CALL DSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,
     $                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
*
*              Copy the current column of A31 back into place
*
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )
     $            CALL DCOPY( NW, WORK31( 1, JJ-J+1 ), 1,
     $                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
*
      RETURN
*
*     End of DGBTRF
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DGETRI computes the inverse of a matrix using the LU factorization
*  computed by DGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by DGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, NB, NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      WORK( 1 ) = MAX( N, 1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DGETRI
*
      END
      SUBROUTINE DGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGBTF2 computes an LU factorization of a real m-by-n band matrix A
*  using partial pivoting with row interchanges.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  Further Details
*  ===============
*
*  The band storage scheme is illustrated by the following example, when
*  M = N = 6, KL = 2, KU = 1:
*
*  On entry:                       On exit:
*
*      *    *    *    +    +    +       *    *    *   u14  u25  u36
*      *    *    +    +    +    +       *    *   u13  u24  u35  u46
*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
*
*  Array elements marked * are not used by the routine; elements marked
*  + need not be set on entry, but are required by the routine to store
*  elements of U, because of fill-in resulting from the row
*  interchanges.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     KV is the number of superdiagonals in the factor U, allowing for
*     fill-in.
*
      KV = KU + KL
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGBTF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Gaussian elimination with partial pivoting
*
*     Set fill-in elements in columns KU+2 to KV to zero.
*
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
*
*     JU is the index of the last column affected by the current stage
*     of the factorization.
*
      JU = 1
*
      DO 40 J = 1, MIN( M, N )
*
*        Set fill-in elements in column J+KV to zero.
*
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
*
*        Find pivot and test for singularity. KM is the number of
*        subdiagonal elements in the current column.
*
         KM = MIN( KL, M-J )
         JP = IDAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
*
*           Apply interchange to columns J to JU.
*
            IF( JP.NE.1 )
     $         CALL DSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1,
     $                     AB( KV+1, J ), LDAB-1 )
*
            IF( KM.GT.0 ) THEN
*
*              Compute multipliers.
*
               CALL DSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
*
*              Update trailing submatrix within the band.
*
               IF( JU.GT.J )
     $            CALL DGER( KM, JU-J, -ONE, AB( KV+2, J ), 1,
     $                       AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ),
     $                       LDAB-1 )
            END IF
         ELSE
*
*           If pivot is zero, set INFO to the index of the pivot
*           unless a zero pivot has already been found.
*
            IF( INFO.EQ.0 )
     $         INFO = J
         END IF
   40 CONTINUE
      RETURN
*
*     End of DGBTF2
*
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IP, IX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSWAP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.EQ.0 )
     $   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, JP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTRI computes the inverse of a real upper or lower triangular
*  matrix A.
*
*  This is the Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*               matrix is singular and its inverse can not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DTRTRI
*
      END
      SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTI2 computes the inverse of a real upper or lower triangular
*  matrix.
*
*  This is the Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading n by n upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DTRMV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DTRTI2
*
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
*     End of ILAENV
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
