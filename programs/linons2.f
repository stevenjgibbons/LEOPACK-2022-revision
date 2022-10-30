C*********************************************************************
C                                                                    C
C Steve Gibbons - Linear Onset 2.                                    C
C Wed Nov 24 13:01:17 GMT 1999                                       C
C                                                                    C
C                                                                    C
C                                                                    C
C                                                                    C
C*********************************************************************
      PROGRAM linons2
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, LHMAX, NHMAX, NPHMAX, NTHMAX, KLMAX,
     1        NBNDMX, NDCS, NDRVM, ISVMAX, NPMAX, LHLH2M, NCFM,
     2        NBN, NCVM, NRUNM
C
      PARAMETER ( NRMAX = 150, LHMAX = 124, NHMAX = 150, NBN = 3,
     1            NTHMAX = 128, NPHMAX = 256, KLMAX = (NBN+1)*NHMAX-1,
     2            NBNDMX = 3*KLMAX+1, NDCS = 3, NDRVM = 4, 
     3            ISVMAX = NRMAX*NHMAX )
      PARAMETER ( NPMAX = (LHMAX+1)*(LHMAX+2)/2, 
     1            LHLH2M = LHMAX*(LHMAX+2), NCFM = 2*NBN + 1,
     2            NCVM = 24, NRUNM = 100 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER MHIBC( 3 ), MHOBC( 3 ), LARR( 3 ), INARR( 3 ),
     1        IPIV( ISVMAX ), MHT( NHMAX ), MHL( NHMAX ),
     2        MHM( NHMAX ), MHP( NHMAX ), MHTR( NHMAX ),
     3        IWORK( NCFM )
C
      INTEGER          NRARR( NRUNM )
      INTEGER          ISPARR( NRUNM )
      INTEGER          LHARR( NRUNM )
      INTEGER          ISYMAR( NRUNM )
      INTEGER          MVALA( NRUNM )
      INTEGER          IOFARR( NRUNM )
      DOUBLE PRECISION CAARR( NRUNM )
      DOUBLE PRECISION CB1ARR( NRUNM )
      DOUBLE PRECISION CB2ARR( NRUNM )
      DOUBLE PRECISION CDARR( NRUNM )
      DOUBLE PRECISION CEARR( NRUNM )
      DOUBLE PRECISION CGARR( NRUNM )
      DOUBLE PRECISION CHARR1( NRUNM )
      DOUBLE PRECISION CHARR2( NRUNM )
      DOUBLE PRECISION CIARR( NRUNM )
C
      DOUBLE PRECISION DPARR( 2 ), XARR( NRMAX ), 
     1                 WVEC( ISVMAX ), SBRVEC( ISVMAX ),
     2                 SVFDC( NCFM, NRMAX, NDRVM+1, NDCS ),
     3                 PARAM( 10 )
C
      DOUBLE PRECISION COEFM1( NCFM, NCFM ), WORK1( NCFM ),
     1                 COEFM2( NCFM, NCFM ), WORK2( NCFM ),
     2                 A( NBNDMX, ISVMAX ),
     3                 GAUX( NTHMAX ), GAUW( NTHMAX )
C
      DOUBLE PRECISION FTF1( 2*NPHMAX ), FTF2( 2*NPHMAX ),
     1                 FTF3( 2*NPHMAX ), QST( LHLH2M, 3 ),
     2                 VF( NPHMAX, NTHMAX, 3 )
C
      DOUBLE PRECISION PA( NPMAX, NTHMAX ), DPA( NPMAX, NTHMAX ),
     1                 DR( NCVM ), DI( NCVM ), D3( NCVM ),
     2                 WORKEV( 3*NCVM ), WORKD( 3*ISVMAX ),
     3                 WORKL( 3*NCVM*NCVM + 6*NCVM )
C
      DOUBLE PRECISION RESID( ISVMAX ), V( ISVMAX, NCVM ),
     1                 W2( ISVMAX )
C
      LOGICAL SELECT( NCVM )
C
      CHARACTER *(200) LINE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER ITHEBC, IFORMF, LU, INSPCF, KL, IRUN, IFORM,
     1        N1, N2, NDRVS, NR, LH, NCV, NEV, MXIT, IEV
C
      INTEGER NH, ISYM, MLOW, MINC, MMAX, IVELBC, NTHP, NPHP,
     1        IWRITE, IAPP, IH, MXATT, IERR, NRUNS
C
      INTEGER I, ILEN, IOF, MVAL, ILNR, IRNR, IZF
C
      INTEGER  LULOG, LUNR, LUISP, LULH, LUSYM,
     1         LUM, LUCA, LUCB1, LUCB2, LUCD, LUCE,
     2         LUCG, LUCH, LUCI, LURES, LURI, LURO, LUBC
C
      CHARACTER *(80) FNLOG, FNAME, FNNR, FNISP, FNLH, FNSYM,
     1                FNM, FNCA, FNCB1, FNCB2, FNCD, FNCE, ROOT,
     2                FNCG, FNCH, FNCI, FNRES, FNRI, FNRO, FNBC
C
      DOUBLE PRECISION RI, RO, LOW, ZERO, CD, CI, CG, 
     1                 X1, X2, CA, CE, CH, CB1, CB2, ARTOL, DRSV,
     2                 GRR, GRI, CTOL, CH1, CH2
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IFORMF = 3, LU = 11, IFORM = 1,
     1            IWRITE = 3, IAPP = 4 )
C
      PARAMETER ( LOW = 1.0d-7, ZERO = 0.0d0, X1 = -1.0d0,
     1            X2 = 1.0d0 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO I = 1, 80
        ROOT(I:I) = ' '
        FNLOG(I:I) = ' '
        FNNR(I:I) = ' '
        FNISP(I:I) = ' '
        FNLH(I:I) = ' '
        FNSYM(I:I) = ' '
        FNM(I:I) = ' '
        FNCA(I:I) = ' '
        FNCB1(I:I) = ' '
        FNCB2(I:I) = ' '
        FNCD(I:I) = ' '
        FNCE(I:I) = ' '
        FNCG(I:I) = ' '
        FNCH(I:I) = ' '
        FNCI(I:I) = ' '
        FNRES(I:I) = ' '
        FNRI(I:I) = ' '
        FNRO(I:I) = ' '
        FNBC(I:I) = ' '
      ENDDO
C
 80   FORMAT(A)
C
C Start reading in input file
C First line contains root only
C
      PRINT *,' Program LINear ONSet 1.'
      PRINT *,' Input file name root.'
      PRINT *,' ======================================== '
 21   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 21
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 49
        ENDIF
      ENDDO
 49   CONTINUE
      ROOT = LINE(1:ILEN)
C
C next line should contain RI, RO, IVELBC, ITHEBC, LULOG
C RI must be atleast zero
C RO must be strictly greater than RI
C IVELBC = 1 --> no slip
C IVELBC = 2 --> stress free
C ITHEBC = 1 --> fixed tm inner and outer boundaries
C ITHEBC = 2 --> fixed tm inner hf outer
C ITHEBC = 3 --> fixed hf inner tm outer
C LULOG = 0 --> no output
C LULOG = 44 --> limited output
C LULOG = 45 --> extensive output
C
      PRINT *,' Enter RI, RO, IVELBC, ITHEBC, LULOG '
      PRINT *,' (ivelbc = 1 --> no slip )'
      PRINT *,' (ivelbc = 2 --> stress free)'
      PRINT *,' (ithebc = 1 --> fixed tm inner and outer)'
      PRINT *,' (ithebc = 2 --> fixed tm inner hf outer)'
      PRINT *,' (ithebc = 3 --> fixed hf inner tm outer)'
      PRINT *,' (lulog = 0  --> no output)'
      PRINT *,' (lulog = 44  --> limited output)'
      PRINT *,' (lulog = 45  --> extensive output)'
      PRINT *,' ======================================== '
 22   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 22
      READ ( LINE, * ) RI, RO, IVELBC, ITHEBC, LULOG
C
      IF ( RI.LT.ZERO ) THEN
        PRINT *,' RI = ', RI
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( RO.LE.RI ) THEN
        PRINT *,' RI = ', RI
        PRINT *,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C ri and ro seem ok
C
      DPARR( 1 ) = RI
      DPARR( 2 ) = RO
C
      IF ( IVELBC.NE.1 .AND. IVELBC.NE.2 ) THEN
        PRINT *,' IVELBC = ', IVELBC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ITHEBC.NE.1 .AND. ITHEBC.NE.2 .AND. ITHEBC.NE.3 ) THEN
        PRINT *,' ITHEBC = ', ITHEBC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LULOG.NE.0 .AND. LULOG.NE.144 .AND. LULOG.NE.145 .AND.
     1          LULOG.NE.44 .AND. LULOG.NE.45  ) THEN
        PRINT *,' LULOG = ', LULOG
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LULOG.EQ.144 ) LULOG = 44
      IF ( LULOG.EQ.145 ) LULOG = 45
C
C Fill in MHIBC, MHOBC and LARR arrays
C First poloidal velocity
C
      LARR( 1 ) = 0
      LARR( 2 ) = 0
      LARR( 3 ) = 0
C
      IF ( IVELBC.EQ.1 ) THEN
        MHIBC( 1 ) = 4
        MHOBC( 1 ) = 4
        MHIBC( 2 ) = 2
        MHOBC( 2 ) = 2
      ELSE
        MHIBC( 1 ) = 5
        MHOBC( 1 ) = 5
        MHIBC( 2 ) = 6
        MHOBC( 2 ) = 6
      ENDIF
C
      IF ( ITHEBC.EQ.1 ) THEN
        MHIBC( 3 ) = 2
        MHOBC( 3 ) = 2
      ENDIF
C
      IF ( ITHEBC.EQ.2 ) THEN
        MHIBC( 3 ) = 2
        MHOBC( 3 ) = 3
      ENDIF
C
      IF ( ITHEBC.EQ.3 ) THEN
        MHIBC( 3 ) = 3
        MHOBC( 3 ) = 2
      ENDIF
C
      FNLOG(1:ILEN) = ROOT(1:ILEN)
      FNLOG(ILEN+1:ILEN+4) = '.log'
      IF ( LULOG.NE.0 ) CALL FOPEN ( LULOG, FNLOG, IWRITE )
C
      LUNR  = 61
      LUISP = 62
      LULH  = 63
      LUSYM = 64
      LUM   = 65
      LUCA  = 66
      LUCB1 = 67
      LUCB2 = 68
      LUCD  = 69
      LUCE  = 70
      LUCG  = 71
      LUCH  = 72
      LUCI  = 73
      LURES = 74
      LURI  = 75
      LURO  = 76
      LUBC  = 77
C
      FNNR(1:ILEN) = ROOT(1:ILEN)
      FNISP(1:ILEN) = ROOT(1:ILEN)
      FNLH(1:ILEN) = ROOT(1:ILEN)
      FNSYM(1:ILEN) = ROOT(1:ILEN)
      FNM(1:ILEN) = ROOT(1:ILEN)
      FNCA(1:ILEN) = ROOT(1:ILEN)
      FNCB1(1:ILEN) = ROOT(1:ILEN)
      FNCB2(1:ILEN) = ROOT(1:ILEN)
      FNCD(1:ILEN) = ROOT(1:ILEN)
      FNCE(1:ILEN) = ROOT(1:ILEN)
      FNCG(1:ILEN) = ROOT(1:ILEN)
      FNCH(1:ILEN) = ROOT(1:ILEN)
      FNCI(1:ILEN) = ROOT(1:ILEN)
      FNRES(1:ILEN) = ROOT(1:ILEN)
      FNRI(1:ILEN) = ROOT(1:ILEN)
      FNRO(1:ILEN) = ROOT(1:ILEN)
      FNBC(1:ILEN) = ROOT(1:ILEN)
      FNAME(1:ILEN) = ROOT(1:ILEN)
C
      FNNR(ILEN+1:ILEN+3) = '.nr'
      FNISP(ILEN+1:ILEN+4) = '.isp'
      FNLH(ILEN+1:ILEN+3) = '.lh'
      FNSYM(ILEN+1:ILEN+5) = '.isym'
      FNM(ILEN+1:ILEN+2) = '.m'
      FNCA(ILEN+1:ILEN+3) = '.ca'
      FNCB1(ILEN+1:ILEN+4) = '.cb1'
      FNCB2(ILEN+1:ILEN+4) = '.cb2'
      FNCD(ILEN+1:ILEN+3) = '.cd'
      FNCE(ILEN+1:ILEN+3) = '.ce'
      FNCG(ILEN+1:ILEN+3) = '.cg'
      FNCH(ILEN+1:ILEN+3) = '.ch'
      FNCI(ILEN+1:ILEN+3) = '.ci'
      FNRES(ILEN+1:ILEN+4) = '.res'
      FNRI(ILEN+1:ILEN+3) = '.ri'
      FNRO(ILEN+1:ILEN+3) = '.ro'
      FNBC(ILEN+1:ILEN+3) = '.bc'
C
      CALL FOPEN( LUNR   , FNNR   , IAPP )
      CALL FOPEN( LUISP  , FNISP  , IAPP )
      CALL FOPEN( LULH   , FNLH   , IAPP )
      CALL FOPEN( LUSYM  , FNSYM  , IAPP )
      CALL FOPEN( LUM    , FNM    , IAPP )
      CALL FOPEN( LUCA   , FNCA   , IAPP )
      CALL FOPEN( LUCB1  , FNCB1  , IAPP )
      CALL FOPEN( LUCB2  , FNCB2  , IAPP )
      CALL FOPEN( LUCD   , FNCD   , IAPP )
      CALL FOPEN( LUCE   , FNCE   , IAPP )
      CALL FOPEN( LUCG   , FNCG   , IAPP )
      CALL FOPEN( LUCH   , FNCH   , IAPP )
      CALL FOPEN( LUCI   , FNCI   , IAPP )
      CALL FOPEN( LURES  , FNRES  , IAPP )
      CALL FOPEN( LURI   , FNRI   , IAPP )
      CALL FOPEN( LURO   , FNRO   , IAPP )
      CALL FOPEN( LUBC   , FNBC   , IAPP )
C
C Next line should contain DRSV, NEV, NCV, MXATT, CTOL
C
      PRINT *,' Enter DRSV, NEV, NCV, MXATT, CTOL '
      PRINT *,' DRSV is real shift value for solving e.system.'
      PRINT *,' NEV = number of requested eigenvalues.'
      PRINT *,' NCV = length of Arnoldi factorisation.'
      PRINT *,' NCV must be at least NEV + 2.'
      PRINT *,' NCV must be less than ', NCVM
      PRINT *,' MXATT = number of attempts allowed '
      PRINT *,' to iterate to Critical Rayleigh number.'
      PRINT *,' CTOL = convergence criterion.'
 23   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 23
      READ ( LINE, * ) DRSV, NEV, NCV, MXATT, CTOL
C
      IF ( NEV.LT.1 .OR. NEV.GT.(NCVM-2) ) THEN
        PRINT *,' NEV  = ', NEV
        PRINT *,' NCVM = ', NCVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NCV.LT.(NEV+2) .OR. NCV.GT.NCVM ) THEN
        PRINT *,' NEV  = ', NEV
        PRINT *,' NCV  = ', NCV
        PRINT *,' NCVM = ', NCVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now loop around until we have solved all the
C required problems ... these parameters are
C
C NR, INSPCF, LH, ISYM, MVAL, IOF
C CA, CB1, CB2, CD, CE, CG, CH, CI
C
C (inspcf = 1 --> equally spaced nodes )
C (inspcf = 2 --> Chebyshev nodes )
C isym = 1 for Equatorially symmetric modes
C isym = 2 for Equatorially anti-symmetric modes
C IOF = 0 --> no write out of eigenfunctions etc.
C IOF = 1 --> standard output of eigenfunctions etc.
C IOF = 2 --> visual output of radial functions.
C
      NRUNS = 0
 24   CONTINUE
      READ ( 5, 80, END=300 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 24
      NRUNS = NRUNS + 1
      IF ( NRUNS.GT.NRUNM ) THEN
        PRINT *,' NRUNM exceeded. Either recompile'
        PRINT *,' or change inputs.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      READ ( LINE, * )  NRARR( NRUNS ), ISPARR( NRUNS ),
     1       LHARR( NRUNS ), ISYMAR( NRUNS ), MVALA( NRUNS ), 
     2       IOFARR( NRUNS ), CAARR( NRUNS ), CB1ARR( NRUNS ),
     3       CB2ARR( NRUNS ), CDARR( NRUNS ), CEARR( NRUNS ),
     4       CGARR( NRUNS ), CHARR1( NRUNS ), CHARR2( NRUNS ),
     5       CIARR( NRUNS )
      GOTO 24
 300  CONTINUE
C
C
C
C o.k. - now we have read in all our parameters
C so let's loop around one by one ....
C
      DO IRUN = 1, NRUNS
C
       NR = NRARR( IRUN )
       INSPCF = ISPARR( IRUN )
       LH = LHARR( IRUN )
       ISYM = ISYMAR( IRUN )
       MVAL = MVALA( IRUN )
       IOF  = IOFARR( IRUN )
       CA   = CAARR( IRUN )
       CB1  = CB1ARR( IRUN )
       CB2  = CB2ARR( IRUN )
       CD   = CDARR( IRUN )
       CE   = CEARR( IRUN )
       CG   = CGARR( IRUN )
       CH1  = CHARR1( IRUN )
       CH2  = CHARR2( IRUN )
       CI   = CIARR( IRUN )
C
C Check values for NR, INSPCF
C
      IF ( NR.LT.10 .OR. NR.GT.NRMAX ) THEN
        PRINT *,' NR = ', NR,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( INSPCF.NE.1 .AND. INSPCF.NE.2 ) THEN
        PRINT *,' INSPCF = ', INSPCF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( INSPCF.EQ.1 ) THEN
        CALL ESNAAS( NR, XARR, RI, RO )
      ELSE
        CALL ZCPAAS( NR, XARR, RI, RO )
      ENDIF
C
      MLOW = MVAL
      MINC = 1
      MMAX = MVAL
C
      IF ( LH.GT.LHMAX ) THEN
        PRINT *,' LH    = ', LH
        PRINT *,' LHMAX = ', LHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ISYM.NE.1 .AND. ISYM.NE.2 .AND. ISYM.NE.3 ) THEN
        PRINT *,' ISYM  = ', ISYM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now read to calculate Legdendre Functions etc.
C
      CALL ONTPPF( LH, MMAX, NTHP, NPHP, NTHMAX, NPHMAX )
      CALL GAUWTS ( X1, X2, GAUX, GAUW, NTHP )
      CALL SCHNLA ( PA, DPA, GAUX, LH, NTHP )
C
C Now calculate the full harmonic sets
C
      CALL VTHMSR( ISYM, NH, MLOW, MINC, MMAX, NHMAX, MHT, MHL,
     1             MHM, LH, LHMAX )
C
      CALL CINDSW ( NH, MHT, MHTR )
      PRINT *,' Total of ', NH,' harmonics. '
C
      INARR( 1 ) = IFORMF
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C
      KL = (NBN+1)*NH - 1
      N1 = 3*KL + 1
      N2 = NH*NR
C
      IF ( CD.LT.LOW .OR. CI.LT.LOW .OR. CG.LT.ZERO ) THEN
        PRINT *,' CD = ', CD
        PRINT *,' CI = ', CI
        PRINT *,' CG = ', CG
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( CA.LT.ZERO .OR. CE.LT.ZERO ) THEN
        PRINT *,' CA = ', CA
        PRINT *,' CE = ', CE
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now set up array MHP
C this is straightforward as we only have 3 types
C
      DO IH = 1, NH
        MHP( IH ) = MHT( IH )
      ENDDO
C
c     FNAME = 'Harmonic_file'
c     CALL HMFWT( NH, MHT, MHL, MHM, MHP, NDCS, MHIBC, MHOBC,
c    1                  LU, FNAME )
C
C Form SVFDC matrix
C
      NDRVS = 1
      CALL SVFDCF( NR, NDCS, NBN, 1, 1, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, WORK1, WORK2 )
C
      NDRVS = 4
      CALL SVFDCF( NR, NDCS, NBN, 2, NR-1, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, WORK1, WORK2 )
C
      NDRVS = 1
      CALL SVFDCF( NR, NDCS, NBN, NR, NR, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, WORK1, WORK2 )
C
 50   CONTINUE
C
C*****************************
C*****************************
C****   Form matrix ....  ****
C*****************************
C*****************************
C
      ARTOL = 0.000001d0
      MXIT  = 400
C
      PARAM(  1 ) = CA
      PARAM(  2 ) = CB1
      PARAM(  3 ) = CB2
      PARAM(  4 ) = CD
      PARAM(  5 ) = CE
      PARAM(  6 ) = CG
      PARAM(  7 ) = CH1
      PARAM(  8 ) = CH2
      PARAM(  9 ) = CI
      PARAM( 10 ) = CTOL
C
      CALL  LOCRBR( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1   NCFM, NDRVM, N1, N2, NDCS, NTHP, NPHP, MMAX, LH, IEV, MHIBC,
     2   MHOBC, NEV, NCV, NCVM, MXIT, IPIV, LULOG, SVFDC, A, XARR,
     3   GAUX, GAUW, PA, DPA, SBRVEC, RESID, W2, WVEC, FTF1, FTF2,
     4   FTF3, VF, QST, SELECT, ARTOL, DRSV, PARAM, MXATT, IERR,
     5   DR, DI, D3, WORKEV, WORKD, WORKL, V )
C
      CA  = PARAM(  1 )
      CB1 = PARAM(  2 )
      CB2 = PARAM(  3 )
      CD  = PARAM(  4 )
      CE  = PARAM(  5 )
      CG  = PARAM(  6 )
      CH  = PARAM(  7 )
      CI  = PARAM(  8 )
      GRR = PARAM(  9 )
      GRI = PARAM( 10 )
C
 701  FORMAT(I4)
 702  FORMAT(I1)
 703  FORMAT(1PD16.8)
 704  FORMAT(I2,I2)
 705  FORMAT(1PD16.8,1PD16.8,1PD16.8)
C
      IF ( IERR.GT.0 ) THEN
        WRITE ( LUNR, 701 ) NR
        WRITE ( LUISP, 702 ) INSPCF
        WRITE ( LULH, 701 ) LH
        WRITE ( LUSYM, 702 ) ISYM
        WRITE ( LUM, 701 ) MVAL
        WRITE ( LUCA, 703 ) CA
        WRITE ( LUCB1 , 703 ) CB1
        WRITE ( LUCB2 , 703 ) CB2
        WRITE ( LUCD  , 703 ) CD
        WRITE ( LUCE  , 703 ) CE
        WRITE ( LUCG  , 703 ) CG
        WRITE ( LUCH  , 703 ) CH
        WRITE ( LUCI  , 703 ) CI
        WRITE ( LURES , 705 ) CH, GRR, GRI
        WRITE ( LULOG , 722 ) NR, LH, MVAL, CG, CH, GRI
        WRITE ( LURI  , 703 ) RI
        WRITE ( LURO  , 703 ) RO
        WRITE ( LUBC  , 704 ) IVELBC, ITHEBC
      ELSE
        PRINT *,' IRUN = ', IRUN,' IERR = ', IERR
        PRINT *,' Solution unsuccessful.'
      ENDIF
C
 722    FORMAT('Result: N:',I4,' L:',I4,' M:',
     1        I4,' G: ',f16.6,' R: ',f16.6,' o: ',f16.6)
C
      CALL FLUSH ( LULOG )
      CALL FLUSH ( LURES )
C
C We have to write out the critical eigenvector
C in the standard format ...
C
      IF ( IOF.EQ.1 ) THEN
C       .
C       . We need to apply a filename
C       . for the .vecs, .xarr and .ints files
C       . FNAME(1:ILEN) contains root
C       .
        IF ( IRUN.GE.1 .AND. IRUN.LT.10 )
     1     WRITE (  FNAME(ILEN+1:ILEN+7), 811 ) IRUN
        IF ( IRUN.GE.10 .AND. IRUN.LT.100 )
     1     WRITE (  FNAME(ILEN+1:ILEN+7), 812 ) IRUN
        IF ( IRUN.GE.100 .AND. IRUN.LT.1000 )
     1     WRITE (  FNAME(ILEN+1:ILEN+7), 813 ) IRUN
C
C Write out the harmonic integers file
C
        FNAME(ILEN+8:ILEN+12) = '.ints'
        FNAME = FNAME(1:ILEN+12)
        CALL HMFWT( NH, MHT, MHL, MHM, MHP, NDCS, MHIBC, MHOBC,
     1            LU, FNAME )
C
C Write out the eigenvectors
C
        FNAME(ILEN+8:ILEN+12) = '.vecr'
        FNAME = FNAME(1:ILEN+12)
        CALL SVFWT( INARR, LU, IFORM, SBRVEC, FNAME )
        FNAME(ILEN+8:ILEN+12) = '.veci'
        FNAME = FNAME(1:ILEN+12)
        CALL SVFWT( INARR, LU, IFORM, W2, FNAME )
C
C Write out radial node data
C
        FNAME(ILEN+8:ILEN+12) = '.xarr'
        FNAME = FNAME(1:ILEN+12)
        CALL XARRWT( NR, XARR, LU, FNAME, IFORM )
C       .
      ENDIF
C
C We have to write out the critical eigenvector
C in the visual display format
C
      IF ( IOF.EQ.2 ) THEN
C       .
C       . We need to apply a filename
C       . for the .disp file
C       . FNAME(1:ILEN) contains root
C       .
        IF ( IRUN.GE.1 .AND. IRUN.LT.10 )
     1     WRITE (  FNAME(ILEN+1:ILEN+7), 811 ) IRUN
        IF ( IRUN.GE.10 .AND. IRUN.LT.100 )
     1     WRITE (  FNAME(ILEN+1:ILEN+7), 812 ) IRUN
        IF ( IRUN.GE.100 .AND. IRUN.LT.1000 )
     1     WRITE (  FNAME(ILEN+1:ILEN+7), 813 ) IRUN
C
C Write out the harmonic integers file
C
        FNAME(ILEN+8:ILEN+12) = '.disp'
        FNAME = FNAME(1:ILEN+12)
        ILNR  = 1
        IRNR  = NR
        IZF   = 2
        CALL FOPEN( LU, FNAME, IWRITE )
        CALL SVPRNT ( SBRVEC, NR, INARR, DPARR, XARR, MHT, MHL,
     1                MHM, ILNR, IRNR, LU, IZF, IFORM )
        CALL FCLOSE ( LU, FNAME, 'Error in closing FNAME' )
C       .
      ENDIF
C
 811  FORMAT('.run00',I1)
 812  FORMAT('.run0',I2)
 813  FORMAT('.run',I3)
C
      ENDDO
      PRINT *,' Program finishing.'
      CALL FCLOSE( LUNR   , FNNR   , 'Error' )
      CALL FCLOSE( LUISP  , FNISP  , 'Error' )
      CALL FCLOSE( LULH   , FNLH   , 'Error' )
      CALL FCLOSE( LUSYM  , FNSYM  , 'Error' )
      CALL FCLOSE( LUM    , FNM    , 'Error' )
      CALL FCLOSE( LUCA   , FNCA   , 'Error' )
      CALL FCLOSE( LUCB1  , FNCB1  , 'Error' )
      CALL FCLOSE( LUCB2  , FNCB2  , 'Error' )
      CALL FCLOSE( LUCD   , FNCD   , 'Error' )
      CALL FCLOSE( LUCE   , FNCE   , 'Error' )
      CALL FCLOSE( LUCG   , FNCG   , 'Error' )
      CALL FCLOSE( LUCH   , FNCH   , 'Error' )
      CALL FCLOSE( LUCI   , FNCI   , 'Error' )
      CALL FCLOSE( LURES  , FNRES  , 'Error' )
      CALL FCLOSE( LURI   , FNRI   , 'Error' )
      CALL FCLOSE( LURO   , FNRO   , 'Error' )
      CALL FCLOSE( LUBC   , FNBC   , 'Error' )
      IF ( LULOG.NE.0 )
     1     CALL FCLOSE ( LULOG, FNLOG, 'Error in closing file.' )
      STOP
 41   FORMAT('Harmonic ',I4,' Type= ',I1,' L= ',I3,' M= ',I3 )
      END
C*********************************************************************

