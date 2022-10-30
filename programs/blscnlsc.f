C*********************************************************************
C                                                                    C
C Steve Gibbons - Boundary-Locked Steady, Non-Linear Solutions Calc. C
C Thu Mar 16 14:40:21 GMT 2000                                       C
C                                                                    C
C                                                                    C
C                                                                    C
C                                                                    C
C*********************************************************************
      PROGRAM blscnlsc
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, LHMAX, NHMAX, NPHMAX, NTHMAX, KLMAX,
     1        NBNDMX, NDCS, NDRVM, ISVMAX, NPMAX, LHLH2M, NCFM,
     2        NBN, NRUNM, NITHMX
C
      PARAMETER ( NRMAX = 40, LHMAX = 62, NHMAX = 162, NBN = 3,
     1            NTHMAX = 64, NPHMAX = 128, KLMAX = (NBN+1)*NHMAX-1,
     2            NBNDMX = 3*KLMAX+1, NDCS = 4, NDRVM = 4, 
     3            ISVMAX = NRMAX*NHMAX, NITHMX = 100 )
      PARAMETER ( NPMAX = (LHMAX+1)*(LHMAX+2)/2, 
     1            LHLH2M = LHMAX*(LHMAX+2), NCFM = 2*NBN + 1,
     2            NRUNM = 100 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER MHIBC( 4 ), MHOBC( 4 ), LARR( 4 ), INARR( 3 ),
     1        IPIV( ISVMAX ), MHT( NHMAX ), MHL( NHMAX ),
     2        MHM( NHMAX ), MHP( NHMAX ), MHTR( NHMAX ),
     3        IWORK( NCFM ), MHI( NHMAX ), MHP2( NHMAX )
C
      INTEGER          NRARR( NRUNM )
      INTEGER          ISPARR( NRUNM )
      INTEGER          LHARR( NRUNM )
      INTEGER          ISYMAR( NRUNM )
      INTEGER          MLOWA( NRUNM )
      INTEGER          MINCA( NRUNM )
      INTEGER          MMAXA( NRUNM )
      INTEGER          IOFARR( NRUNM )
      DOUBLE PRECISION CCARR( NRUNM )
      DOUBLE PRECISION CB1ARR( NRUNM )
      DOUBLE PRECISION CB2ARR( NRUNM )
      DOUBLE PRECISION CDARR( NRUNM )
      DOUBLE PRECISION CFARR( NRUNM )
      DOUBLE PRECISION CGARR( NRUNM )
      DOUBLE PRECISION CHARR( NRUNM )
      DOUBLE PRECISION CIARR( NRUNM )
      DOUBLE PRECISION EPSIAR( NRUNM )
      DOUBLE PRECISION EPSOAR( NRUNM )
C
      DOUBLE PRECISION DPARR( 2 ), XARR( NRMAX ), VEC( ISVMAX ),
     1                 SVFDC( NCFM, NRMAX, NDRVM+1, NDCS ),
     2                 PARAM( 12 ), FDCM( NCFM, NRMAX, 1 )
C
      DOUBLE PRECISION COEFM1( NCFM, NCFM ), WORK1( NCFM ),
     1                 COEFM2( NCFM, NCFM ), WORK2( NCFM ),
     2                 A( NBNDMX, ISVMAX ), VEC2( ISVMAX ),
     3                 GAUX( NTHMAX ), GAUW( NTHMAX )
C
      DOUBLE PRECISION PA( NPMAX, NTHMAX ), DPA( NPMAX, NTHMAX ),
     1                 CAFIT( 3, NITHMX ), SHCI( LHLH2M ), SHCIM,
     2                 SHCO( LHLH2M ), SHCOM, SHCI2( LHLH2M ),
     3                 SHCO2( LHLH2M )
C
      CHARACTER *(200) LINE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER ITHEBC, IFORMF, LU, INSPCF, KL, IRUN, IFORM,
     1        N1, N2, NDRVS, NR, LH, IOP, INDFUN
C
      INTEGER NH, ISYM, MLOW, MINC, MMAX, IVELBC, NTHP, NPHP,
     1        IWRITE, IAPP, IH, MXATT, IERR, NRUNS, IR, IITH
C
      INTEGER I, ILEN, IOF, ILNR, IRNR, IZF, IREAD,
     1        L, M, ICS, IND, INDSHC, KIB, KOB, NITH, IHD
C
      INTEGER  LULOG, LUNR, LUISP, LULH, LUSYM,
     1         LUM, LUCC, LUCB1, LUCB2, LUCD, LUCF,
     2         LUCG, LUCH, LUCI, LURES, LURI, LURO, LUBC,
     3         LUEPSI, LUEPSO
C
      CHARACTER *(2)  BCH
      CHARACTER *(80) FNLOG, FNAME, FNNR, FNISP, FNLH, FNSYM,
     1                FNM, FNCC, FNCB1, FNCB2, FNCD, FNCF, ROOT,
     2                FNCG, FNCH, FNCI, FNRES, FNRI, FNRO, FNBC,
     3                FNEPSI, FNEPSO
C
      DOUBLE PRECISION RI, RO, LOW, ZERO, CD, CI, CG, CAK, CBK, CCK,
     1                 X1, X2, CC, CF, CH, CB1, CB2, RAD,
     2                 CTOL, EPSI, EPSO, COEF, DERV( 1 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IFORMF = 3, LU = 12, IFORM = 1, IOP = 0,
     1            IWRITE = 3, IAPP = 4, IREAD = 1 )
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
        FNCC(I:I) = ' '
        FNCB1(I:I) = ' '
        FNCB2(I:I) = ' '
        FNCD(I:I) = ' '
        FNCF(I:I) = ' '
        FNCG(I:I) = ' '
        FNCH(I:I) = ' '
        FNCI(I:I) = ' '
        FNRES(I:I) = ' '
        FNRI(I:I) = ' '
        FNRO(I:I) = ' '
        FNBC(I:I) = ' '
        FNEPSI(I:I) = ' '
        FNEPSO(I:I) = ' '
      ENDDO
C
 80   FORMAT(A)
C
C Start reading in input file
C First line contains root only
C
      PRINT *,' Program Boundary Locked Steady Solution Solve.'
      PRINT *,' Enter name which contains boundary coef.s '
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
      FNAME = LINE(1:ILEN)
C
C Zero arrays of boundary coefficients
C
      SHCIM = ZERO
      SHCOM = ZERO
      CALL VECOP( SHCI, ZERO, LHLH2M, IOP )
      CALL VECOP( SHCO, ZERO, LHLH2M, IOP )
C
C Now let's read in boundary coefficients
C
      CALL FOPEN( LU, FNAME, IREAD )
 92   CONTINUE
      READ ( LU, 80, END = 93 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 92
C
      BCH = LINE(1:2)
      READ ( LINE(3:80), * ) L, M, ICS, COEF
      IF ( BCH.NE.'IB' .AND. BCH.NE.'OB' ) THEN
        PRINT *,' Start of boundary coefficient line'
        PRINT *,' must either be IB (inner) or '
        PRINT *,' OB (outer) boundary.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( L.GT.LHMAX ) THEN
        PRINT *,' Boundary coefficient with L = ', L
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IND = INDSHC( L, M, ICS )
      IF ( BCH.EQ.'IB' ) THEN
        IF ( IND.EQ.0 ) THEN
          SHCIM = COEF
        ELSE
          SHCI( IND ) = COEF
        ENDIF
      ELSE
        IF ( IND.EQ.0 ) THEN
          SHCOM = COEF
        ELSE
          SHCO( IND ) = COEF
        ENDIF
      ENDIF
C
C zeros mean heat flux variation ...
C
      SHCIM = ZERO
      SHCOM = ZERO
C
      GOTO 92
 93   CONTINUE
      CALL FCLOSE( LU, FNAME, 'Error' )
C
C next line contains name of file which contains
C spherical harmonic coefficients
C
      PRINT *,' Input file name root.'
      PRINT *,' ======================================== '
 121  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 121
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 149
        ENDIF
      ENDDO
 149  CONTINUE
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
      PRINT *,' (lulog = 6  --> output to screen)'
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
      IF ( LULOG.NE.0 .AND. LULOG.NE.6 .AND. LULOG.NE.145 .AND.
     1         LULOG.NE.45 ) THEN
        PRINT *,' LULOG = ', LULOG
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Fill in MHIBC, MHOBC and LARR arrays
C First poloidal velocity
C
      LARR( 1 ) = 0
      LARR( 2 ) = 0
      LARR( 3 ) = 0
      LARR( 4 ) = 0
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
        KIB        = 1
        KOB        = 1
      ENDIF
C
      IF ( ITHEBC.EQ.2 ) THEN
        MHIBC( 3 ) = 2
        MHOBC( 3 ) = 3
        KIB        = 1
        KOB        = 2
      ENDIF
C
      IF ( ITHEBC.EQ.3 ) THEN
        MHIBC( 3 ) = 3
        MHOBC( 3 ) = 2
        KIB        = 2
        KOB        = 1
      ENDIF
C
      MHIBC( 4 ) = 1
      MHOBC( 4 ) = 1
C
      FNLOG(1:ILEN) = ROOT(1:ILEN)
      FNLOG(ILEN+1:ILEN+4) = '.log'
      IF ( LULOG.NE.0 .AND. LULOG.NE.6 ) 
     1              CALL FOPEN ( LULOG, FNLOG, IWRITE )
C
      LUNR   = 61
      LUISP  = 62
      LULH   = 63
      LUSYM  = 64
      LUM    = 65
      LUCC   = 66
      LUCB1  = 67
      LUCB2  = 68
      LUCD   = 69
      LUCF   = 70
      LUCG   = 71
      LUCH   = 72
      LUCI   = 73
      LURES  = 74
      LURI   = 75
      LURO   = 76
      LUBC   = 77
      LUEPSI = 78
      LUEPSO = 79
C
      FNNR(1:ILEN) = ROOT(1:ILEN)
      FNISP(1:ILEN) = ROOT(1:ILEN)
      FNLH(1:ILEN) = ROOT(1:ILEN)
      FNSYM(1:ILEN) = ROOT(1:ILEN)
      FNM(1:ILEN) = ROOT(1:ILEN)
      FNCC(1:ILEN) = ROOT(1:ILEN)
      FNCB1(1:ILEN) = ROOT(1:ILEN)
      FNCB2(1:ILEN) = ROOT(1:ILEN)
      FNCD(1:ILEN) = ROOT(1:ILEN)
      FNCF(1:ILEN) = ROOT(1:ILEN)
      FNCG(1:ILEN) = ROOT(1:ILEN)
      FNCH(1:ILEN) = ROOT(1:ILEN)
      FNCI(1:ILEN) = ROOT(1:ILEN)
      FNRES(1:ILEN) = ROOT(1:ILEN)
      FNRI(1:ILEN) = ROOT(1:ILEN)
      FNRO(1:ILEN) = ROOT(1:ILEN)
      FNBC(1:ILEN) = ROOT(1:ILEN)
      FNEPSI(1:ILEN) = ROOT(1:ILEN)
      FNEPSO(1:ILEN) = ROOT(1:ILEN)
      FNAME(1:ILEN) = ROOT(1:ILEN)
C
      FNNR(ILEN+1:ILEN+3) = '.nr'
      FNISP(ILEN+1:ILEN+4) = '.isp'
      FNLH(ILEN+1:ILEN+3) = '.lh'
      FNSYM(ILEN+1:ILEN+5) = '.isym'
      FNM(ILEN+1:ILEN+2) = '.m'
      FNCC(ILEN+1:ILEN+3) = '.cc'
      FNCB1(ILEN+1:ILEN+4) = '.cb1'
      FNCB2(ILEN+1:ILEN+4) = '.cb2'
      FNCD(ILEN+1:ILEN+3) = '.cd'
      FNCF(ILEN+1:ILEN+3) = '.cf'
      FNCG(ILEN+1:ILEN+3) = '.cg'
      FNCH(ILEN+1:ILEN+3) = '.ch'
      FNCI(ILEN+1:ILEN+3) = '.ci'
      FNRES(ILEN+1:ILEN+4) = '.res'
      FNRI(ILEN+1:ILEN+3) = '.ri'
      FNRO(ILEN+1:ILEN+3) = '.ro'
      FNBC(ILEN+1:ILEN+3) = '.bc'
      FNEPSI(ILEN+1:ILEN+5) = '.epsi'
      FNEPSO(ILEN+1:ILEN+5) = '.epso'
C
      CALL FOPEN( LUNR   , FNNR   , IAPP )
      CALL FOPEN( LUISP  , FNISP  , IAPP )
      CALL FOPEN( LULH   , FNLH   , IAPP )
      CALL FOPEN( LUSYM  , FNSYM  , IAPP )
      CALL FOPEN( LUM    , FNM    , IAPP )
      CALL FOPEN( LUCC   , FNCC   , IAPP )
      CALL FOPEN( LUCB1  , FNCB1  , IAPP )
      CALL FOPEN( LUCB2  , FNCB2  , IAPP )
      CALL FOPEN( LUCD   , FNCD   , IAPP )
      CALL FOPEN( LUCF   , FNCF   , IAPP )
      CALL FOPEN( LUCG   , FNCG   , IAPP )
      CALL FOPEN( LUCH   , FNCH   , IAPP )
      CALL FOPEN( LUCI   , FNCI   , IAPP )
      CALL FOPEN( LURES  , FNRES  , IAPP )
      CALL FOPEN( LURI   , FNRI   , IAPP )
      CALL FOPEN( LURO   , FNRO   , IAPP )
      CALL FOPEN( LUBC   , FNBC   , IAPP )
      CALL FOPEN( LUEPSI , FNEPSI , IAPP )
      CALL FOPEN( LUEPSO , FNEPSO , IAPP )
C
C Next line should contain MXATT, CTOL
C
      PRINT *,' Enter MXATT, CTOL '
      PRINT *,' MXATT = number of attempts allowed '
      PRINT *,' to iterate to non-linear solution.'
      PRINT *,' CTOL = convergence criterion.'
 23   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 23
      READ ( LINE, * ) MXATT, CTOL
C
C Now loop around until we have solved all the
C required problems ... these parameters are
C
C NR, INSPCF, LH, ISYM, MLOW, MINC, MMAX, IOF
C CC, CB1, CB2, CD, CF, CG, CH, CI
C
C (inspcf = 1 --> equally spaced nodes )
C (inspcf = 2 --> Chebyshev nodes )
C isym = 1 for Equatorially symmetric modes
C isym = 2 for Equatorially anti-symmetric modes
C isym = 3 for both symmetries
C IOF = 0 --> no write out of eigenfunctions etc.
C IOF = 1 --> standard output of homogeneous eigenfunctions
C IOF = 2 --> visual output of radial functions.
C IOF = 3 --> standard output of full solution (with inhom. bdry)
C IOF = 4 --> standard output of hom and inhom. boundary)
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
     1       LHARR( NRUNS ), ISYMAR( NRUNS ), MLOWA( NRUNS ), 
     2       MINCA( NRUNS ), MMAXA( NRUNS ), IOFARR( NRUNS ),
     3       CCARR( NRUNS ), CB1ARR( NRUNS ), CB2ARR( NRUNS )
 25   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 25
      READ ( LINE, * ) CDARR( NRUNS ), CFARR( NRUNS ),
     1       CGARR( NRUNS ), CHARR( NRUNS ), CIARR( NRUNS ),
     2       EPSIAR( NRUNS ), EPSOAR( NRUNS )
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
C Zero VEC for new run ...
C
      CALL VECOP( VEC, ZERO, ISVMAX, IOP )
C
       NR = NRARR( IRUN )
       INSPCF = ISPARR( IRUN )
       LH = LHARR( IRUN )
       ISYM = ISYMAR( IRUN )
       MLOW = MLOWA( IRUN )
       MINC = MINCA( IRUN )
       MMAX = MMAXA( IRUN )
       IOF  = IOFARR( IRUN )
       CC   = CCARR( IRUN )
       CB1  = CB1ARR( IRUN )
       CB2  = CB2ARR( IRUN )
       CD   = CDARR( IRUN )
       CF   = CFARR( IRUN )
       CG   = CGARR( IRUN )
       CH   = CHARR( IRUN )
       CI   = CIARR( IRUN )
       EPSI = EPSIAR( IRUN )
       EPSO = EPSOAR( IRUN )
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
C Now set up array MHP
C this is straightforward as we only have 3 types
C
      DO IH = 1, NH
        MHP( IH ) = MHT( IH )
      ENDDO
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
      CALL FDCMBD( NR, NBN, 2, NR-1, 2, NR-1, NCFM,
     1             NCFM, 1, 1, IWORK, XARR, FDCM,
     2             COEFM1, WORK1, WORK2 )
C
      NDRVS = 1
      CALL SVFDCF( NR, NDCS, NBN, NR, NR, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, WORK1, WORK2 )
C
C*****************************
C Here we must form the array
C CAFIT which contains the
C coefficients which represent
C the inhomogeneous boundary
C condition.
C*****************************
      CALL SHCANC( LH, SHCI, SHCI2, EPSI )
      CALL SHCANC( LH, SHCO, SHCO2, EPSO )
      NITH = 0
      CALL ITHCAR( KIB, KOB, NITH, NITHMX, NH, MHT, MHL, MHM,
     1             MHI, LH, RI, RO, SHCI2, SHCIM, SHCO2, SHCOM,
     2             CAFIT )
C*****************************
      PARAM(  2 ) = CB1
      PARAM(  3 ) = CB2
      PARAM(  4 ) = CC 
      PARAM(  5 ) = CD 
      PARAM(  7 ) = CF 
      PARAM(  8 ) = CG 
      PARAM(  9 ) = CH 
      PARAM( 10 ) = CI    
      PARAM( 12 ) = CTOL
C
      CALL BLCNRS( VEC, INARR, NR, PARAM, MXATT, IERR, A, N1,
     1  N2, KL, MHT, MHL, MHM, MHP, MHTR, NBN, NCFM, NDRVM, NDCS,
     2  NTHP, NPHP, MMAX, LH, MHIBC, MHOBC, IPIV, LULOG, SVFDC, XARR,
     3  GAUX, GAUW, PA, DPA, FDCM, MHI, CAFIT )
C
 701  FORMAT(I4)
 702  FORMAT(I1)
 703  FORMAT(1PD16.8)
 704  FORMAT(I2,I2)
 705  FORMAT(1PD16.8)
 706  FORMAT(I4,I4,I4)
C
      IF ( IERR.GT.0 ) THEN
        WRITE ( LUNR, 701 ) NR
        WRITE ( LUISP, 702 ) INSPCF
        WRITE ( LULH, 701 ) LH
        WRITE ( LUSYM, 702 ) ISYM
        WRITE ( LUM, 706 ) MLOW, MINC, MMAX
        WRITE ( LUCC, 703 ) CC
        WRITE ( LUCB1 , 703 ) CB1
        WRITE ( LUCB2 , 703 ) CB2
        WRITE ( LUCD  , 703 ) CD
        WRITE ( LUCF  , 703 ) CF
        WRITE ( LUCG  , 703 ) CG
        WRITE ( LUCH  , 703 ) CH
        WRITE ( LUCI  , 703 ) CI
        WRITE ( LURES , 705 ) CH
        WRITE ( LURI  , 703 ) RI
        WRITE ( LURO  , 703 ) RO
        WRITE ( LUBC  , 704 ) IVELBC, ITHEBC
      ELSE
        PRINT *,' IRUN = ', IRUN,' IERR = ', IERR
        PRINT *,' Solution unsuccessful.'
      ENDIF
C
C We have to write out the critical eigenvector
C in the standard format ...
C
      IF ( IOF.EQ.1 .OR. IOF.EQ.4 ) THEN
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
        FNAME(ILEN+8:ILEN+12) = '.vecs'
        FNAME = FNAME(1:ILEN+12)
        CALL SVFWT( INARR, LU, IFORM, VEC, FNAME )
C
C Write out radial node data
C
        FNAME(ILEN+8:ILEN+12) = '.xarr'
        FNAME = FNAME(1:ILEN+12)
        CALL XARRWT( NR, XARR, LU, FNAME, IFORM )
C       .
      ENDIF
C
C We have to write out the eigenvector
C WITH THE INHOMOGENEOUS TEMPERATURE
C
      IF ( IOF.EQ.3 .OR. IOF.EQ.4 ) THEN
        IHD = 0
C       .
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.3 ) THEN
            MHP2( IH ) = 4
            IITH   = MHI( IH )
            CAK    = CAFIT( 1, IITH )
            CBK    = CAFIT( 2, IITH )
            CCK    = CAFIT( 3, IITH )
            DO IR = 1, NR
              RAD = XARR( IR )
C             .
              IND = INDFUN( IR, IH, INARR )
              DERV( 1 ) = 0.0d0
C
              CALL ITFA( RAD, RI, RO, CAK, CBK, CCK, DERV, IHD )
C             .
              VEC2( IND ) = VEC( IND ) + DERV( 1 )
C             .
            ENDDO
         ELSE
            MHP2( IH ) = MHP( IH )
            DO IR = 1, NR
              IND = INDFUN( IR, IH, INARR )
              VEC2( IND ) = VEC( IND )
            ENDDO
         ENDIF
        ENDDO
C       .
C       . We need to apply a filename
C       . for the .inhom, .xarr and .ints files
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
        CALL HMFWT( NH, MHT, MHL, MHM, MHP2, NDCS, MHIBC, MHOBC,
     1            LU, FNAME )
C
C Write out the eigenvectors
C
        FNAME(ILEN+8:ILEN+13) = '.inhom'
        FNAME = FNAME(1:ILEN+13)
        CALL SVFWT( INARR, LU, IFORM, VEC2, FNAME )
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
        CALL SVPRNT ( VEC, NR, INARR, DPARR, XARR, MHT, MHL,
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
      CALL FCLOSE( LUCC   , FNCC   , 'Error' )
      CALL FCLOSE( LUCB1  , FNCB1  , 'Error' )
      CALL FCLOSE( LUCB2  , FNCB2  , 'Error' )
      CALL FCLOSE( LUCD   , FNCD   , 'Error' )
      CALL FCLOSE( LUCF   , FNCF   , 'Error' )
      CALL FCLOSE( LUCG   , FNCG   , 'Error' )
      CALL FCLOSE( LUCH   , FNCH   , 'Error' )
      CALL FCLOSE( LUCI   , FNCI   , 'Error' )
      CALL FCLOSE( LURES  , FNRES  , 'Error' )
      CALL FCLOSE( LURI   , FNRI   , 'Error' )
      CALL FCLOSE( LURO   , FNRO   , 'Error' )
      CALL FCLOSE( LUBC   , FNBC   , 'Error' )
      CALL FCLOSE( LUEPSI , FNEPSI , 'Error' )
      CALL FCLOSE( LUEPSO , FNEPSO , 'Error' )
      IF ( LULOG.NE.0 .AND. LULOG.NE.6  )
     1     CALL FCLOSE ( LULOG, FNLOG, 'Error in closing file.' )
      STOP
      END
C*********************************************************************
