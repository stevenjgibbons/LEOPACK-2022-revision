C*********************************************************************
C                                                                    C
C Steve Gibbons - Boundary-Locked Steady, Non-Linear Solutions       C
C                  Instability Calculate.                            C
C                                                                    C
C   with EigenVECtors                                                C
C                                                                    C
C Wed Jun  6 16:22:08 WEST 2001                                      C
C                                                                    C
C (New version Thu Aug  3 13:21:20 BST 2000)                         C
C                                                                    C
C*********************************************************************
      PROGRAM blscnlsic_evecs
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, LHMAX, NHMAX, NPHMAX, NTHMAX, KLMAX,
     1        NBNDMX, NDCS, NDRVM, ISVMAX, NPMAX, LHLH2M, NCFM,
     2        NBN, NRUNM, NITHMX, NCVM
C
      PARAMETER ( NRMAX = 60, LHMAX = 62, NHMAX = 250, NBN = 3,
     1            NTHMAX = 64, NPHMAX = 128, KLMAX = (NBN+1)*NHMAX-1,
     2            NBNDMX = 3*KLMAX+1, NDCS = 4, NDRVM = 4,
     3            ISVMAX = NRMAX*NHMAX, NITHMX = 100 )
      PARAMETER ( NPMAX = (LHMAX+1)*(LHMAX+2)/2, 
     1            LHLH2M = LHMAX*(LHMAX+2), NCFM = 2*NBN + 1,
     2            NRUNM = 100, NCVM = 30 )
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
C 'X' arrays are for instability.
C
      INTEGER MHTX( NHMAX ), MHLX( NHMAX ), INARRX( 3 ),
     2        MHMX( NHMAX ), MHPX( NHMAX ), MHTRX( NHMAX )
C
      INTEGER          NRARR( NRUNM )
      INTEGER          ISPARR( NRUNM )
      INTEGER          LHARR( NRUNM )
      INTEGER          ISYMAR( NRUNM )
      INTEGER          MLOWA( NRUNM )
      INTEGER          MINCA( NRUNM )
      INTEGER          MMAXA( NRUNM )
      INTEGER          IOFARR( NRUNM )
      INTEGER          NFLOQA( NRUNM )
      INTEGER          NCHARR( NRUNM )
      DOUBLE PRECISION CCARR( NRUNM )
      DOUBLE PRECISION CB1ARR( NRUNM )
      DOUBLE PRECISION CB2ARR( NRUNM )
      DOUBLE PRECISION CAARR( NRUNM )
      DOUBLE PRECISION CEARR( NRUNM )
      DOUBLE PRECISION CDARR( NRUNM )
      DOUBLE PRECISION CFARR( NRUNM )
      DOUBLE PRECISION CGARR( NRUNM )
      DOUBLE PRECISION CHARR( NRUNM )
      DOUBLE PRECISION DCHARR( NRUNM )
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
C
      DOUBLE PRECISION V( ISVMAX, NCVM ), RESID( ISVMAX ),
     1                 DR( NCVM ), DI( NCVM ), D3( NCVM ),
     2                 WORKEV( 3*NCVM ), WORKD( 3*ISVMAX ),
     3                 WORKL( 3*NCVM*NCVM + 6*NCVM ),
     4                 WVEC( ISVMAX ), SBRVEC( ISVMAX ), W2( ISVMAX )
C
      LOGICAL SELECT( NCVM )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER ITHEBC, IFORMF, LU, INSPCF, KL, IRUN, IFORM,
     1        N1, N2, NDRVS, NR, LH, IOP, INDFUN, NFLOQ, NCH,
     2        N1X, N2X, KLX
C
      INTEGER NH, ISYM, MLOW, MINC, MMAX, IVELBC, NTHP, NPHP,
     1        IWRITE, IAPP, IH, MXATT, IERR, NRUNS, IR, IITH
C
      INTEGER I, ILEN, IOF, ILNR, IRNR, IZF, IREAD, NEV, NCV,
     1        L, M, ICS, IND, INDSHC, KIB, KOB, NITH, IHD,
     2        NHX, MFLOQ
C
      INTEGER  LULOG, LURES, ICH, NCE, IEV, MXIT
C
      CHARACTER *(2)  BCH
      CHARACTER *(80) FNLOG, FNAME, ROOT, FNRES
C
      DOUBLE PRECISION RI, RO, LOW, ZERO, CD, CI, CG, CAK, CBK, CCK,
     1                 X1, X2, CC, CF, CH, CB1, CB2, RAD, DCH, TOTKE,
     2                 CTOL, EPSI, EPSO, COEF, DERV( 1 ), DKE( 2 ),
     3                 DRSV, CA, CE, GRR, GRI, ARTOL, GRMAX, GIMAX
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IFORMF = 3, LU = 12, IFORM = 1, IOP = 0,
     1            IWRITE = 3, IAPP = 4, IREAD = 1, MXIT = 400 )
C
      PARAMETER ( LOW = 1.0d-7, ZERO = 0.0d0, X1 = -1.0d0,
     1            X2 = 1.0d0, ARTOL = 1.0d-9 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO I = 1, 80
        ROOT(I:I) = ' '
        FNLOG(I:I) = ' '
        FNRES(I:I) = ' '
        FNAME(I:I) = ' '
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
      LURES  = 74
C
      FNRES(1:ILEN) = ROOT(1:ILEN)
      FNAME(1:ILEN) = ROOT(1:ILEN)
C
      FNRES(ILEN+1:ILEN+4) = '.res'
C
      CALL FOPEN( LURES  , FNRES  , IAPP )
C
C Next line should contain MXATT, CTOL, DRSV, NEV, NCV
C
      PRINT *,' Enter MXATT, CTOL, DRSV, NEV, NCV '
      PRINT *,' MXATT = number of attempts allowed '
      PRINT *,' to iterate to non-linear solution.'
      PRINT *,' CTOL = convergence criterion.'
      PRINT *,' DRSV = real shift for eigenproblem.'
      PRINT *,' NEV = requested number of eigenvalues.'
      PRINT *,' NCV = length of Arnoldi factorisation.'
      PRINT *,' NCV must be at least NEV + 2.'
      PRINT *,' NCV must be less than ', NCVM
 23   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 23
      READ ( LINE, * ) MXATT, CTOL, DRSV, NEV, NCV
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
C NR, INSPCF, LH, ISYM, MLOW, MINC, MMAX, IOF
C CC, CB1, CB2, CA, CE, CD, CF, CG, CH, CI
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
     3       CCARR( NRUNS ), CB1ARR( NRUNS ), CB2ARR( NRUNS ),
     4       CAARR( NRUNS ), CEARR( NRUNS )
 25   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 25
      READ ( LINE, * ) CDARR( NRUNS ), CFARR( NRUNS ),
     1       CGARR( NRUNS ), CHARR( NRUNS ), CIARR( NRUNS ),
     2       EPSIAR( NRUNS ), EPSOAR( NRUNS ), DCHARR( NRUNS ),
     3       NCHARR( NRUNS ), NFLOQA( NRUNS )
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
       NR     = NRARR( IRUN )
       INSPCF = ISPARR( IRUN )
       LH     = LHARR( IRUN )
       ISYM   = ISYMAR( IRUN )
       MLOW   = MLOWA( IRUN )
       MINC   = MINCA( IRUN )
       MMAX   = MMAXA( IRUN )
       IOF    = IOFARR( IRUN )
       CC     = CCARR( IRUN )
       CA     = CAARR( IRUN )
       CE     = CEARR( IRUN )
       CB1    = CB1ARR( IRUN )
       CB2    = CB2ARR( IRUN )
       CD     = CDARR( IRUN )
       CF     = CFARR( IRUN )
       CG     = CGARR( IRUN )
       CH     = CHARR( IRUN )
       CI     = CIARR( IRUN )
       EPSI   = EPSIAR( IRUN )
       EPSO   = EPSOAR( IRUN )
       NFLOQ  = NFLOQA( IRUN )
       IF ( NFLOQ.GT.MINC/2 ) NFLOQ = MINC/2
       DCH    = DCHARR( IRUN )
       NCH    = NCHARR( IRUN )
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
C OK: this is the point within the RUN loop where we
C know the radial grid nodes AND the harmonic sets for
C the boundary locked flow: we can therefore output
C the respective files and need not do so again ...
C       . FNAME(1:ILEN) contains root
C 
      IF ( IRUN.GE.1 .AND. IRUN.LT.10 )
     1   WRITE (  FNAME(ILEN+1:ILEN+7), 811 ) IRUN
      IF ( IRUN.GE.10 .AND. IRUN.LT.100 )
     1   WRITE (  FNAME(ILEN+1:ILEN+7), 812 ) IRUN
      IF ( IRUN.GE.100 .AND. IRUN.LT.1000 )
     1   WRITE (  FNAME(ILEN+1:ILEN+7), 813 ) IRUN
C
C Now set up array MHP
C this is straightforward as we only have 3 types
C
      DO IH = 1, NH
        MHP( IH ) = MHT( IH )
      ENDDO
C
C Write out integers file:
C
C                              8901234567
      FNAME(ILEN+8:ILEN+17) = '.main.ints'
      FNAME                 = FNAME(1:ILEN+17)
      CALL HMFWT( NH, MHT, MHL, MHM, MHP, NDCS, MHIBC, MHOBC,
     1            LU, FNAME )
C
C Write out radial node data
C
      FNAME(ILEN+8:ILEN+17) = '.main.xarr'
      FNAME                 = FNAME(1:ILEN+17)
      CALL XARRWT( NR, XARR, LU, FNAME, IFORM )
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
C Now we loop ICH from 0 to NCH
C and repeat the calculation at
C many values of CH.
C*****************************
      DO ICH = 0, NCH
C       :----------------------
      PARAM(  1 ) = CA
      PARAM(  2 ) = CB1
      PARAM(  3 ) = CB2
      PARAM(  4 ) = CC 
      PARAM(  5 ) = CD 
      PARAM(  6 ) = CE 
      PARAM(  7 ) = CF 
      PARAM(  8 ) = CG 
      PARAM(  9 ) = CH + DCH*DBLE( ICH )
      PARAM( 10 ) = CI    
      PARAM( 12 ) = CTOL
C
      CALL BLCNRS( VEC, INARR, NR, PARAM, MXATT, IERR, A, N1,
     1  N2, KL, MHT, MHL, MHM, MHP, MHTR, NBN, NCFM, NDRVM, NDCS,
     2  NTHP, NPHP, MMAX, LH, MHIBC, MHOBC, IPIV, LULOG, SVFDC, XARR,
     3  GAUX, GAUW, PA, DPA, FDCM, MHI, CAFIT )
C
C Calculate kinetic energy of solution
C
      TOTKE = 0.0d0
      DO IH = 1, NH
        CALL SHKEER( IH, NDCS, NR, INARR, MHT, MHL, MHP, NBN,
     1               NDRVM, NDRVM, NCFM, VEC, XARR, DKE, SVFDC )
        TOTKE = TOTKE + DKE( 1 )
      ENDDO
C
 701  FORMAT(I4)
 702  FORMAT(I1)
 703  FORMAT(1PD16.8)
 704  FORMAT(I2,I2)
 705  FORMAT(1PD16.8)
 706  FORMAT(I4,I4,I4)
 707  FORMAT(1PD16.8,1PD16.8)
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
        IF ( ICH.GE.0 .AND. ICH.LT.10 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 814 ) ICH
        IF ( ICH.GE.10 .AND. ICH.LT.100 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 815 ) ICH
        IF ( ICH.GE.100 .AND. ICH.LT.1000 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 816 ) ICH
C
C Write out the harmonic integers file
C
c       FNAME(ILEN+14:ILEN+18) = '.ints'
c       FNAME = FNAME(1:ILEN+18)
c       CALL HMFWT( NH, MHT, MHL, MHM, MHP, NDCS, MHIBC, MHOBC,
c    1            LU, FNAME )
C
C Write out the eigenvectors
C
        FNAME(ILEN+14:ILEN+18) = '.vecs'
        FNAME = FNAME(1:ILEN+18)
        CALL SVFWT( INARR, LU, IFORM, VEC, FNAME )
C
C Write out radial node data
C
c       FNAME(ILEN+14:ILEN+18) = '.xarr'
c       FNAME = FNAME(1:ILEN+18)
c       CALL XARRWT( NR, XARR, LU, FNAME, IFORM )
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
        IF ( ICH.GE.0 .AND. ICH.LT.10 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 814 ) ICH
        IF ( ICH.GE.10 .AND. ICH.LT.100 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 815 ) ICH
        IF ( ICH.GE.100 .AND. ICH.LT.1000 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 816 ) ICH
C
C Write out the harmonic integers file
C
        FNAME(ILEN+14:ILEN+22) = '.inh.ints'
        FNAME = FNAME(1:ILEN+22)
        CALL HMFWT( NH, MHT, MHL, MHM, MHP2, NDCS, MHIBC, MHOBC,
     1            LU, FNAME )
C
C Write out the eigenvectors
C
        FNAME(ILEN+14:ILEN+19) = '.inhom'
        FNAME = FNAME(1:ILEN+19)
        CALL SVFWT( INARR, LU, IFORM, VEC2, FNAME )
C
C Write out radial node data
C
c       FNAME(ILEN+14:ILEN+18) = '.xarr'
c       FNAME = FNAME(1:ILEN+18)
c       CALL XARRWT( NR, XARR, LU, FNAME, IFORM )
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
        IF ( ICH.GE.0 .AND. ICH.LT.10 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 814 ) ICH
        IF ( ICH.GE.10 .AND. ICH.LT.100 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 815 ) ICH
        IF ( ICH.GE.100 .AND. ICH.LT.1000 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 816 ) ICH
C
C Write out the harmonic integers file
C
        FNAME(ILEN+14:ILEN+18) = '.disp'
        FNAME = FNAME(1:ILEN+18)
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
 814  FORMAT('.ch00',I1)
 815  FORMAT('.ch0',I2)
 816  FORMAT('.ch',I3)
C
C Now we must calculate the linear stability of our
C boundary locked solution.
C
       DO MFLOQ = 0, NFLOQ
C        .
C        . Must calculate the harmonic sets.
C        .
         CALL VTFHSR( ISYM, NHX, MLOW, MINC, MMAX, NHMAX, MHTX, MHLX,
     1             MHMX, LH, LHMAX, MFLOQ )
C
         CALL CINDSW ( NHX, MHTX, MHTRX )
         PRINT *,' MFLOQ = ', MFLOQ,': Have ', NHX,' harmonics. '
C
         INARRX( 1 ) = IFORMF
         INARRX( 2 ) = NR
         INARRX( 3 ) = NHX
C
         KLX = (NBN+1)*NHX - 1
         N1X = 3*KLX + 1
         N2X = NHX*NR
C        .
         DO IH = 1, NHX
           MHPX( IH ) = MHTX( IH )
         ENDDO
C        .
C        . Write out instability harmonic files
C        .
         IF ( IOF.NE.0 .AND. ICH.EQ.0 ) THEN
C          .
C          . FNAME(1:ILEN) contains root
C          . FNAME(ILEN+1:ILEN+7) contains '.run000'
C          .
           IF ( MFLOQ.GE.0 .AND. MFLOQ.LT.10 )
     1        WRITE (  FNAME(ILEN+8:ILEN+21), 817 ) MFLOQ
           IF ( MFLOQ.GE.10 .AND. MFLOQ.LT.11 )
     1        WRITE (  FNAME(ILEN+8:ILEN+21), 818 ) MFLOQ
           FNAME = FNAME(1:ILEN+21)
C          .
C             890   1  2345678901
 817  FORMAT('.M0',I1,'.inst.ints')
 818  FORMAT('.M',I2,'.inst.ints')
C          .
           CALL HMFWT( NHX, MHTX, MHLX, MHMX, MHPX, NDCS,
     1               MHIBC, MHOBC, LU, FNAME )
C          .
         ENDIF
C        .
         CALL ITSLSR( INARR, MHT, MHL, MHM, MHI, MHP, NR, N1X, N2X,
     1  KLX, INARRX, MHTX, MHLX, MHMX, MHPX, MHTRX, NBN, NCFM, NDRVM,
     2      NDCS, NTHP, NPHP, MMAX, LH, MHIBC, MHOBC, IPIV, NEV, NCV,
     3   NCVM, MXIT, VEC, PARAM, A, SVFDC, XARR, GAUX, GAUW, PA, DPA,
     4  CAFIT, SBRVEC, DR, DI, D3, DRSV, ARTOL, WORKEV, WORKD, RESID,
     5               WVEC, WORKL, V, W2, SELECT, NCE, IEV, GRR, GRI )
C
C Ouput instability eigenvectors: first the harmonic files
C
         DO IH = 1, NCE
           IF ( IOF.NE.0 ) THEN
C       .
C       . FNAME(1:ILEN) contains root
C       . FNAME(ILEN+1:ILEN+7) contains '.run000'
C       .
        IF ( ICH.GE.0 .AND. ICH.LT.10 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 814 ) ICH
        IF ( ICH.GE.10 .AND. ICH.LT.100 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 815 ) ICH
        IF ( ICH.GE.100 .AND. ICH.LT.1000 )
     1     WRITE (  FNAME(ILEN+8:ILEN+13), 816 ) ICH
        IF ( MFLOQ.GE.0 .AND. MFLOQ.LT.10 )
     1     WRITE (  FNAME(ILEN+14:ILEN+22), 821 ) MFLOQ
        IF ( MFLOQ.GE.10 .AND. MFLOQ.LT.11 )
     1     WRITE (  FNAME(ILEN+14:ILEN+22), 822 ) MFLOQ
        IF ( IH.GE.0 .AND. IH.LT.10 )
     1     WRITE (  FNAME(ILEN+23:ILEN+27), 819 ) IH
        IF ( IH.GE.10 .AND. IH.LT.11 )
     1     WRITE (  FNAME(ILEN+23:ILEN+27), 820 ) IH
        FNAME = FNAME(1:ILEN+27)
C       .
        CALL EVECEX( N2, NCE, IH, V, SBRVEC )
        CALL SVFWT( INARRX, LU, IFORM, SBRVEC, FNAME )
C       .
      ENDIF
C
C             3   4  567    
 819  FORMAT('0',I1,'.sv')
 820  FORMAT(I2,'.sv')
C             456   7  89012
 821  FORMAT('.M0',I1,'.inst')
 822  FORMAT('.M',I2,'.inst')
C
C Display eigenvalues to LULOG
C
          WRITE ( LULOG, 18 ) MFLOQ, IH, DR( IH ), DI( IH )
          WRITE ( LURES, 912 ) MFLOQ, IH, DR( IH ), DI( IH ),
     1         CH + DCH*DBLE( ICH ), EPSO, TOTKE
         ENDDO
 17   FORMAT('----------------------------------------------')
 18   FORMAT('Mfl: ',I2,' Eval ',I2,' (',f16.7,',',f16.7,')')
 912  FORMAT(I2,I2,1PD16.8,1PD16.8,1PD16.8,1PD16.8,1PD16.8)
C
         PRINT *,' GRR = ', GRR,' GRI = ', GRI
         IF ( MFLOQ.EQ.0 .OR. GRR.GT.GRMAX ) THEN
c          MMFLOQ = MFLOQ
           GRMAX  = GRR
           GIMAX  = GRI
         ENDIF
         PRINT *,' GRMAX = ', GRMAX,' GIMAX = ', GIMAX
       WRITE ( 6, 708 ) CH + DCH*DBLE( ICH ), TOTKE, MFLOQ,
     1        GRR, GRI, EPSO
       ENDDO
c      WRITE ( LURES , 708 ) CH + DCH*DBLE( ICH ), TOTKE, MMFLOQ,
c    1        GRMAX, GIMAX, EPSO
 708    FORMAT(1PD16.8,1PD16.8,I3,1PD16.8,1PD16.8,1PD16.8)
C        .
C      .     end of MFLOQ loop
C
       ENDDO
C      .     end of ICH loop
      ENDDO
C     .     end of IRUN loop
      PRINT *,' Program finishing.'
      CALL FCLOSE( LURES  , FNRES  , 'Error' )
      IF ( LULOG.NE.0 .AND. LULOG.NE.6  )
     1     CALL FCLOSE ( LULOG, FNLOG, 'Error in closing file.' )
      STOP
      END
C*********************************************************************
