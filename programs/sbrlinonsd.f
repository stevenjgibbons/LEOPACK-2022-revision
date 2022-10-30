C*********************************************************************
C                                                                    C
C Steve Gibbons - Tue Nov  6 11:56:48 WET 2001                       C
C                                                                    C
C Solid Body Rotation Linear Onset Drifting frame solve              C
C -     -    -        -      -     -                                 C
C                                                                    C
C*********************************************************************
      PROGRAM sbrlinonsd
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, LHMAX, NHMAX, NPHMAX, NTHMAX, KLMAX,
     1        NBNDMX, NDCS, NDRVM, ISVMAX, NPMAX, LHLH2M, NCFM,
     2        NBN, NCVM, NRUNM, NH0MAX, IS0MAX
C
      PARAMETER ( NRMAX = 60, LHMAX = 62, NHMAX = 160, NBN = 3,
     1            NTHMAX = 64, NPHMAX = 128, KLMAX = (NBN+1)*NHMAX-1,
     2            NBNDMX = 3*KLMAX+1, NDCS = 4, NDRVM = 4, 
     3            ISVMAX = NRMAX*NHMAX, NH0MAX = 1,
     4            IS0MAX = NRMAX*NH0MAX )
      PARAMETER ( NPMAX = (LHMAX+1)*(LHMAX+2)/2, 
     1            LHLH2M = LHMAX*(LHMAX+2), NCFM = 2*NBN + 1,
     2            NCVM = 24, NRUNM = 100 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER MHIBC( 4 ), MHOBC( 4 ), LARR( 4 ), INARR( 3 ),
     1        IPIV( ISVMAX ), MHT( NHMAX ), MHL( NHMAX ),
     2        MHM( NHMAX ), MHP( NHMAX ), MHTR( NHMAX ),
     3        IWORK( NCFM ), IN0( 3 )
      INTEGER MT0( NH0MAX ), ML0( NH0MAX ),
     1        MM0( NH0MAX ), MP0( NH0MAX )
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
      DOUBLE PRECISION ROTARR( NRUNM )
      DOUBLE PRECISION CHARR1( NRUNM )
      DOUBLE PRECISION CHARR2( NRUNM )
      DOUBLE PRECISION CIARR( NRUNM )
      DOUBLE PRECISION SBRAR1( NRUNM )
      DOUBLE PRECISION SBRAR2( NRUNM )
C
      DOUBLE PRECISION DPARR( 2 ), XARR( NRMAX ), 
     1                 WVEC( ISVMAX ), SBRVEC( ISVMAX ),
     2                 SVFDC( NCFM, NRMAX, NDRVM+1, NDCS ),
     3                 PARAM( 14 )
C
      DOUBLE PRECISION COEFM1( NCFM, NCFM ), WORK1( NCFM ),
     1                 COEFM2( NCFM, NCFM ), WORK2( NCFM ),
     2                 A( NBNDMX, ISVMAX ),
     3                 GAUX( NTHMAX ), GAUW( NTHMAX )
C
      DOUBLE PRECISION FTF1( 2*NPHMAX ), FTF2( 2*NPHMAX ),
     1                 FTF3( 2*NPHMAX ), QST( LHLH2M, 3 ),
     2                 SHC( LHLH2M ), SF( NPHMAX, NTHMAX )
      DOUBLE PRECISION
     1                 VF1( NPHMAX, NTHMAX, 3 ),
     2                 VF2( NPHMAX, NTHMAX, 3 ),
     3                 VF3( NPHMAX, NTHMAX, 3 )
C
      DOUBLE PRECISION PA( NPMAX, NTHMAX ), DPA( NPMAX, NTHMAX ),
     1                 DR( NCVM ), DI( NCVM ), D3( NCVM ),
     2                 WORKEV( 3*NCVM ), WORKD( 3*ISVMAX ),
     3                 WORKL( 3*NCVM*NCVM + 6*NCVM )
C
      DOUBLE PRECISION RESID( ISVMAX ), V( ISVMAX, NCVM ),
     1                 W2( ISVMAX ), VEC0( IS0MAX ), VEC0M( IS0MAX )
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
      INTEGER I, ILEN, IOF, MVAL, ILNR, IRNR, IZF, IA, IOP
C
      INTEGER  LULOG, LURES
C
      CHARACTER *(80) FNLOG, FNAME, FNRES, ROOT
C
      DOUBLE PRECISION RI, RO, LOW, ZERO, CD, CI, ROT, REY,
     1                 X1, X2, CA, CE, CH, CB1, CB2, ARTOL, DRSV,
     2                 GRR, GRI, CTOL, CH1, CH2, CC, CF, REYSB1,
     3                 REYSB2
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IFORMF = 3, LU = 12, IFORM = 1,
     1            IWRITE = 3, IAPP = 3 )
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
        FNRES(I:I) = ' '
      ENDDO
C
 80   FORMAT(A)
C
C Start reading in input file
C First line contains root only
C
      PRINT *,' Program Solid Body Rotation LINear ONSet 1.'
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
     1         LULOG.NE.44 .AND. LULOG.NE.45 ) THEN
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
      MHIBC( 4 ) = 1
      MHOBC( 4 ) = 1
C
      FNLOG(1:ILEN) = ROOT(1:ILEN)
      FNLOG(ILEN+1:ILEN+4) = '.log'
      IF ( LULOG.NE.0 ) CALL FOPEN ( LULOG, FNLOG, IWRITE )
C
      LURES = 77
C
      FNRES(1:ILEN) = ROOT(1:ILEN)
      FNAME(1:ILEN) = ROOT(1:ILEN)
C
      FNRES(ILEN+1:ILEN+4) = '.res'
C
      CALL FOPEN( LURES  , FNRES  , IAPP )
C
C Next line should contain DRSV, NEV, NCV, MXATT, CTOL, IA
C
      PRINT *,' Enter DRSV, NEV, NCV, MXATT, CTOL, IA '
      PRINT *,' DRSV is real shift value for solving e.system.'
      PRINT *,' NEV = number of requested eigenvalues.'
      PRINT *,' NCV = length of Arnoldi factorisation.'
      PRINT *,' NCV must be at least NEV + 2.'
      PRINT *,' NCV must be less than ', NCVM
      PRINT *,' MXATT = number of attempts allowed '
      PRINT *,' to iterate to Critical Rayleigh number.'
      PRINT *,' CTOL = convergence criterion.'
      PRINT *,' IA = 1 --> minimum number of theta points.'
      PRINT *,' IA = 2 --> NTHP = 1.5*LH to avoid aliassing.'
 23   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 23
      READ ( LINE, * ) DRSV, NEV, NCV, MXATT, CTOL, IA
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
      IF ( IA.NE.1 .AND. IA.NE.2 ) THEN
        PRINT *,' IA = ', IA
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now loop around until we have solved all the
C required problems ... these parameters are
C
C NR, INSPCF, LH, ISYM, MVAL, IOF, REYSB1, REYSB2
C CA, CB1, CB2, CC, CD, CE, CF, ROT, CH1, CH2, CI
C
C In this program, parameters are stored on two lines
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
     2       IOFARR( NRUNS ), SBRAR1( NRUNS ), SBRAR2( NRUNS )
 25   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 25
      READ ( LINE, * ) CAARR( NRUNS ), CB1ARR( NRUNS ),
     1       CB2ARR( NRUNS ), CDARR( NRUNS ),
     2       CEARR( NRUNS ), ROTARR( NRUNS ),
     3       CHARR1( NRUNS ), CHARR2( NRUNS ), CIARR( NRUNS )
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
       LH     = LHARR( IRUN )
       ISYM   = ISYMAR( IRUN )
       MVAL   = MVALA( IRUN )
       IOF    = IOFARR( IRUN )
       REYSB1 = SBRAR1( IRUN )
       REYSB2 = SBRAR2( IRUN )
       CA     = CAARR( IRUN )
       CB1    = CB1ARR( IRUN )
       CB2    = CB2ARR( IRUN )
       CC     = CA
       CD     = CDARR( IRUN )
       CE     = CEARR( IRUN )
       CF     = CE
       ROT    = ROTARR( IRUN )
       CH1    = CHARR1( IRUN )
       CH2    = CHARR2( IRUN )
       CI     = CIARR( IRUN )
       REY    = ZERO
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
c     CALL NNTPPF( LH, MMAX, NTHP, NPHP, NTHMAX, NPHMAX, IA )
      IF ( IA.EQ.1 ) NTHP = LH + 1
      IF ( IA.EQ.2 ) NTHP = LH + LH/2
      IF ( NTHP/2*2.NE.NTHP ) NTHP = NTHP + 1
      CALL NPHPF( MVAL, MVAL, NPHP, NPHMAX )
      CALL GAUWTS ( X1, X2, GAUX, GAUW, NTHP )
      CALL SCHNLA ( PA, DPA, GAUX, LH, NTHP )
C
C Now calculate the full harmonic sets
C
      CALL VTHMSR( ISYM, NH, MLOW, MINC, MMAX, NHMAX, MHT, MHL,
     1             MHM, LH, LHMAX )
C
      CALL CINDSW ( NH, MHT, MHTR )
c     PRINT *,' Total of ', NH,' harmonics. '
C
      INARR( 1 ) = IFORMF
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C
      KL = (NBN+1)*NH - 1
      N1 = 3*KL + 1
      N2 = NH*NR
C
      IF ( CD.LT.LOW .OR. CI.LT.LOW ) THEN
        PRINT *,' CD = ', CD
        PRINT *,' CI = ', CI
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
C Now fill the VEC0 vector and defining arrays according
C
      IN0( 1 ) = IFORMF
      IN0( 2 ) = NR
      IN0( 3 ) = 1
C
      MT0( 1 )   =  2
      ML0( 1 )   =  1
      MM0( 1 )   =  0
      MP0( 1 )   =  4
C
      IOP = 0
      CALL VECOP( VEC0, ZERO, IS0MAX, IOP )
C
C Form SVFDC matrix
C
      NDRVS = 1
      CALL SVFDCF( NR, NDCS, NBN, 1, NR, MHIBC, MHOBC,
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
      PARAM(  6 ) = ROT
      PARAM(  7 ) = CH1
      PARAM(  8 ) = CH2
      PARAM(  9 ) = CI
      PARAM( 10 ) = CTOL
      PARAM( 11 ) = CC
      PARAM( 12 ) = CC
      PARAM( 13 ) = CF
      PARAM( 14 ) = CF
C
      CALL IOCRBD( NR,INARR,MHT,MHL,MHM,MHP,MHTR,NBN,KL,NCFM,
     1   NDRVM,N1,N2,NDCS,NTHP,NPHP,MMAX,LH,IEV,MHIBC,MHOBC,NEV,NCV,
     2   NCVM,MXIT,IPIV,LULOG,SVFDC,A,XARR,GAUX,GAUW,PA,DPA,SBRVEC,
     3   RESID,W2,WVEC,FTF1,FTF2,FTF3,VF1,VF2,VF3,QST,SF,SHC,SELECT,
     4   ARTOL,DRSV,PARAM,MXATT,IERR,DR,DI,D3,WORKEV,WORKD,WORKL,V,
     5   IN0,MT0,ML0,MM0,MP0,VEC0,VEC0M,MVAL,REY,REYSB1,REYSB2 )
C
      CA  = PARAM(  1 )
      CB1 = PARAM(  2 )
      CB2 = PARAM(  3 )
      CD  = PARAM(  4 )
      CE  = PARAM(  5 )
      CH  = PARAM(  7 )
      CI  = PARAM(  8 )
      GRR = PARAM(  9 )
      GRI = PARAM( 10 )
C
c701  FORMAT(I4)
c702  FORMAT(I1)
c703  FORMAT(1PD16.8)
c704  FORMAT(I2,I2)
c705  FORMAT(1PD16.8,1PD16.8,1PD16.8,1PD16.8)
 705  FORMAT(I4,I4,I4,1PD16.8,1PD16.8,1PD14.6,
     1       1PD16.8,1PD16.8)
C
      IF ( IERR.GT.0 ) THEN
        WRITE ( LURES , 705 ) MVAL, LH, NR, 
     1                        ROT, CH, GRR, GRI, REYSB1
      ELSE
        PRINT *,' IRUN = ', IRUN,' IERR = ', IERR
        PRINT *,' Solution unsuccessful.'
      ENDIF
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
      CALL FCLOSE( LURES  , FNRES  , 'Error' )
      IF ( LULOG.NE.0 )
     1     CALL FCLOSE ( LULOG, FNLOG, 'Error in closing file.' )
      STOP
 41   FORMAT('Harmonic ',I4,' Type= ',I1,' L= ',I3,' M= ',I3 )
      END
C*********************************************************************
