C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C                                                                    C
C Dudley James  Instability Eigenvalue Problem Growth Rate Find      C
C -      -      -           -          -       -      -    -         C
C Tue Mar 13 12:56:18 MET 2001                                       C
C                                                                    C
C*********************************************************************
      PROGRAM djiepgrf
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, LHMAX, NHMAX, NPHMAX, NTHMAX, KLMAX,
     1        NBNDMX, NDCS, NDRVM, ISVMAX, NPMAX, LHLH2M, NCFM,
     2        NBN, NBN0MX, NRUNM, NCFM0M, NVMAX, IVELMX, NDCS0,
     3        MAXNVI, NCVM
C
      PARAMETER ( NRMAX = 200, LHMAX = 50, MAXNVI = 30000, NRUNM = 50,
     1            LHLH2M = LHMAX*(LHMAX+2), NBN = 2, NCVM = 25,
     2            NHMAX = LHMAX*2, NPHMAX = 64, NTHMAX = LHMAX + 2,
     3            KLMAX = (NBN+1)*NHMAX-1, NBNDMX = 3*KLMAX+1 )
C
      PARAMETER ( NBN0MX = 5, NDCS = LHMAX + 1, NDRVM = 2,
     1            NCFM = 2*NBN + 1, NCFM0M = 2*NBN0MX + 1,
     2            ISVMAX = NRMAX*NHMAX, NVMAX = 2,
     3            IVELMX = NVMAX*NRMAX, NDCS0 = 1,
     4            NPMAX = (LHMAX+1)*(LHMAX+2)/2 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER IQ( 2, LHMAX, 2*(LHMAX+1) ), IHNALP( MAXNVI ),
     1        IHNBET( MAXNVI ), IHNGAM( MAXNVI ),
     2        IWORK( NCFM ), IWORK0( NCFM0M )
C
      DOUBLE PRECISION COFM1( NCFM, NCFM ), WORK1( NCFM ),
     1                 COFM2( NCFM, NCFM ), WORK2( NCFM )
C
      DOUBLE PRECISION COFM10( NCFM0M, NCFM0M ), WORK10( NCFM0M ),
     1                 COFM20( NCFM0M, NCFM0M ), WORK20( NCFM0M )
C
      INTEGER MHBCS( NDCS ), LARR( NDCS ), INARR( 3 ),
     1        IN0( 3 ), MT0( NVMAX ), ML0( NVMAX ), MM0( NVMAX ),
     2        MP0( NVMAX ), MHBCS0( NDCS0 ), LARR0( NDCS0 ),
     3        MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ), MHP( NHMAX )
C
      DOUBLE PRECISION XARR( NRMAX ), CVI( MAXNVI ),
     1                 GAUX( NTHMAX ), GAUW( NTHMAX ), 
     2                 PA( NPMAX, NTHMAX ), DPA( NPMAX, NTHMAX )
C
      DOUBLE PRECISION VEC0( IVELMX ), FF1( 2*NPHMAX ),
     1                 FF2( 2*NPHMAX ), FF3( 2*NPHMAX ),
     2                 QST( LHLH2M, 3 )
C
      DOUBLE PRECISION VF1( NPHMAX, NTHMAX, 3 ),
     1                 VF2( NPHMAX, NTHMAX, 3 ),
     2                 VF3( NPHMAX, NTHMAX, 3 )
C
      INTEGER NRARR( NRUNM ), LHARR( NRUNM ),
     1        IOARR( NRUNM ), IPARR( NRUNM ), 
     2        NBARR( NRUNM ), IPIV( ISVMAX )
C
      DOUBLE PRECISION EPSARR( NRUNM ), LAMARR( NRUNM ), 
     1                 RMARR( NRUNM )
C
      INTEGER SFTARR( NRUNM ), SFLARR( NRUNM ), SFMARR( NRUNM )
C
      CHARACTER *(4) TVHI( MAXNVI )
C
      CHARACTER *(80) ROOT, FNLOG, FNAME
      CHARACTER *(200) LINE
C
      DOUBLE PRECISION SVFDC( NCFM, NRMAX, NDRVM+1, NDCS ),
     1                 SVFDC0( NCFM0M, NRMAX, NDRVM+1, NDCS0 )
C
      DOUBLE PRECISION RESID( ISVMAX ), V( ISVMAX, NCVM ),
     1                 W2( ISVMAX ), A( NBNDMX, ISVMAX ),
     2                 DR( NCVM ), DI( NCVM ), D3( NCVM ),
     3                 WORKEV( 3*NCVM ), WORKD( 3*ISVMAX ),
     4                 WORKL( 3*NCVM*NCVM + 6*NCVM )
C
      LOGICAL SELECT( NCVM )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER LH, NR, ITRI, IH, NH0, SFIT, SFIL, SFIM, NTHP, NPHP,
     1        NH, KL, NBN0, NCFM0, I, ILEN, NEV, NCV, NRUNS, NEG,
     2        LULOG, IWRITE, N1, N2, IRUN, IS, ILN, IRN, MXIT, NVI,
     3        NDRVS, IMF, IOF, NCE, ISP, LU, NEGRP, IFORM, NEGIP
C
      DOUBLE PRECISION RI, RO, X1, X2, CMTD, ARTOL, ZERO, DRSV,
     1                 LAMBDA, DIAGEL, FAC, EPS, RM, GRMX
C
      LOGICAL FIRST
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( X1 = -1.0d0, X2 = 1.0d0, ARTOL = 1.0d-8,
     1            ZERO = 0.0d0, IWRITE = 3, IMF = 1, IFORM = 1 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DIAGEL = -200.0d0
      LU     = 43
      LULOG  = 44
C
      DO I = 1, 80
        ROOT(I:I) = ' '
      ENDDO
C
C Start to read in input file
C First line contains ROOT only.
C
 80   FORMAT(A)
      PRINT *,' Enter ROOT '
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
      FNLOG(1:ILEN) = ROOT(1:ILEN)
      FNAME(1:ILEN) = ROOT(1:ILEN)
      FNLOG(ILEN+1:ILEN+4) = '.log'
C
C Next line should contain RI, RO, DRSV, NEV, NCV
C RI must be atleast zero
C RO must be strictly greater than RI
C
      PRINT *,' Enter RI, RO, DRSV, NEV, NCV '
      PRINT *,' DRSV is real shift value for solving e.system.'
      PRINT *,' NEV = number of requested eigenvalues.'
      PRINT *,' NCV = length of Arnoldi factorisation.'
      PRINT *,' NCV must be at least NEV + 2.'
      PRINT *,' NCV must be less than ', NCVM
 23   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 23
      READ ( LINE, * ) RI, RO, DRSV, NEV, NCV
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
      READ ( LINE, * )  NRARR( NRUNS ), LHARR( NRUNS ),
     1       SFTARR( NRUNS ), SFLARR( NRUNS ), SFMARR( NRUNS ),
     2       IOARR( NRUNS ), IPARR( NRUNS ), NBARR( NRUNS ),
     3       EPSARR( NRUNS ), LAMARR( NRUNS ), RMARR( NRUNS )
      GOTO 24
 300  CONTINUE
C
      CALL FOPEN ( LULOG, FNLOG, IWRITE )
C
C o.k. - now have read all the parameters: so let's
C loop around one by one.
C
      ITRI  = 0
      DO IRUN = 1, NRUNS
C
        NR     = NRARR( IRUN )
        LH     = LHARR( IRUN )
        SFIT   = SFTARR( IRUN )
        SFIL   = SFLARR( IRUN )
        SFIM   = SFMARR( IRUN )
        IOF    = IOARR( IRUN )
        ISP    = IPARR( IRUN )
        NBN0   = NBARR( IRUN )
        EPS    = EPSARR( IRUN )
        LAMBDA = LAMARR( IRUN )
        RM     = RMARR( IRUN )
C
C*********************************************************************
C*********************************************************************
C*****                                                           *****
C***** This is where we actually begin calculations for a given  *****
C***** run. Variables needed are NR, LH, RI, RO, ITRI, IOF,      *****
C***** ISP, EPS, LAMBDA, RM and NBN0.                            *****
C*****                                                           *****
C*********************************************************************
C*********************************************************************
C
      IF ( NR.GT.NRMAX ) THEN
        PRINT *,' NR = ',NR,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LH.GT.LHMAX ) THEN
        PRINT *,' LH = ',LH,' LHMAX = ', LHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ISP.NE.1 .AND. ISP.NE.2 ) THEN
        PRINT *,' ISP  = ', ISP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C IOF = 0 --> no output of eigenfunctions
C IOF = 1 --> output only eigenfunctions corresponding
C             to the largest real eigenvalue
C IOF = 2 --> output all eigenfunctions
C
      IF ( IOF.NE.0 .AND. IOF.NE.1 .AND. IOF.NE.2 ) THEN
        PRINT *,' IOF = ', IOF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NH0 = NVMAX 
      IN0( 1 ) = 4
      IN0( 2 ) = NR
      IN0( 3 ) = NH0
      I        = 1
      CALL DJVHMF( IN0, MT0, ML0, MM0, I )
C
      DO IH = 1, NH0
        MP0( IH ) = 1
      ENDDO
C
C Fill array with radius values
C
      IF ( ISP.EQ.1 ) THEN
        CALL ESNAAS( NR, XARR, RI, RO )
      ELSE
        CALL ZCPAAS( NR, XARR, RI, RO )
      ENDIF
C
      ILN = 1
      IRN = NR
      CALL SVDJVF( NR, IN0, MT0, ML0, MM0, ILN, IRN, VEC0,
     1             XARR, EPS, LAMBDA )
      ILN = 2
      IRN = NR - 1
C
      CALL KDTHSR( SFIT, SFIL, SFIM, MT0, ML0, MM0,
     1             MHT, MHL, MHM, IQ, LHMAX, NHMAX,
     2             ITRI, NVMAX, LH, NH, NH0 )
C
      INARR( 1 ) = 3
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C
C Prepare for finite difference schemes
C We will set boundary conditions for LH+1
C different finite difference schemes
C for IS = 1, LH we will set MHBCS( is ) = 7
C (for insulating magnetic field and LARR( is ) = is (=L)
C for IS = LH + 1, set MHBCS( is ) = 2
C (function must vanish at boundaries - this is
C the poloidal magnetic field function)
C - can set LARR( is ) to zero
C for IS = LH + 2, NDCS, simply set LARR( is ) = -1
C to be ignored.
C
      DO IS = 1, NDCS
C
C poloidal magnetic field finite diff. schemes
C
        IF ( IS.LE.LH ) THEN
          LARR( IS )  = IS
          MHBCS( IS ) = 7
        ENDIF
C
C toroidal magnetic field finite diff. scheme
C
        IF ( IS.EQ.(LH+1) ) THEN
          LARR( IS )  = 0
          MHBCS( IS ) = 2
        ENDIF
C
C unused elements of SVFDC
C
        IF ( IS.GT.(LH+1) ) LARR( IS ) = -1
C
      ENDDO
C
C Now set up array MHP
C
      DO IH = 1, NH
        IF ( MHT( IH ).EQ.4 ) MHP( IH ) = MHL( IH )
        IF ( MHT( IH ).EQ.5 ) MHP( IH ) = LH + 1
      ENDDO
C
      LARR0( 1 )  = 0
      MHBCS0( 1 ) = 1
C
      IF ( NBN0.GT.NBN0MX ) THEN
        PRINT *,' NBN0 = ',NBN0,' NBN0MX = ',NBN0MX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      NCFM0 = 2*NBN0 + 1
C
C Calculate finite difference coeff.s for VEC0
C
      NDRVS = 0
      CALL SVFDCF( NR, NDCS0, NBN0, 1, 1, MHBCS0, MHBCS0,
     1             LARR0, NCFM0, NCFM0, NDRVS, NDRVM, XARR,
     2             IWORK0, SVFDC0, COFM10, COFM20, WORK10, WORK20 )
C
      NDRVS = 2
      CALL SVFDCF( NR, NDCS0, NBN0, 2, NR-1, MHBCS, MHBCS,
     1             LARR0, NCFM0, NCFM0, NDRVS, NDRVM, XARR,
     2             IWORK0, SVFDC0, COFM10, COFM20, WORK10, WORK20 )
C
      NDRVS = 0
      CALL SVFDCF( NR, NDCS0, NBN0, NR, NR, MHBCS, MHBCS,
     1             LARR0, NCFM0, NCFM0, NDRVS, NDRVM, XARR,
     2             IWORK0, SVFDC0, COFM10, COFM20, WORK10, WORK20 )
C
C Calculate finite difference coeff.s for magnetic field
C
      NDRVS = 0
      CALL SVFDCF( NR, NDCS, NBN, 1, 1, MHBCS, MHBCS,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COFM1, COFM2, WORK1, WORK2 )
C
      NDRVS = 2
      CALL SVFDCF( NR, NDCS, NBN, 2, NR-1, MHBCS, MHBCS,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COFM1, COFM2, WORK1, WORK2 )
C
      NDRVS = 0
      CALL SVFDCF( NR, NDCS, NBN, NR, NR, MHBCS, MHBCS,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COFM1, COFM2, WORK1, WORK2 )
C
C Calculate Legendre Polynomials etc.
C
      CALL ONTPPF( LH, LH, NTHP, NPHP, NTHMAX, NPHMAX )
      CALL GAUWTS ( X1, X2, GAUX, GAUW, NTHP )
      CALL SCHNLA ( PA, DPA, GAUX, LH, NTHP )
C
      KL = (NBN+1)*NH - 1
      N1 = 3*KL + 1
      N2 = NH*NR
C
C Zero matrix A
C
      I = 0
      CALL MATOP( A, ZERO, N1, N2, I )
C
      NVI = 0
      CALL VCPCC( NVI, MAXNVI, IHNALP, IHNBET, IHNGAM, NH0,
     1    MT0, ML0, MM0, NH, MHT, MHL, MHM, NH, MHT, MHL, MHM,
     2    LH, NTHP, NPHP, LH, TVHI, CVI, GAUX, GAUW, PA, DPA,
     3    FF1, FF2, FF3, VF1, VF2, VF3, QST )
C
      FAC   = RM
      NDRVS = 2
      CALL RV0MFA( NR, N1, N2, KL, KL, KL, IMF, ILN, IRN,
     1    ILN, IRN, INARR, MHT, MHL, MHM, MHP, MHT, MHL, MHM, IN0,
     2    MT0, ML0, MM0, MP0, NBN, NDCS, NDRVS, NDRVM, NCFM, NBN0,
     3    NDCS0, NDRVS, NDRVM, NCFM0, A, FAC, XARR, VEC0, SVFDC,
     4    SVFDC0, NVI, IHNALP, IHNBET, IHNGAM, TVHI, CVI )
C
C Add diffusion coefficients
C
      FAC = 1.0d0
      I   = 4
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, I, MHT,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILN, IRN, SVFDC, XARR, NDCS, A, N1,
     3           N2, IMF, KL, KL, KL )
C
      I   = 5
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, I, MHT,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILN, IRN, SVFDC, XARR, NDCS, A, N1,
     3           N2, IMF, KL, KL, KL )
C
      CMTD = 1.0d0
      MXIT = 400
      CALL MFSEPS( NR, INARR, MHT, MHL, MHM, MHP, NBN, KL,
     1             NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR,
     2             DIAGEL, MHBCS, MHBCS, NEV, NCV, NCVM, MXIT, NCE,
     3             DRSV, SELECT, DR, DI, D3, WORKEV, WORKD, WORKL,
     4             IPIV, ARTOL, RESID, V, W2, CMTD )
C
C Output results to the log file
C
      WRITE ( LULOG, 120 ) RI, RO, NR, ISP, NBN0
 120  FORMAT('ri: ',1PD16.7,' ro: ',1PD16.7,' nr: ',I5,
     1       ' isp: ',I2,' nbn0: ',I2)
      WRITE ( LULOG, 121 ) LH, SFIT, SFIL, SFIM
 121  FORMAT('lh: ',I3,' sft: ',I3,' sfl: ',I3,' sfm: ',I3)
      WRITE ( LULOG, 122 ) EPS, LAMBDA, RM
 122  FORMAT('eps: ',1PD16.7,' La: ',1PD16.7,' rm: ',1PD16.7)
      WRITE ( LULOG, 124 )  NCE
 124  FORMAT(I4,' eigenvalues converged.')
      GRMX   = -1.0d6
      NEGRP  = -1
      NEGIP  = -1
      FIRST  = .TRUE.
      DO I = 1, NCE
        NEG = -1
        WRITE ( LULOG, 125 ) I, DR( I ), DI( I ), D3( I )
        IF ( DR( I ).GE.GRMX ) THEN
          NEG   = I
          GRMX  = DR( I )
        ENDIF
        IF ( NEG.GE.0 ) THEN
          IF ( FIRST ) THEN
            NEGRP = NEG
            NEGIP = NEG + 1
            FIRST = .FALSE.
          ELSE
            FIRST = .TRUE.
          ENDIF
        ENDIF
      ENDDO
 125  FORMAT('eval: ',I3,' (',1PD16.7,',',1PD16.7,') res: ',
     1       1PD16.7)
C
C Now prepare to output eigenvectors if they
C are requested.
C
      IF ( IOF.EQ.0 ) GOTO 500
C
C Atleast one eigenvector is requested, so
C let's output the xarr and ints files
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
        CALL HMFWT( NH, MHT, MHL, MHM, MHP, NDCS, MHBCS, MHBCS,
     1            LU, FNAME )
C
C Write out radial node data
C
        FNAME(ILEN+8:ILEN+12) = '.xarr'
        FNAME = FNAME(1:ILEN+12)
        CALL XARRWT( NR, XARR, LU, FNAME, IFORM )
C
C Write out the eigenvectors
C
        DO I = 1, NCE
          IF (    I.NE.NEGRP .AND. I.NE.NEGIP
     1               .AND. IOF.EQ.1 ) GOTO 510
C
C Let's extract the appropriate eigenvalue from
C V and put it into RESID
C
          CALL EVECEX( N2, NCE, I, V, RESID )
C
          IF ( I.GE.1 .AND. I.LT.10 )
     1     WRITE (  FNAME(ILEN+8:ILEN+12), 711 ) I
          IF ( I.GE.10 .AND. I.LT.100 )
     1     WRITE (  FNAME(ILEN+8:ILEN+12), 712 ) I
          IF ( I.GE.100 .AND. I.LT.1000 )
     1     WRITE (  FNAME(ILEN+8:ILEN+12), 713 ) I
          FNAME(ILEN+13:ILEN+17) = '.vecs'
          IF ( I.EQ.NEGRP .AND. IOF.EQ.1 )
     1           FNAME(ILEN+13:ILEN+17) = '.vecr'
          IF ( I.EQ.NEGIP .AND. IOF.EQ.1 )
     1           FNAME(ILEN+13:ILEN+17) = '.veci'
          FNAME = FNAME(1:ILEN+17)
          CALL SVFWT( INARR, LU, IFORM, RESID, FNAME )
 510    CONTINUE
        ENDDO       
C
C-------------------------------------------------------------
C End of IRUN calculation.
C-------------------------------------------------------------
 500  CONTINUE
      ENDDO
C
 711  FORMAT('.e00',I1)
 712  FORMAT('.e0',I2)
 713  FORMAT('.e',I3)
C
 811  FORMAT('.run00',I1)
 812  FORMAT('.run0',I2)
 813  FORMAT('.run',I3)
C
      CALL FCLOSE ( LULOG, FNLOG, 'Error for FNLOG.' )
      STOP
      END
C*********************************************************************
