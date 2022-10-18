C*********************************************************************
C subroutine Thermal and Shear Inst. Independent Growth Rate Find ****
C            -           -     -     -           -      -    -    ****
C Steve Gibbons Sat Jan 22 14:58:07 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Solves the linear onset of convection problem for a growth rate.   C
C All the large arrays are contained within and so recompilation     C
C will be necessitated by large change of parameters.                C
C                                                                    C
C This is purely for convenience and the routine is designed for     C
C programs where optimisation over M, LH and NR is needed.           C
C For this reason it is easiest to recalculate the entire finite     C
C difference and transform arrays.                                   C
C____________________________________________________________________C
C                                                                    C
C Finds the growth rate (real part GRR, imaginary part GRI) for the  C
C linear convective and shear instability problem :-                 C
C                                                                    C
C  c_a d \Theta/ dt = CD \nabla^2 \Theta                             C
C                     + v . ( CB1 r + CB2 r^{-2} , 0 , 0 )           C
C                     - REY. CC  v_0 . Grad ( \Theta )               C
C                                                                    C
C  c_e \curl dv/dt  = CI \nabla^2 \curl v                            C
C                     - CG \curl ( k \times v )                      C
C                     + RAY \curl ( \Theta {\bm r } )                C
C                     - REY. CF \curl ( v_0 . \nabla ) v             C
C                     - REY. CF \curl ( v . \nabla ) v_0             C
C                                                                    C
C The exact nature of the flow depends upon the integer flag IMODE.  C
C If IMODE = 1, the flow is simply set as a solid body rotation      C
C                                                                    C
C  i.e. v_{\phi} = s                                                 C
C                                                                    C
C If IMODE = 2, the flow is equatorially symmetric and depends       C
C upon the parameters NN1, NN2, ICS, RI and DPNF which are input     C
C in the arrays DPRARR and INTARR (see ZVDF also).                   C
C                                                                    C
C The zonal velocity is defined by                                   C
C                                                                    C
C   v_r         = 0                                                  C
C   v_{\theta}  = 0                                                  C
C   v_{\phi}    = s^{nn1} {cos/sin} ( nn2 * \pi ( s - r_i ) )        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                 (NR is limited by NRMAX: internal parameter)       C
C     LH        : Maximum spherical harmonic degree, l.              C
C                 (LH is limited by LHMAX: internal parameter)       C
C     M         : Wavenumber in phi.                                 C
C     INTPAR    : Array containing integer parameters:               C
C                                                                    C
C         INTPAR( 1 ) = INSPCF. This flag sets radial spacing of     C
C                       grid nodes. The only options are             C
C                       (inspcf = 1 --> equally spaced nodes )       C
C                       (inspcf = 2 --> Chebyshev nodes )            C
C                                                                    C
C         INTPAR( 2 ) = NEV.    Number of requested eigenvalues.     C
C                       (GRR and GRI are only given for the e.val    C
C                       with largest real part.)                     C
C         INTPAR( 3 ) = NCV.    Length of Arnoldi factorisation.     C
C                       NCV must be atleast 2 greater than NEV.      C
C         INTPAR( 4 ) = ISYM. ( isym = 1 --> equatorial symmetry)    C
C                             ( isym = 2 --> equatorial anti-symm.)  C
C         INTPAR( 5 ) = IVELBC. Velocity boundary condition.         C
C                       (ivelbc = 1 --> no slip )                    C
C                       (ivelbc = 2 --> stress free)                 C
C         INTPAR( 6 ) = ITHEBC. Temperature boundary condition.      C
C                       (ithebc = 1 --> fixed tm inner and outer)    C
C                       (ithebc = 2 --> fixed tm inner hf outer)     C
C                       (ithebc = 3 --> fixed hf inner tm outer)     C
C         INTPAR( 7 ) = MXIT. Maximum number of iterations in IRAM.  C
C                       (Suggested value ~200 - 400)                 C
C         INTPAR( 8 ) = IMODE. Selects form of v_0.                  C
C                       (If IMODE = 1, v_0 is solid body rotation)   C
C                       (If IMODE = 2, v_0 is as defined above)      C
C         INTPAR( 9 ) = NN1: See definition of zonal velocity.       C
C         INTPAR( 10) = NN2: See definition of zonal velocity.       C
C         INTPAR( 11) = IICS: Selects cos/sin dependence of v_phi    C
C                       (iics = 1 --> cos ( nn2 * \pi ( s - r_i ) )  C
C                       (iics = 2 --> sin ( nn2 * \pi ( s - r_i ) )  C
C                                                                    C
C (Note that NN1, NN2 and IICS are not referred to if IMODE = 1.)    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     GRR       : Real part of growth rate (output)                  C
C     GRI       : Imaginary part of growth rate (output)             C
C     RAY       : Rayleigh number - see above equations.             C
C     REY       : Reynolds number - see above equations.             C
C     DPRPAR    : Array containing double precision parameters:      C
C                                                                    C
C         DPRPAR(  1 ) = RI.     Radius of inner boundary.           C
C         DPRPAR(  2 ) = RO.     Radius of outer boundary.           C
C         DPRPAR(  3 ) = CA     (coefficient of dTheta/dt)           C
C         DPRPAR(  4 ) = CB1    (coefficient of v . ( r, 0, 0 )      C
C         DPRPAR(  5 ) = CB2    (coefficient of v . ( r^{-2}, 0, 0 ) C
C         DPRPAR(  6 ) = CD     (coefficient of nabla^2 Theta)       C
C         DPRPAR(  7 ) = CE     (coefficient of curl dv/dt)          C
C         DPRPAR(  8 ) = CG     (coefficient of curl ( k times v )   C
C         DPRPAR(  9 ) = CI     (coefficient of curl nabla^2 v)      C
C         DPRPAR( 10 ) = DRSV    Real shift for eigensolution.       C
C         DPRPAR( 11 ) = ARTOL   Convergence parameter for Arnoldi.  C
C         DPRPAR( 12 ) = DSRSV  (Output only). Suggested real shift. C
C         DPRPAR( 13 ) = CC     (coefficient of v_0.grad( theta )   )C
C         DPRPAR( 14 ) = CF     (coefficient of curl ( v.\nabla ) v )C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE TSIGRF( NR, LH, M, GRR, GRI, RAY, REY, DPRPAR,
     1                   INTPAR)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, LH, M, INTPAR( * )
      DOUBLE PRECISION GRR, GRI, RAY, REY, DPRPAR( * )
C____________________________________________________________________C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, LHMAX, NHMAX, NPHMAX, NTHMAX, KLMAX,
     1        NBNDMX, NDCS, NDRVM, ISVMAX, NPMAX, LHLH2M, NCFM,
     2        NBN, NCVM
C
      PARAMETER ( NRMAX = 50, LHMAX = 64, NHMAX = 100, NBN = 3,
     1            NTHMAX = 68, NPHMAX = 128, KLMAX = (NBN+1)*NHMAX-1,
     2            NBNDMX = 3*KLMAX+1, NDCS = 4, NDRVM = 4, 
     3            ISVMAX = NRMAX*NHMAX )
      PARAMETER ( NPMAX = (LHMAX+1)*(LHMAX+2)/2, 
     1            LHLH2M = LHMAX*(LHMAX+2), NCFM = 2*NBN + 1,
     2            NCVM = 40 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER MHIBC( 4 ), MHOBC( 4 ), LARR( 4 ), INARR( 3 ),
     1        IPIV( ISVMAX ), MHT( NHMAX ), MHL( NHMAX ),
     2        MHM( NHMAX ), MHP( NHMAX ), MHTR( NHMAX ),
     3        IWORK( NCFM ), IN0( 3 ), INTARR( 3 )
C
      INTEGER MT0( NHMAX ), ML0( NHMAX ),
     1        MM0( NHMAX ), MP0( NHMAX )
C
      DOUBLE PRECISION XARR( NRMAX ),
     1                 WVEC( ISVMAX ), SBRVEC( ISVMAX ),
     2                 SVFDC( NCFM, NRMAX, NDRVM+1, NDCS ),
     3                 PARAM( 12 )
C
      DOUBLE PRECISION COEFM1( NCFM, NCFM ), WORK1( NCFM ),
     1                 COEFM2( NCFM, NCFM ), WORK2( NCFM ),
     2                 A( NBNDMX, ISVMAX ),
     3                 GAUX( NTHMAX ), GAUW( NTHMAX )
C
      DOUBLE PRECISION FTF1( 2*NPHMAX ), FTF2( 2*NPHMAX ),
     1                 FTF3( 2*NPHMAX ), QST( LHLH2M, 3 ),
     2                 VF1( NPHMAX, NTHMAX, 3 ),
     3                 VF2( NPHMAX, NTHMAX, 3 ),
     4                 VF3( NPHMAX, NTHMAX, 3 )
C
      DOUBLE PRECISION SHC( LHLH2M ), SF( NPHMAX, NTHMAX ),
     1                 VEC0( ISVMAX )
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
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER ITHEBC, IFORMF, INSPCF, KL, ITGN,
     1        N1, N2, NDRVS, ICLS, NCE, NCV, NEV, MXIT, IEV,
     2        IMODE, NN1, NN2, IICS, NH0
C
      INTEGER NH, ISYM, MLOW, MINC, MMAX, IVELBC, NTHP, NPHP, IH,
     1        NP
C
      LOGICAL OSBR
C
      DOUBLE PRECISION RI, RO, LOW, ZERO, CD, CI, CG, DIAGEL, 
     1                 X1, X2, CA, CE, CB1, CB2, ARTOL, DRSV,
     2                 DSRSV, RESTOL, CC, CF, DPNF, DPRARR( 2 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IFORMF = 3 )
C
      PARAMETER ( LOW = 1.0d-7, ZERO = 0.0d0, X1 = -1.0d0,
     1            X2 = 1.0d0, DIAGEL = -200.0d0 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the validity of input parameters
C
      IF ( NR.LT.10 .OR. NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' NR = ', NR,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LH.GT.LHMAX ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' LH    = ', LH
        PRINT *,' LHMAX = ', LHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INSPCF = INTPAR( 1 )
      IF ( INSPCF.NE.1 .AND. INSPCF.NE.2 ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' intpar(1): INSPCF = ', INSPCF
        PRINT *,' Must be either 1 or 2.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RI     = DPRPAR( 1 )
      RO     = DPRPAR( 2 )
C
      IF ( RI.LT.ZERO ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' RI = ', RI
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( RO.LE.RI ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' RI = ', RI
        PRINT *,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NEV    = INTPAR( 2 )
      IF ( NEV.LT.1 .OR. NEV.GT.(NCVM-2) ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' NEV  = ', NEV
        PRINT *,' NCVM = ', NCVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NCV    = INTPAR( 3 )
      IF ( NCV.GT.NCVM ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' NCV  = ', NCV
        PRINT *,' NCVM = ', NCVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      ISYM   = INTPAR( 4 )
      IF ( ISYM.NE.1 .AND. ISYM.NE.2 ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' ISYM  = ', ISYM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      CA     = DPRPAR( 3 )
      CB1    = DPRPAR( 4 )
      CB2    = DPRPAR( 5 )
      CD     = DPRPAR( 6 )
      CE     = DPRPAR( 7 )
      CG     = DPRPAR( 8 )
      CI     = DPRPAR( 9 )
C
      IF ( CD.LT.LOW .OR. CI.LT.LOW ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' CD = ', CD
        PRINT *,' CI = ', CI
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IVELBC = INTPAR( 5 )
      IF ( IVELBC.NE.1 .AND. IVELBC.NE.2 ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' IVELBC = ', IVELBC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      ITHEBC = INTPAR( 6 )
      IF ( ITHEBC.NE.1 .AND. ITHEBC.NE.2 .AND. ITHEBC.NE.3 ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' ITHEBC = ', ITHEBC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      MXIT   = INTPAR( 7 )
C
      DRSV   = DPRPAR( 10 )
      ARTOL  = DPRPAR( 11 )
C
      IMODE  = INTPAR( 8 )
      IF ( IMODE.NE.1 .AND. IMODE.NE.2 ) THEN
        PRINT *,' Subroutine TSIGRF.'
        PRINT *,' IMODE = ', IMODE
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NN1    = INTPAR(  9 )
      NN2    = INTPAR( 10 )
      IICS   = INTPAR( 11 )
C
      IF ( IMODE.EQ.2 ) THEN
        IF ( IICS.NE.1 .AND. IICS.NE.2 ) THEN
          PRINT *,' Subroutine TSIGRF.'
          PRINT *,' IMODE = ', IMODE,' and IICS = ', IICS
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDIF
C
      CC     = DPRPAR( 13 )
      CF     = DPRPAR( 14 )
C
C ------------------------
C INPUT PARAMETERS ARE OK.
C ------------------------
C
C itgn is used to aid the solution routine VMEPS
C
      ITGN = NR/2
C
C inspcf is valid, so we can set the radial nodes.
C
      IF ( INSPCF.EQ.1 ) THEN
        CALL ESNAAS( NR, XARR, RI, RO )
      ELSE
        CALL ZCPAAS( NR, XARR, RI, RO )
      ENDIF
C
      MLOW = M
      MINC = 1
      MMAX = M
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
C Calculate the harmonic sets for the vorticity
C
      CALL CINDSW ( NH, MHT, MHTR )
C
      INARR( 1 ) = IFORMF
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C
C Calculate the matrix dimensions
C
      KL = (NBN+1)*NH - 1
      N1 = 3*KL + 1
      N2 = NH*NR
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
C Now set up array MHP
C this is straightforward as we only have 3 types
C
      DO IH = 1, NH
        MHP( IH ) = MHT( IH )
      ENDDO
C
C Prepare velocity, v_0.
C
      INTARR( 1 ) = NN1
      INTARR( 2 ) = NN2
      INTARR( 3 ) = IICS
C
C This will calculate the maximum value of the shear function
C
      NP = 100
      IF ( IMODE.EQ.2 ) CALL ZVNFF( INTARR, NP, RI, RO, DPNF )
C
      DPRARR( 1 ) = RI
      DPRARR( 2 ) = DPNF
C
      IN0( 1 ) = IFORMF
      IN0( 2 ) = NR
C
      CALL ASFGR( VEC0, IN0, MT0, ML0, MM0, NHMAX, LH, IMODE,
     1            NTHP, DPRARR, INTARR, QST, PA, DPA, GAUX,
     2            GAUW, VF1, XARR )
C
C Set boundary flags to 'no' boundary condition
C
      NH0 = IN0( 3 )
      DO IH = 1, NH0
        MP0( IH ) = 4
      ENDDO
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
C*****************************
C*****************************
C****   Form matrix ....  ****
C*****************************
C*****************************
C
      ICLS = 1
C
      PARAM(  1 ) = ZERO
      PARAM(  2 ) = CB1
      PARAM(  3 ) = CB2
      PARAM(  4 ) = CD
      PARAM(  5 ) = ZERO
      PARAM(  6 ) = CG
      PARAM(  7 ) = RAY
      PARAM(  8 ) = CI
      PARAM(  9 ) = REY*CC
      PARAM( 10 ) = REY*CC
      PARAM( 11 ) = REY*CF
      PARAM( 12 ) = REY*CF
C
      CALL AVMATA( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1          NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR, NTHP,
     2          NPHP, MMAX, LH, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     3          FTF3, VF1, VF2, VF3, QST, SF, SHC, PARAM, ICLS,
     4          IN0, MT0, ML0, MM0, MP0, VEC0 )
C
C***********************************
C***********************************
C*****  Solve eigen system ...  ****
C***********************************
C***********************************
C
      CALL SBRRFC( NR, INARR, MHT, MHL, MHM, MHP, MHIBC, MHOBC,
     1             NCFM, XARR, SBRVEC, OSBR, WORK1, WORK2,
     2             IWORK, COEFM1 )
C
      CALL VMEPS( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1       NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR, OSBR, SBRVEC,
     2       DIAGEL, MHIBC, MHOBC, ITGN, NEV, NCV, NCVM, MXIT, NCE,
     3       DRSV, SELECT, DR, DI, D3, WORKEV, WORKD, WORKL, IPIV,
     4       ARTOL, RESID, V, W2, WVEC, CA, CE )
C
      RESTOL = 1.0d-5
      CALL EVALAS( NCE, IEV, DRSV, RESTOL, DR, DI, D3,
     1                   GRR, GRI, DSRSV )
C
      DPRPAR( 12 ) = DSRSV
C
      RETURN
      END
C*********************************************************************
