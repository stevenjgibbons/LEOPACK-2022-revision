C*********************************************************************
C subroutine Purely Shear Instability Independent Growth Rate Find ***
C            -      -     -           -           -      -    -    ***
C Steve Gibbons Mon Jan 24 14:52:08 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Solves the shear instability problem (Rayleigh number = 0) for a   C
C growth rate. All the large arrays are contained within and so      C
C recompilation  will be necessitated by large change of parameters. C
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
C      \curl dv/dt  = \nabla^2 \curl v                               C
C                     - CG \curl ( k \times v )                      C
C                     - REY. \curl ( v_0 . \nabla ) v                C
C                     - REY. \curl ( v . \nabla ) v_0                C
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
C         INTPAR( 6 ) = MXIT. Maximum number of iterations in IRAM.  C
C                       (Suggested value ~200 - 400)                 C
C         INTPAR( 7 ) = IMODE. Selects form of v_0.                  C
C                       (If IMODE = 1, v_0 is solid body rotation)   C
C                       (If IMODE = 2, v_0 is as defined above)      C
C         INTPAR( 8 ) = NN1: See definition of zonal velocity.       C
C         INTPAR( 9 ) = NN2: See definition of zonal velocity.       C
C         INTPAR( 10) = IICS: Selects cos/sin dependence of v_phi    C
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
C     REY       : Reynolds number - see above equations.             C
C     DPRPAR    : Array containing double precision parameters:      C
C                                                                    C
C         DPRPAR(  1 ) = RI.     Radius of inner boundary.           C
C         DPRPAR(  2 ) = RO.     Radius of outer boundary.           C
C         DPRPAR(  3 ) = CG     (coefficient of curl ( k times v )   C
C         DPRPAR(  4 ) = DRSV    Real shift for eigensolution.       C
C         DPRPAR(  5 ) = ARTOL   Convergence parameter for Arnoldi.  C
C         DPRPAR(  6 ) = DSRSV  (Output only). Suggested real shift. C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PSIGRF( NR, LH, M, GRR, GRI, REY, DPRPAR, INTPAR)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, LH, M, INTPAR( * )
      DOUBLE PRECISION GRR, GRI, REY, DPRPAR( * )
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
      INTEGER MHIBC( 3 ), MHOBC( 3 ), LARR( 3 ), INARR( 3 ),
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
      INTEGER IFORMF, INSPCF, KL, ITGN,
     1        N1, N2, NDRVS, ICLS, NCE, NCV, NEV, MXIT, IEV,
     2        IMODE, NN1, NN2, IICS, NH0
C
      INTEGER NH, ISYM, MLOW, MINC, MMAX, IVELBC, NTHP, NPHP, IH,
     1        NP
C
      LOGICAL OSBR
C
      DOUBLE PRECISION RI, RO, ZERO, CG, DIAGEL, 
     1                 X1, X2, CA, CE, ARTOL, DRSV,
     2                 DSRSV, RESTOL, DPNF, DPRARR( 2 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IFORMF = 3 )
C
      PARAMETER ( ZERO = 0.0d0, X1 = -1.0d0,
     1            X2 = 1.0d0, DIAGEL = -200.0d0 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the validity of input parameters
C
      IF ( NR.LT.10 .OR. NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine PSIGRF.'
        PRINT *,' NR = ', NR,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LH.GT.LHMAX ) THEN
        PRINT *,' Subroutine PSIGRF.'
        PRINT *,' LH    = ', LH
        PRINT *,' LHMAX = ', LHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INSPCF = INTPAR( 1 )
      IF ( INSPCF.NE.1 .AND. INSPCF.NE.2 ) THEN
        PRINT *,' Subroutine PSIGRF.'
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
        PRINT *,' Subroutine PSIGRF.'
        PRINT *,' RI = ', RI
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( RO.LE.RI ) THEN
        PRINT *,' Subroutine PSIGRF.'
        PRINT *,' RI = ', RI
        PRINT *,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NEV    = INTPAR( 2 )
      IF ( NEV.LT.1 .OR. NEV.GT.(NCVM-2) ) THEN
        PRINT *,' Subroutine PSIGRF.'
        PRINT *,' NEV  = ', NEV
        PRINT *,' NCVM = ', NCVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NCV    = INTPAR( 3 )
      IF ( NCV.GT.NCVM ) THEN
        PRINT *,' Subroutine PSIGRF.'
        PRINT *,' NCV  = ', NCV
        PRINT *,' NCVM = ', NCVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      ISYM   = INTPAR( 4 )
      IF ( ISYM.NE.1 .AND. ISYM.NE.2 ) THEN
        PRINT *,' Subroutine PSIGRF.'
        PRINT *,' ISYM  = ', ISYM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      CG     = DPRPAR( 3 )
C
      IVELBC = INTPAR( 5 )
      IF ( IVELBC.NE.1 .AND. IVELBC.NE.2 ) THEN
        PRINT *,' Subroutine PSIGRF.'
        PRINT *,' IVELBC = ', IVELBC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      MXIT   = INTPAR( 6 )
C
      DRSV   = DPRPAR( 4 )
      ARTOL  = DPRPAR( 5 )
C
      IMODE  = INTPAR( 7 )
      IF ( IMODE.NE.1 .AND. IMODE.NE.2 ) THEN
        PRINT *,' Subroutine PSIGRF.'
        PRINT *,' IMODE = ', IMODE
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NN1    = INTPAR(  8 )
      NN2    = INTPAR(  9 )
      IICS   = INTPAR( 10 )
C
      IF ( IMODE.EQ.2 ) THEN
        IF ( IICS.NE.1 .AND. IICS.NE.2 ) THEN
          PRINT *,' Subroutine PSIGRF.'
          PRINT *,' IMODE = ', IMODE,' and IICS = ', IICS
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDIF
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
      CALL VOHMSR( ISYM, NH, MLOW, MINC, MMAX, NHMAX, MHT, MHL,
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
      MHIBC( 3 ) = 1
      MHOBC( 3 ) = 1
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
        MP0( IH ) = 3
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
      PARAM(  2 ) = ZERO
      PARAM(  3 ) = ZERO
      PARAM(  4 ) = ZERO
      PARAM(  5 ) = ZERO
      PARAM(  6 ) = CG
      PARAM(  7 ) = ZERO
      PARAM(  8 ) = 1.0d0
      PARAM(  9 ) = ZERO
      PARAM( 10 ) = ZERO
      PARAM( 11 ) = REY
      PARAM( 12 ) = REY
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
      CA = 0.0d0
      CE = 1.0d0
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
      DPRPAR( 6 ) = DSRSV
C
      RETURN
      END
C*********************************************************************
