C*********************************************************************
C subroutine Linear Convection Growth Rate Convergence Check *********
C            -      -          -      -    -           -     *********
C Steve Gibbons Mon Jan 24 19:48:49 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Takes a set of defining parameters for the linear convective       C
C instability problem                                                C
C                                                                    C
C  c_a d \Theta/ dt = CD \nabla^2 \Theta                             C
C                     + v . ( CB1 r + CB2 r^{-2} , 0 , 0 )           C
C                                                                    C
C  c_e \curl dv/dt  = CI \nabla^2 \curl v                            C
C                     - CG \curl ( k \times v )                      C
C                     + RAY \curl ( \Theta {\bm r } )                C
C                                                                    C
C and calculates how fine the resolution needs to be for the growth  C
C rate to be converged. An initial lateral and radial resolution is  C
C supplied by the two numbers LH1 and NH1.                           C
C                                                                    C
C A character string called CHRES is then supplied which must be     C
C either 'LH' or 'NR'. This determines whether convergence is to be  C
C determined laterally or radially.                                  C
C                                                                    C
C If CHRES = 'LH' then LH is first set to LH1 and NR to NR1.         C
C The values GRR, GRI and DNORM = DSQRT( GRR**2 + GRI**2 )           C
C are then evaluated by the routine LOIGRF. One of these values      C
C (call it VAL - selected by SHVAL) is then chosen as the property   C
C of the growth rate for a given truncation.                         C
C                                                                    C
C LCGRCC then tries different values of either LH or NR (denote      C
C selected variable with LHNR) in order to find a critical number    C
C of LHNR (denoted ICRP) such that                                   C
C                                                                    C
C  VAL( ICRP ) - VAL( ICRP + INC )                                   C
C  ------------------------------- = ERR  .lt.  ERRTOL               C
C           VAL( ICRP )                                              C
C                                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH1       : Initial maximum spherical harmonic degree, l.      C
C     NR1       : Initial number of radial grid nodes.               C
C     M         : Wavenumber in phi.                                 C
C     INTPAR    : Array containing integer parameters for LOIGRF.    C
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
C                                                                    C
C     INC       : Size of increment for either LH or NR, depending   C
C                 upon the setting of CHRES (hereafter denoted LHNR) C
C                                                                    C
C     MAXP      : The maximum value LHNR can achieve. It is sensible C
C                 to set this to either LHMAX or NRMAX as defined    C
C                 in LOIGRF to avoid a breakdown of the program.     C
C                                                                    C
C     ICRP      : The lowest value of LHNR such that                 C
C                 ERR (defined above) is less than ERRTOL.           C
C                 If convergence is never achieved (IERR = 1), then  C
C                 ICRP is returned with the highest value of LHNR    C
C                 which was possible with ICRP + INC .le. MAXP       C
C                                                                    C
C     LULOG     : If this is zero, no output is written.             C
C                 Otherwise, the growth rate is output with the      C
C                 truncation at each attempt.                        C
C                                                                    C
C     IERR      : = 0 for normal exit.                               C
C                 = 1 means that within the bounds given for LHNR    C
C                   by MAXP, convergence was never achieved.         C
C                   ERR is returned as the best achieved.            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DPRPAR    : Parameters for LOIGRF.                             C
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
C                                                                    C
C     ERRTOL    : Acceptable limit for ERR.                          C
C                                                                    C
C     CRGRR     : Growth rate achieved (real part) for LHNR = ICRP.  C
C     CRGRI     : Growth rate achieved (imag part) for LHNR = ICRP.  C
C                                                                    C
C     ERR       : See above for definition. If convergence is not    C
C                 achieved, ERR is simply returned as ERR with       C
C                 LHNR at a maximum.                                 C
C                                                                    C
C     RAY       : Rayleigh number.                                   C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHRES     : 'NR' selects LHNR = NR                             C
C                 'LH' selects LHNR = LH                             C
C                                                                    C
C     CHVAL     : 'GRR' selects VAL = GRR                            C
C                 'GRI' selects VAL = GRI                            C
C                 'DNM' selects VAL = DNORM                          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LCGRCC( LH1, NR1, M, DPRPAR, INTPAR, INC, ERRTOL,
     1                   MAXP, CRGRR, CRGRI, ICRP, CHRES, CHVAL,
     2                   LULOG, IERR, ERR, RAY )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH1, NR1, M, INTPAR( * ), INC, MAXP, ICRP, LULOG, IERR
      DOUBLE PRECISION DPRPAR( * ), ERRTOL, CRGRR, CRGRI, ERR, RAY
      CHARACTER *(2) CHRES
      CHARACTER *(3) CHVAL
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER LH, NR, IJUNK
      DOUBLE PRECISION LOW, GRR, GRI, DNORM, DLAPY2, VALLO,
     1                 VALHI, OLDGRR, OLDGRI
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      IF ( CHVAL.NE.'GRR' .AND. CHVAL.NE.'GRI' .AND.
     1     CHVAL.NE.'DNM' ) THEN
        PRINT *,' Subroutine LCGRCC.'
        PRINT *,' CHVAL = ', CHVAL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( CHRES.NE.'NR' .AND. CHRES.NE.'LH' ) THEN
        PRINT *,' Subroutine LCGRCC.'
        PRINT *,' CHRES = ', CHRES
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( ERRTOL.LT.LOW ) THEN
        PRINT *,' Subroutine LCGRCC.'
        PRINT *,' ERRTOL = ', ERRTOL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( INC.LT.1 ) THEN
        PRINT *,' Subroutine LCGRCC.'
        PRINT *,' INC = ', INC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( ( CHRES.EQ.'NR' .AND. MAXP.LT.(NR1+INC) ) .OR.
     1     ( CHRES.EQ.'LH' .AND. MAXP.LT.(LH1+INC) ) ) THEN
        PRINT *,' Subroutine LCGRCC.'
        PRINT *,' CHRES = ', CHRES,' MAXP = ', MAXP
        PRINT *,' INC   = ', INC,' NR1 = ', NR1,' LH1 = ', LH1
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Do initial run at LH1 and NR1
C     .
      LH = LH1
      NR = NR1
      CALL LOIGRF( NR, LH, M, GRR, GRI, RAY, DPRPAR, INTPAR)
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 91 ) LH, NR, M, GRR, GRI
      DNORM = DLAPY2( GRR, GRI )
      IF ( CHVAL.EQ.'GRR' ) VALLO = GRR
      IF ( CHVAL.EQ.'GRI' ) VALLO = GRI
      IF ( CHVAL.EQ.'DNM' ) VALLO = DNORM
C     .
C     . Do second run with increased resolution
C     .
      IF ( CHRES.EQ.'NR' ) NR = NR + INC
      IF ( CHRES.EQ.'LH' ) LH = LH + INC
      CALL LOIGRF( NR, LH, M, GRR, GRI, RAY, DPRPAR, INTPAR)
      IF ( LULOG.NE.0 ) WRITE ( LULOG, 91 ) LH, NR, M, GRR, GRI
      DNORM = DLAPY2( GRR, GRI )
      IF ( CHVAL.EQ.'GRR' ) VALHI = GRR
      IF ( CHVAL.EQ.'GRI' ) VALHI = GRI
      IF ( CHVAL.EQ.'DNM' ) VALHI = DNORM
C     .
      CALL CHRSEF( VALLO, VALHI, IJUNK, ERR )
C     .
C     . We now know ERR for the first run.
C     . Either ERR is less than or greater than ERRTOL.
C     . In the first case, LH1,NR1 is adequate resolution
C     . and in the second case, it is not
C     .
      IF ( DABS( ERR ).LT.ERRTOL ) THEN
C       .
C       . We are converged so let's decrease resolution
C       .
        NR = NR1
        LH = LH1
 40     CONTINUE
        VALHI = VALLO
        IF ( CHRES.EQ.'NR' ) NR = NR - INC
        IF ( CHRES.EQ.'LH' ) LH = LH - INC
        OLDGRR = GRR
        OLDGRI = GRI
        CALL LOIGRF( NR, LH, M, GRR, GRI, RAY, DPRPAR, INTPAR)
        IF ( LULOG.NE.0 ) WRITE ( LULOG, 91 ) LH, NR, M, GRR, GRI
        DNORM = DLAPY2( GRR, GRI )
        IF ( CHVAL.EQ.'GRR' ) VALLO = GRR
        IF ( CHVAL.EQ.'GRI' ) VALLO = GRI
        IF ( CHVAL.EQ.'DNM' ) VALLO = DNORM
        CALL CHRSEF( VALLO, VALHI, IJUNK, ERR )
        IF ( DABS( ERR ).LT.ERRTOL ) THEN
          GOTO 40
        ELSE
          IF ( CHRES.EQ.'NR' ) ICRP = NR + INC
          IF ( CHRES.EQ.'LH' ) ICRP = LH + INC
          IERR = 0
          GOTO 60
        ENDIF
C       .
      ELSE
C       .
C       . We are not converged so let's increase resolution
C       .
 41     CONTINUE
        VALLO = VALHI
        IF ( CHRES.EQ.'NR' ) NR = NR + INC
        IF ( CHRES.EQ.'LH' ) LH = LH + INC
        IF ( (CHRES.EQ.'NR' .AND. NR.GT.MAXP ) .OR.
     1       (CHRES.EQ.'LH' .AND. LH.GT.MAXP ) ) THEN
          IERR = 1
          RETURN
        ENDIF
        OLDGRR = GRR
        OLDGRI = GRI
        CALL LOIGRF( NR, LH, M, GRR, GRI, RAY, DPRPAR, INTPAR)
        IF ( LULOG.NE.0 ) WRITE ( LULOG, 91 ) LH, NR, M, GRR, GRI
        DNORM = DLAPY2( GRR, GRI )
        IF ( CHVAL.EQ.'GRR' ) VALHI = GRR
        IF ( CHVAL.EQ.'GRI' ) VALHI = GRI
        IF ( CHVAL.EQ.'DNM' ) VALHI = DNORM
        CALL CHRSEF( VALLO, VALHI, IJUNK, ERR )
        IF ( DABS( ERR ).LT.ERRTOL ) THEN
          IF ( CHRES.EQ.'NR' ) ICRP = NR - INC
          IF ( CHRES.EQ.'LH' ) ICRP = LH - INC
          IERR = 0
          GOTO 60
        ELSE
          GOTO 41
        ENDIF
C       .
      ENDIF
C     .
 60   CONTINUE
      CRGRR = OLDGRR
      CRGRI = OLDGRI
C     .
      RETURN
 91   FORMAT('LH: ',I4,' NR: ',I4,' M: ',I4,' (',1PD16.7,',',
     1       1PD16.7,')')
      END
C*********************************************************************

