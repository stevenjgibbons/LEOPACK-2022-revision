C*********************************************************************
C subroutine Recalled coefficient adapted Vorticity Matrix ***********
C            -                            -         -      ***********
C                                     Linear Term Addition ***********
C                                     -      -    -        ***********
C Steve Gibbons Thu Jan 20 15:23:25 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Adds to a matrix A the linear terms to the vorticity equation such C
C that multiplying A by VI will give an output vector consisting     C
C of ( p, \tau, \Theta ) such that the output is                     C
C                                                                    C
C        CD \nabla^2 \Theta                                          C
C                                                                    C
C      + CI \nabla^2 \curl v                                         C
C                                                                    C
C      - CG \curl ( k \times v )                                     C
C                                                                    C
C      - CA d theta / d t                                            C
C                                                                    C
C      - CE \curl d v / d t                                          C
C                                                                    C
C      + CH \curl ( \Theta {\bm r } )                                C
C                                                                    C
C      + v . ( CB1 r + CB2 r^{-2} , 0 , 0 )                          C
C                                                                    C
C Serves exactly the same purpose as AVMLTA except that the          C
C coefficients for the Coriolis term must be pre-calculated and      C
C stored in the arrays IHNALP, IHNGAM, TVHI, CVI by a prior call     C
C to CFVICC.                                                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARR     : Int. parameter array corresponding to VI.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR      See INDFUN for details       C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     MHT      : Array length ( * ) - atleast length NH              C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. degree, l.                             C
C     MHM      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHP      : Array length ( * ) - atleast length NH              C
C                  Pointer array to finite difference coefficients.  C
C                  MPI( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     MHTR     : Array length ( * ) - atleast length NH              C
C                 Type of function in equation rows.                 C
C                 Under normal circumstances, MHTR is formed by      C
C                                                                    C
C                 CALL CINDSW( NH, MHT, MHTR )                       C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     KL        : Number of lower diagonals in matrix                C
C                                                                    C
C     NCFM      : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C      (NDRVM must be atleast 4 and NDRVS must be 4 for atleast      C
C       grid nodes IR = 2, NR - 1 ).                                 C
C                                                                    C
C     N1        : First dimension of matrix, A.                      C
C     N2        : Second dimension of matrix, A.                     C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     ICLS      : Zeros the matrix A on entry if and only if ICLS=1. C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 )                             C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C     PARAM     : Array dim ( 8 ) containing above param.s           C
C                                                                    C
C            PARAM(  1 ) = CA                                        C
C            PARAM(  2 ) = CB1                                       C
C            PARAM(  3 ) = CB2                                       C
C            PARAM(  4 ) = CD                                        C
C            PARAM(  5 ) = CE                                        C
C            PARAM(  6 ) = CG                                        C
C            PARAM(  7 ) = CH                                        C
C            PARAM(  8 ) = CI                                        C
C                                                                    C
C     CVI       : Coefficients of vector interaction. Dim ( * )      C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. ( * )        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RVMLTA( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1                   NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR,
     2                   PARAM, ICLS, NVI, IHNALP, IHNGAM, TVHI, CVI )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MHT( * ), MHL( * ), MHM( * ), MHP( * ),
     1        MHTR( * ), NBN, KL, NCFM, NDRVM, N1, N2, NDCS, ICLS,
     2        NVI, IHNALP( * ), IHNGAM( * )
      CHARACTER *(4) TVHI( * )
      DOUBLE PRECISION CVI( * ), PARAM( 8 ), XARR( NR ), A( N1, N2 ),
     1                 SVFDC( NCFM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IMF, NH, NDRVS, ICMP, ICMPO, IOP, IFORMF, NRR,
     1        ILNPOL, ILNTOR, ILNTHE, IRNPOL, IRNTOR, IRNTHE, KLE
      DOUBLE PRECISION CA, CB1, CB2, CD, CE, CG, CH, CI, FAC, ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IMF = 1
C
      CA     = PARAM( 1 )
      CB1    = PARAM( 2 )
      CB2    = PARAM( 3 )
      CD     = PARAM( 4 )
      CE     = PARAM( 5 )
      CG     = PARAM( 6 )
      CH     = PARAM( 7 )
      CI     = PARAM( 8 )
C
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C
      IF ( IFORMF.NE.3 ) THEN
        PRINT *,' Subroutine RVMLTA.'
        PRINT *,' IFORMF   = ', IFORMF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NRR.NE.NR ) THEN
        PRINT *,' Subroutine RVMLTA.'
        PRINT *,' NR       = ', NR
        PRINT *,' INARR(2) = ', NRR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( KL.NE.((NBN+1)*NH - 1) ) THEN
        PRINT *,' Subroutine RVMLTA'
        PRINT *,' KL  = ', KL
        PRINT *,' NH  = ', NH
        PRINT *,' NBN = ', NBN
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( N1.NE.(3*KL+1) .AND. N1.NE.(2*KL+1) ) THEN
        PRINT *,' Subroutine RVMLTA'
        PRINT *,' N1 = ', N1
        PRINT *,' KL = ', KL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( N1.EQ.(3*KL+1) ) THEN
        KLE = KL
      ELSE
        KLE = 0
      ENDIF
C
      IF ( N2.NE.NR*NH ) THEN
        PRINT *,' Subroutine RVMLTA'
        PRINT *,' N2 = ', N2
        PRINT *,' NR = ', NR
        PRINT *,' NH = ', NH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Clear matrix if requested
C
      IF ( ICLS.EQ.1 ) THEN
        IOP = 0
        CALL MATOP( A, ZERO, N1, N2, IOP )
      ENDIF
C
      ILNPOL = 2
      IRNPOL = NR - 1
      ILNTOR = 3
      IRNTOR = NR - 2
      ILNTHE = 2
      IRNTHE = NR - 1
C
      NDRVS = 4
C
C First; Laplacian of temperature
C
      ICMP  = 3
      FAC   = CD
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILNTHE, IRNTHE, SVFDC, XARR, NDCS, A,
     3           N1, N2, IMF, KL, KL, KLE )
C
C Curl of Laplacian for velocity
C
      FAC = CI
      CALL AMLC( NR, INARR, MHT, MHL, MHM, MHP, MHTR, MHL,
     1           MHM, FAC, NBN, NDRVS, NDRVM, ILNPOL, IRNPOL,
     2           ILNTOR, IRNTOR, XARR, NCFM, A, N1, N2, IMF,
     3           KL, KL, KLE, NDCS, SVFDC )
C
C Curl ( k x V )
C
      FAC = (-1.0d0)*CG
      CALL RMCCFA( NR, INARR, MHT, MHL, MHM, MHP, MHTR, MHL, MHM,
     1      FAC, NBN, NDRVS, NDRVM, ILNPOL, IRNPOL, ILNTOR, IRNTOR,
     2      XARR, NCFM, A, N1, N2, IMF, KL, KL, KLE, NDCS, SVFDC,
     3      NVI, IHNALP, IHNGAM, TVHI, CVI )
C
C d theta / dt
C
      FAC = (-1.0d0)*CA
      ICMP = 3
      CALL AMTA( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1           MHTR, MHL, MHM, ICMP, FAC, ILNTHE, IRNTHE, A,
     2           N1, N2, IMF, KL, KL, KLE, NDRVS, NDRVM, NCFM,
     3           NDCS, SVFDC, XARR, NBN )
C
C curl dv/dt
C
      FAC   = (-1.0d0)*CE
      ICMP  = 1
      ICMPO = 2
      CALL AMCL( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, ICMPO, FAC, NBN, NDRVS, NDRVM,
     2           ILNTOR, IRNTOR, XARR, NCFM, SVFDC, A, N1, N2,
     3           IMF, KL, KL, KLE, NDCS )
C
      FAC   = (-1.0d0)*CE
      ICMP  = 2
      ICMPO = 1
      CALL AMCL( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, ICMPO, FAC, NBN, NDRVS, NDRVM,
     2           ILNPOL, IRNPOL, XARR, NCFM, SVFDC, A, N1, N2,
     3           IMF, KL, KL, KLE, NDCS )
C
C buoyancy term
C
      FAC   = CH
      ICMP  = 3
      ICMPO = 2
      CALL AMTA( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1           MHTR, MHL, MHM, ICMPO, FAC, ILNTOR, IRNTOR, A,
     2           N1, N2, IMF, KL, KL, KLE, NDRVS, NDRVM, NCFM,
     3           NDCS, SVFDC, XARR, NBN )
C
C heat source terms
C
      FAC = 1.0d0
      CALL AMHST( NR, INARR, MHT, MHL, MHM, MHP, MHTR,
     1            MHL, MHM, FAC, NBN, NDRVS, NDRVM, ILNTHE, IRNTHE,
     2            XARR, NCFM, SVFDC, A, N1, N2, CB1, CB2,
     3            IMF, KL, KL, KLE, NDCS )
C
      RETURN
      END
C*********************************************************************
