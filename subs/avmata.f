C*********************************************************************
C subroutine Adapted Vorticity Matrix All Term Addition **************
C            -       -         -      -   -    -        **************
C Steve Gibbons Thu Nov 18 19:00:33 GMT 1999                         C
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
C      - CC1 v_0 . Grad ( theta )                                    C
C                                                                    C
C      - CC2 v . Grad ( theta_0 )                                    C
C                                                                    C
C      - CF1 \curl ( v_0 . \nabla ) v                                C
C                                                                    C
C      - CF2 \curl ( v . \nabla ) v_0                                C
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
C     NTHP      : The number of theta points.                        C
C     NPHP      : The number of phi points.                          C
C                                                                    C
C     MMAX      : Maximum sph. harmonic order, m.                    C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                                                                    C
C     ICLS      : Zeros the matrix A on entry if and only if ICLS=1. C
C                                                                    C
C     IN0       : Equivalent to INARR for existing solution ( * ).   C
C     MT0       : Equivalent to MHT for existing solution ( * ).     C
C     ML0       : Equivalent to MHL for existing solution ( * ).     C
C     MM0       : Equivalent to MHM for existing solution ( * ).     C
C     MP0       : Equivalent to MHP for existing solution ( * ).     C
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
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     FTF1      : Work array - dim. (2*NPHP)                         C
C     FTF2      : Work array - dim. (2*NPHP)                         C
C     FTF3      : Work array - dim. (2*NPHP)                         C
C     VF1       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF2       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF3       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     QST       : Work array - dim. ( LH*(LH+2), 3)                  C
C     SF        : Work array - dim. ( NPHP, NTHP )                   C
C     SHC       : Work array - dim. ( LH*(LH+2) )                    C
C                                                                    C
C     PARAM     : Array dim ( * )  containing above param.s          C
C                                                                    C
C            PARAM(  1 ) = CA                                        C
C            PARAM(  2 ) = CB1                                       C
C            PARAM(  3 ) = CB2                                       C
C            PARAM(  4 ) = CD                                        C
C            PARAM(  5 ) = CE                                        C
C            PARAM(  6 ) = CG                                        C
C            PARAM(  7 ) = CH                                        C
C            PARAM(  8 ) = CI                                        C
C            PARAM(  9 ) = CC1                                       C
C            PARAM( 10 ) = CC2                                       C
C            PARAM( 11 ) = CF1                                       C
C            PARAM( 12 ) = CF2                                       C
C                                                                    C
C     VEC0      : Vector containing existing solution. Dim ( * ).    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AVMATA( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN, KL,
     1          NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR, NTHP, NPHP,
     2          MMAX, LH, GAUX, GAUW, PA, DPA, FTF1, FTF2, FTF3,
     3          VF1, VF2, VF3, QST, SF, SHC, PARAM, ICLS,
     4          IN0, MT0, ML0, MM0, MP0, VEC0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MHT( * ), MHL( * ), MHM( * ), MHP( * ),
     1        MHTR( * ), NBN, KL, NCFM, NDRVM, N1, N2, NDCS, NTHP,
     2        NPHP, MMAX, LH, ICLS, IN0( * ), MT0( * ), ML0( * ),
     3        MM0( * ), MP0( * )
      DOUBLE PRECISION PARAM( * ), XARR( NR ), A( N1, N2 ),
     1                 SVFDC( NCFM, NR, NDRVM+1, NDCS ), VEC0( * )
      DOUBLE PRECISION
     1                 QST( LH*(LH+2), 3),
     2                 VF1( NPHP, NTHP, 3), VF2( NPHP, NTHP, 3),
     3                 VF3( NPHP, NTHP, 3), FTF1( 2*NPHP ),
     4                 FTF2( 2*NPHP ), FTF3( 2*NPHP )
      DOUBLE PRECISION
     1                 SF( NPHP, NTHP ), SHC( LH*(LH+2) )
      DOUBLE PRECISION
     1                 GAUX( NTHP ), GAUW( NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IMF, NH, NDRVS, ICMP, ICMPO, IOP, IFORMF, NRR, 
     1        ILNPOL, ILNTOR, ILNTHE, IRNPOL, IRNTOR, IRNTHE,
     2        KLE
      DOUBLE PRECISION CA, CB1, CB2, CD, CE, CG, CH, CI, FAC, ZERO,
     1                 CC1, CC2, CF1, CF2
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IMF = 1
C
      CA     = PARAM(  1 )
      CB1    = PARAM(  2 )
      CB2    = PARAM(  3 )
      CD     = PARAM(  4 )
      CE     = PARAM(  5 )
      CG     = PARAM(  6 )
      CH     = PARAM(  7 )
      CI     = PARAM(  8 )
      CC1    = PARAM(  9 )
      CC2    = PARAM( 10 )
      CF1    = PARAM( 11 )
      CF2    = PARAM( 12 )
C
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C
      IF ( IFORMF.NE.3 ) THEN
        PRINT *,' Subroutine AVMATA.'
        PRINT *,' IFORMF   = ', IFORMF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NRR.NE.NR ) THEN
        PRINT *,' Subroutine AVMATA.'
        PRINT *,' NR       = ', NR
        PRINT *,' INARR(2) = ', NRR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( KL.NE.((NBN+1)*NH - 1) ) THEN
        PRINT *,' Subroutine AVMATA'
        PRINT *,' KL  = ', KL
        PRINT *,' NH  = ', NH
        PRINT *,' NBN = ', NBN
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( N1.NE.(3*KL+1) .AND. N1.NE.(2*KL+1) ) THEN
        PRINT *,' Subroutine AVMATA'
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
        PRINT *,' Subroutine AVMATA'
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
     3           N1, N2, IMF, KL, KL, KL )
C
C Curl of Laplacian for velocity
C
      FAC = CI
      CALL AMLC( NR, INARR, MHT, MHL, MHM, MHP, MHTR, MHL,
     1           MHM, FAC, NBN, NDRVS, NDRVM, ILNPOL, IRNPOL,
     2           ILNTOR, IRNTOR, XARR, NCFM, A, N1, N2, IMF,
     3           KL, KL, KL, NDCS, SVFDC )
C
C Curl ( k x V )
C
      FAC = (-1.0d0)*CG
      CALL AMCCFA( NR, INARR, MHT, MHL, MHM, MHP, MHTR, MHL, MHM,
     1      FAC, NBN, NDRVS, NDRVM, ILNPOL, IRNPOL, ILNTOR, IRNTOR,
     2      LH, NTHP, NPHP, XARR, NCFM, A, N1, N2, IMF, KL, KL, KLE,
     3      NDCS, SVFDC, MMAX, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     4      FTF3, VF1, QST )
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
C Now the non-linear terms
C First, v_0 . Grad (theta)
C
      FAC = (-1.0d0)*CC1
      CALL IV0GTA( NR,N1,N2,KL,KL,KLE,IMF,ILNTHE,IRNTHE,INARR,
     1    MHT,MHL,MHM,MHP,MHTR,MHL,MHM,IN0,MT0,ML0,MM0,MP0,
     2    NBN,NDCS,NDRVS,NDRVM,NCFM,NBN,NDCS,NDRVS,NDRVM,
     3    NCFM,LH,NTHP,NPHP,MMAX,A,FAC,XARR,VEC0,SVFDC,
     4    SVFDC,GAUX,GAUW,PA,DPA,FTF1,VF1,VF2,SF,SHC )
C
C Now v . Grad (theta_0)
C
      FAC = (-1.0d0)*CC2
      CALL IVGT0A( NR,N1,N2,KL,KL,KLE,IMF,ILNTHE,IRNTHE,INARR,
     1    MHT,MHL,MHM,MHP,MHTR,MHL,MHM,IN0,MT0,ML0,MM0,MP0,
     2    NBN,NDCS,NDRVS,NDRVM,NCFM,NBN,NDCS,NDRVS,NDRVM,
     3    NCFM,LH,NTHP,NPHP,MMAX,A,FAC,XARR,VEC0,SVFDC,
     4    SVFDC,GAUX,GAUW,PA,DPA,FTF1,VF1,VF2,SF,SHC )
C
C Now \curl ( v_0 . Grad ) v
C
      FAC = (-1.0d0)*CF1
      CALL IV0CVA( NR,N1,N2,KL,KL,KLE,IMF,ILNPOL,IRNPOL,
     1    ILNTOR,IRNTOR,INARR,MHT,MHL,MHM,MHP,MHTR,MHL,MHM,IN0,
     2    MT0,ML0,MM0,MP0,NBN,NDCS,NDRVS,NDRVM,NCFM,NBN,
     3    NDCS,NDRVS,NDRVM,NCFM,LH,NTHP,NPHP,MMAX,A,FAC,
     4    XARR,VEC0,SVFDC,SVFDC,GAUX,GAUW,PA,DPA,FTF1,FTF2,
     5    FTF3,VF1,VF2,VF3,QST )
C
C Now \curl ( v . Grad ) v_0
C
      FAC = (-1.0d0)*CF2
      CALL IVCV0A( NR,N1,N2,KL,KL,KLE,IMF,ILNPOL,IRNPOL,
     1    ILNTOR,IRNTOR,INARR,MHT,MHL,MHM,MHP,MHTR,MHL,MHM,IN0,
     2    MT0,ML0,MM0,MP0,NBN,NDCS,NDRVS,NDRVM,NCFM,NBN,
     3    NDCS,NDRVS,NDRVM,NCFM,LH,NTHP,NPHP,MMAX,A,FAC,
     4    XARR,VEC0,SVFDC,SVFDC,GAUX,GAUW,PA,DPA,FTF1,FTF2,
     5    FTF3,VF1,VF2,VF3,QST )
C
      RETURN
      END
C*********************************************************************
