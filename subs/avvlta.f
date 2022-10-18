C*********************************************************************
C subroutine Adapted Vorticity Vector Linear Term Addition ***********
C            -       -         -      -      -    -        ***********
C Steve Gibbons Fri Nov 19 10:37:44 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Adds to a vector RHS the linear terms to the vorticity equation    C
C such that AVVLTA will act upon VEC to give an output vector        C
C consisting of ( p, \tau, \Theta ) such that the output is          C
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
C     ILN       : Left most node for general differences (see FDCM)  C
C     IRN       : Right most node for general differences (see FDCM) C
C                                                                    C
C     NCFM      : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C      (NDRVM must be atleast 4 and NDRVS must be 4 for atleast      C
C       grid nodes IR = 2, NR - 1 ).                                 C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     NTHP      : The number of theta points.                        C
C     NPHP      : The number of phi points.                          C
C                                                                    C
C     MMAX      : Maximum sph. harmonic order, m.                    C
C                                                                    C
C     M0        : Lowest non-zero mode if simple multiplicity        C
C                  applies. (Put 1 if this is safest).               C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                                                                    C
C     ICLS      : Zeros the vec. RHS on entry if and only if ICLS=1. C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     FDCM      : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM ).                    C
C                  Formed by a call to FDCM with                     C
C                  NLMN = ILN                                        C
C                  NRMN = IRN                                        C
C                  NLMC = ILN                                        C
C                  NRMC = IRN                                        C
C                                                                    C
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
C     VF        : Work array - dim. ( NPHP, NTHP, 3)                 C
C     RQST1     : Work array - dim. ( LH*(LH+2), 3, NR ).            C
C     RQST2     : Work array - dim. ( LH*(LH+2), 3, NR ).            C
C     ZCF       : Work array - dim. ( NR ).                          C
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
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AVVLTA( NR, INARR, MHT, MHL, MHM, MHP, MHTR, NBN,
     1                   ILN, IRN, NCFM, NDRVM, SVFDC, FDCM, NDCS,
     2                   XARR, NTHP, NPHP, MMAX, M0, LH, GAUX, GAUW,
     3                   PA, DPA, FTF1, FTF2, FTF3, VF, RQST1, RQST2,
     4                   PARAM, ZCF, ICLS, VEC, RHS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MHT( * ), MHL( * ), MHM( * ), MHP( * ),
     1        MHTR( * ), NBN, ILN, IRN, NCFM, NDRVM, NDCS, NTHP,
     2        NPHP, MMAX, M0, LH, ICLS
      DOUBLE PRECISION PARAM( 8 ), XARR( NR ), RHS( * ), VEC( * ),
     1                 SVFDC( NCFM, NR, NDRVM+1, NDCS ),
     2                 FDCM( NCFM, NR, NDRVM ), ZCF( NR )
      DOUBLE PRECISION RQST1( LH*(LH+2) ,3 ,NR ),
     1                 RQST2( LH*(LH+2) ,3 ,NR ),
     2                 VF( NPHP, NTHP, 3), FTF1( 2*NPHP ),
     3                 FTF2( 2*NPHP ), FTF3( 2*NPHP )
      DOUBLE PRECISION
     1                 GAUX( NTHP ), GAUW( NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NH, NDRVS, ICMP, ICMPO, IOP, IFORMF, NRR,
     1        ILNPOL, ILNTOR, ILNTHE, IRNPOL, IRNTOR, IRNTHE,
     2        IRES, N2
      DOUBLE PRECISION CA, CB1, CB2, CD, CE, CG, CH, CI, FAC, ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IRES   = 2
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
      N2     = NR*NH
C
      IF ( IFORMF.NE.3 ) THEN
        PRINT *,' Subroutine AVVLTA.'
        PRINT *,' IFORMF   = ', IFORMF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NRR.NE.NR ) THEN
        PRINT *,' Subroutine AVVLTA.'
        PRINT *,' NR       = ', NR
        PRINT *,' INARR(2) = ', NRR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Clear rhs if requested
C
      IF ( ICLS.EQ.1 ) THEN
        IOP = 0
        CALL VECOP( RHS, ZERO, N2, IOP )
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
      CALL ASVLP( NR, NDCS, VEC, INARR, MHT, MHL, MHM, ICMP,
     1            MHP, RHS, INARR, MHTR, MHL, MHM, FAC, NBN,
     2            NDRVS, NDRVM, ILNTHE, IRNTHE, SVFDC, XARR, NCFM )
C
C Curl of Laplacian for velocity
C
      FAC = CI
      CALL ASVLC( NR, NDCS, VEC, INARR, MHT, MHL, MHM, MHP,
     1            RHS, INARR, MHTR, MHL, MHM, FAC, NBN, NDRVS,
     2            NDRVM, ILNPOL, IRNPOL, SVFDC, XARR, NCFM )
C
C Curl ( k x V )
C
      FAC = (-1.0d0)*CG
      CALL ASVCCF( NDCS, NTHP, NPHP, LH, IRES, MMAX, M0, NBN,
     1     NCFM, NR, NDRVS, NDRVM, ILN, IRN, INARR, MHT, MHL,
     2     MHM, MHP, INARR, MHTR, MHL, MHM, VEC, RHS, SVFDC, FDCM,
     3     GAUX, GAUW, PA, DPA, RQST1, RQST2, ZCF, XARR, VF,
     4     FTF1, FTF2, FTF3, FAC )
C
C d theta / dt
C
      FAC = (-1.0d0)*CA
      ICMP = 3
      CALL ASVTA( NR, NDCS, VEC, INARR, MHT, MHL, MHM, ICMP,
     1            MHP, RHS, INARR, MHTR, MHL, MHM, ICMP,
     2            FAC, ILNTHE, IRNTHE, NBN, NDRVS, NDRVM, NCFM,
     3            SVFDC )
C
C curl dv/dt
C
      FAC   = (-1.0d0)*CE
      ICMP  = 1
      ICMPO = 2
      CALL ASVCL( NR, NDCS, VEC, INARR, MHT, MHL, MHM, ICMP,
     1            MHP, RHS, INARR, MHTR, MHL, MHM, ICMPO,
     2            FAC, NBN, NDRVS, NDRVM, ILNTOR, IRNTOR, SVFDC,
     3            XARR, NCFM )
C
      FAC   = (-1.0d0)*CE
      ICMP  = 2
      ICMPO = 1
      CALL ASVCL( NR, NDCS, VEC, INARR, MHT, MHL, MHM, ICMP,
     1            MHP, RHS, INARR, MHTR, MHL, MHM, ICMPO,
     2            FAC, NBN, NDRVS, NDRVM, ILNPOL, IRNPOL, SVFDC,
     3            XARR, NCFM )
C
C buoyancy term
C
      FAC   = CH
      ICMP  = 3
      ICMPO = 2
      CALL ASVTA( NR, NDCS, VEC, INARR, MHT, MHL, MHM, ICMP,
     1            MHP, RHS, INARR, MHTR, MHL, MHM, ICMPO,
     2            FAC, ILNTHE, IRNTHE, NBN, NDRVS, NDRVM, NCFM,
     3            SVFDC )
C
C heat source terms
C
      FAC = 1.0d0
      CALL ASVHST( NR, NDCS, VEC, INARR, MHT, MHL, MHM, MHP,
     1             RHS, INARR, MHTR, MHL, MHM, FAC, NBN,
     2             NDRVS, NDRVM, ILNTHE, IRNTHE, SVFDC, XARR,
     3             CB1, CB2, NCFM)
C
      RETURN
      END
C*********************************************************************
