C*********************************************************************
C subroutine Adapted Matrix Curl of Coriolis Force Add ***************
C            -       -      -       -        -     -   ***************
C Steve Gibbons Mon Nov 15 08:07:49 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Adds FAC* all the elements of the matrix A which, when             C
C multiplied by the vector VI adds the curl of the Coriolis Force    C
C to the appropriate elements of vector VO.                          C
C                                                                    C
C Multiplying VI by A after running AMCCFO should be equivalent to   C
C calling ASVCCF with VI.                                            C
C                                                                    C
C Now MTI defines what each scalar function in a solution vector     C
C represents.                                                        C
C                                                                    C
C         MTI( IH ) = 1 --> harmonic is poloidal velocity            C
C         MTI( IH ) = 2 --> harmonic is toroidal velocity            C
C         MTI( IH ) = 3 --> harmonic is temperature.                 C
C         MTI( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MTI( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C See the documentation to the subroutine AMCCFT for details of the  C
C different interactions. Alternatively, they are given in appendix  C
C B of my thesis ('Dynamo Models of the Earth's Magnetic Field')     C
C     (University of Leeds 1998)                                     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
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
C     MTI      : Array length ( * ) - atleast length NHI             C
C                  See above for key. (corresponds to input vec.)    C
C     MLI      : Array length ( * ) - atleast length NHI             C
C                  Sph. harm. degree, l.                             C
C     MMI      : Array length ( * ) - atleast length NHI             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MPI      : Array length ( * ) - atleast length NHI             C
C                  Pointer array to finite difference coefficients.  C
C                  MPI( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     Note that MTO, MLO and MMO correspond to                       C
C    exactly the above variables but for the output and not input    C
C    vectors. NRI must equal NRO and must both equal NR.             C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NDRVS     : Highest derivative stored in SVFDC.                C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     ILP       : First radial node to act on.                       C
C                  when contributing to the poloidal vorticity eqn.s C
C     IRP       : Last radial node to act on.                        C
C                  when contributing to the poloidal vorticity eqn.s C
C     ILT       : First radial node to act on.                       C
C                  when contributing to the toroidal vorticity eqn.s C
C     IRT       : Last radial node to act on.                        C
C                  when contributing to the toroidal vorticity eqn.s C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NTHP      : The number of theta points.                        C
C     NPHP      : The number of phi points.                          C
C                                                                    C
C     N1        : First dimension of matrix, A.                      C
C     N2        : Second dimension of matrix, A.                     C
C     IMF       : Matrix format flag. (See MATIND)                   C
C     KL        : Number of lower diagonals in matrix                C
C     KU        : Number of upper diagonals in matrix                C
C     KLE       : Number of additional lower diagonals in matrix     C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     MMAX      : Maximum sph. harmonic order, m.                    C
C     M0        : Minimum non-zero wavenumber.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 )                             C
C     FAC       : Multiplier of curl to be added.                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
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
C     QST       : Work array - dim. ( LH*(LH+2), 3)                  C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AMCCFO( NR, INARR, MTI, MLI, MMI, MPI, MTO, MLO, MMO,
     1          FAC, NBN, NDRVS, NDRVM, ILP, IRP, ILT, IRT, LH, NTHP,
     2          NPHP, XARR, NFDCM, A, N1, N2, IMF, KL, KU, KLE, NDCS,
     3          SVFDC, MMAX, GAUX, GAUW, PA, DPA, FTF1, FTF2, FTF3,
     4          VF, QST, M0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MTI( * ), MLI( * ), MMI( * ), MPI( * ),
     1        MTO( * ), MLO( * ), MMO( * ), NBN, NDRVS, NDRVM, ILP,
     2        IRP, ILT, IRT, LH, NTHP, NPHP, NFDCM, N1, N2, IMF, KL,
     3        KU, KLE, NDCS, MMAX, M0
      DOUBLE PRECISION FAC, XARR( NR ), A( N1, N2 ),
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
      DOUBLE PRECISION
     1                 QST( LH*(LH+2), 3),
     2                 VF( NPHP, NTHP, 3), FTF1( 2*NPHP ),
     3                 FTF2( 2*NPHP ), FTF3( 2*NPHP )
      DOUBLE PRECISION
     1                 GAUX( NTHP ), GAUW( NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IFORMF, NRR, NH, IHI, IHO, ILNR, IRNR,
     1        IHD, IPARS( 4 ), IS, LI, MI, ICSI, IQSTI,
     2        LO, MO, ICSO, IQSTO, IHARM, INDSHC
      DOUBLE PRECISION LOW, DPARS( 1 ), WORK( 3 ),
     1        QAB, SAB, TAB
      PARAMETER ( LOW = 1.0d-9 )
      EXTERNAL AMCCFT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Don't bother checking ILNR, IRNR, NDRVS or NBN
C as these errors will all be trapped by AMGLICA
C
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 ) THEN
         PRINT *,' Subroutine AMCCFO.'
         PRINT *,' IFORMF = ', IFORMF
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine AMCCFO.'
         PRINT *,' NR   = ', NR
         PRINT *,' NRR  = ', NRR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C     . early exit?
C
      IF ( ABS(FAC).LT.LOW ) RETURN
C     .
C     . Loop around in harmonics
C     .
      DO IHI = 1, NH
        IF ( MTI( IHI ).NE.1 .AND. MTI( IHI ).NE.2 ) GOTO 50
C
C Allocate whichever finite difference scheme
C belongs to the input harmonic
C
        IS = MPI( IHI )
C
C Consider separately the cases of input harmonic
C being poloidal and toroidal.
C
        IF ( MTI( IHI ).EQ.1 ) THEN
C         .
C         . Input harmonic is poloidal ...
C         .
          LI = MLI( IHI )
          IF ( MMI( IHI ).LT.0 ) THEN
            MI   = -MMI( IHI )
            ICSI = 2
          ELSE
            MI   = MMI( IHI )
            ICSI = 1
          ENDIF
C         .
C         . Call CORCOO with IQSTI = 1.
C         . This will give us K_{qs} and K_{qt} 
C         . (c.f. eq B.42 in thesis) - K_{qq} should be zero!
C         .
          IQSTI = 1
          CALL CORCOO( IQSTI, LI, MI, ICSI, QST, LH, NPHP, NTHP,
     1           FTF1, FTF2, FTF3, VF, GAUX, GAUW, PA, DPA, MMAX, M0 )
          IPARS( 1 ) = IQSTI
          IPARS( 2 ) = LI
C         .
C         . Now loop around the output harmonics
C         .
          DO IHO = 1, NH
C           .
            LO = MLO( IHO )
            IF ( MMO( IHO ).LT.0 ) THEN
              MO   = -MMO( IHO )
              ICSO = 2
            ELSE
              MO   = MMO( IHO )
              ICSO = 1
            ENDIF
            IHARM = INDSHC( LO, MO, ICSO )
            IF( IHARM.GT.0 ) SAB   = QST( IHARM, 2 )
            IF( IHARM.GT.0 ) TAB   = QST( IHARM, 3 )
C           .
C           . Case of poloidal harmonic
C           .
            IF ( MTO( IHO ).EQ.1 ) THEN
C              .
C              . Need to do K_{qt} interaction
C              .
               IHD        = 0
               ILNR       = ILP
               IRNR       = IRP
               IQSTO      = 3
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = TAB
C              .
               IF ( ABS( TAB ).GT.LOW ) THEN
C              .
               CALL AMLICA( N1, N2, KL, KU, KLE, IMF, IHI, IHO, INARR,
     1                   IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2                   NDRVS, NDRVM, IPARS, AMCCFT, A, FAC, XARR,
     3                   WORK, DPARS, SVFDC )
C              .
               ENDIF
C              .
            ENDIF
C           .
C           . Case of toroidal harmonic
C           .
            IF ( MTO( IHO ).EQ.2 ) THEN
C              .
C              . Need to do K_{qs} interaction
C              .
               IHD        = 1
               ILNR       = ILT
               IRNR       = IRT
               IQSTO      = 2
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = SAB
C              .
               IF ( ABS( SAB ).GT.LOW ) THEN
C              .
               CALL AMLICA( N1, N2, KL, KU, KLE, IMF, IHI, IHO, INARR,
     1                   IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2                   NDRVS, NDRVM, IPARS, AMCCFT, A, FAC, XARR,
     3                   WORK, DPARS, SVFDC )
C              .
               ENDIF
C              .
            ENDIF
C           .
          ENDDO
C         .
C         . Call CORCOO with IQSTI = 2.
C         . This will give us K_{sq}, K_{ss} and K_{st}
C         . (c.f. eq B.43 in thesis)
C         .
          IQSTI = 2
          CALL CORCOO( IQSTI, LI, MI, ICSI, QST, LH, NPHP, NTHP,
     1           FTF1, FTF2, FTF3, VF, GAUX, GAUW, PA, DPA, MMAX, M0 )
          IPARS( 1 ) = IQSTI
          IPARS( 2 ) = LI
C         .
C         . Now loop around the output harmonics
C         .
          DO IHO = 1, NH
C           .
            LO = MLO( IHO )
            IF ( MMO( IHO ).LT.0 ) THEN
              MO   = -MMO( IHO )
              ICSO = 2
            ELSE
              MO   = MMO( IHO )
              ICSO = 1
            ENDIF
            IHARM = INDSHC( LO, MO, ICSO )
            IF( IHARM.GT.0 ) QAB   = QST( IHARM, 1 )
            IF( IHARM.GT.0 ) SAB   = QST( IHARM, 2 )
            IF( IHARM.GT.0 ) TAB   = QST( IHARM, 3 )
C           .
C           . Case of poloidal harmonic
C           .
            IF ( MTO( IHO ).EQ.1 ) THEN
C              .
C              . Need to do K_{st} interaction
C              .
               IHD        = 1
               ILNR       = ILP
               IRNR       = IRP
               IQSTO      = 3
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = TAB
C              .
               IF ( ABS( TAB ).GT.LOW ) THEN
C              .
               CALL AMLICA( N1, N2, KL, KU, KLE, IMF, IHI, IHO, INARR,
     1                   IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2                   NDRVS, NDRVM, IPARS, AMCCFT, A, FAC, XARR,
     3                   WORK, DPARS, SVFDC )
C              .
               ENDIF
C              .
            ENDIF
C           .
C           . Case of toroidal harmonic
C           .
            IF ( MTO( IHO ).EQ.2 ) THEN
C              .
C              . Need to do K_{sq} interaction
C              .
               IHD        = 1
               ILNR       = ILT
               IRNR       = IRT
               IQSTO      = 1
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = QAB
C              .
               IF ( ABS( QAB ).GT.LOW ) THEN
C              .
               CALL AMLICA( N1, N2, KL, KU, KLE, IMF, IHI, IHO, INARR,
     1                   IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2                   NDRVS, NDRVM, IPARS, AMCCFT, A, FAC, XARR,
     3                   WORK, DPARS, SVFDC )
C              .
               ENDIF
C              .
C              . Need to do K_{ss} interaction
C              .
               IHD        = 2
               ILNR       = ILT
               IRNR       = IRT
               IQSTO      = 2
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = SAB
C              .
               IF ( ABS( SAB ).GT.LOW ) THEN
C              .
               CALL AMLICA( N1, N2, KL, KU, KLE, IMF, IHI, IHO, INARR,
     1                   IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2                   NDRVS, NDRVM, IPARS, AMCCFT, A, FAC, XARR,
     3                   WORK, DPARS, SVFDC )
C              .
               ENDIF
C              .
            ENDIF
C           .
          ENDDO
C         .
        ENDIF
C       endif case of ihi poloidal
C    
        IF ( MTI( IHI ).EQ.2 ) THEN
C         .
C         . Input harmonic is toroidal ...
C         .
          LI = MLI( IHI )
          IF ( MMI( IHI ).LT.0 ) THEN
            MI   = -MMI( IHI )
            ICSI = 2
          ELSE
            MI   = MMI( IHI )
            ICSI = 1
          ENDIF
C         .
C         . Call CORCOO with IQSTI = 3.
C         . This will give us K_{tq}, K_{ts} and K_{tt}
C         . (c.f. eq B.44 in thesis)
C         .
          IQSTI = 3
          CALL CORCOO( IQSTI, LI, MI, ICSI, QST, LH, NPHP, NTHP,
     1           FTF1, FTF2, FTF3, VF, GAUX, GAUW, PA, DPA, MMAX, M0 )
          IPARS( 1 ) = IQSTI
          IPARS( 2 ) = LI
C         .
C         . Now loop around the output harmonics
C         .
          DO IHO = 1, NH
C           .
            LO = MLO( IHO )
            IF ( MMO( IHO ).LT.0 ) THEN
              MO   = -MMO( IHO )
              ICSO = 2
            ELSE
              MO   = MMO( IHO )
              ICSO = 1
            ENDIF
            IHARM = INDSHC( LO, MO, ICSO )
            IF( IHARM.GT.0 ) QAB   = QST( IHARM, 1 )
            IF( IHARM.GT.0 ) SAB   = QST( IHARM, 2 )
            IF( IHARM.GT.0 ) TAB   = QST( IHARM, 3 )
C           .
C           . Case of poloidal harmonic
C           .
            IF ( MTO( IHO ).EQ.1 ) THEN
C              .
C              . Need to do K_{tt} interaction
C              .
               IHD        = 0
               ILNR       = ILP
               IRNR       = IRP
               IQSTO      = 3
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = TAB
C              .
               IF ( ABS( TAB ).GT.LOW ) THEN
C              .
               CALL AMLICA( N1, N2, KL, KU, KLE, IMF, IHI, IHO, INARR,
     1                   IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2                   NDRVS, NDRVM, IPARS, AMCCFT, A, FAC, XARR,
     3                   WORK, DPARS, SVFDC )
C              .
               ENDIF
C              .
            ENDIF
C           .
C           . Case of toroidal harmonic
C           .
            IF ( MTO( IHO ).EQ.2 ) THEN
C              .
C              . Need to do K_{tq} interaction
C              .
               IHD        = 0
               ILNR       = ILT
               IRNR       = IRT
               IQSTO      = 1
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = QAB
C              .
               IF ( ABS( QAB ).GT.LOW ) THEN
C              .
               CALL AMLICA( N1, N2, KL, KU, KLE, IMF, IHI, IHO, INARR,
     1                   IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2                   NDRVS, NDRVM, IPARS, AMCCFT, A, FAC, XARR,
     3                   WORK, DPARS, SVFDC )
C              .
               ENDIF
C              .
C              . Need to do K_{ts} interaction
C              .
               IHD        = 1
               ILNR       = ILT
               IRNR       = IRT
               IQSTO      = 2
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = SAB
C              .
               IF ( ABS( SAB ).GT.LOW ) THEN
C              .
               CALL AMLICA( N1, N2, KL, KU, KLE, IMF, IHI, IHO, INARR,
     1                   IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2                   NDRVS, NDRVM, IPARS, AMCCFT, A, FAC, XARR,
     3                   WORK, DPARS, SVFDC )
C              .
               ENDIF
C              .
            ENDIF
C           .
          ENDDO
C         .
        ENDIF
C       endif case of ihi toroidal
C    
 50   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
