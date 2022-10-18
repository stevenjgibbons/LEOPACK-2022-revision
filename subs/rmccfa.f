C*********************************************************************
C subroutine Recalled coefficient Matrix Curl of Coriolis Force Add **
C            -                    -      -       -        -     -   **
C Steve Gibbons Tue Jan 18 08:01:30 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Adds FAC* all the elements of the matrix A which, when             C
C multiplied by the vector VI adds the curl of the Coriolis Force    C
C to the appropriate elements of vector VO.                          C
C                                                                    C
C Multiplying VI by A after running RMCCFA should be equivalent to   C
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
C Unlike AMCCFA which evaluates the vector interactions in situ,     C
C RMCCFA searches through the arrays IHNALP, IHNGAM and TVHI (as     C
C pre-calculated by CFVICC) with the help of VLICER for the non-zero C
C interactions.                                                      C
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
C     NVI      : Number of vector interactions.                      C
C     IHNALP   : Number of alpha harmonics. Dim ( * )                C
C     IHNGAM   : Number of gamma harmonics. Dim ( * )                C
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
      SUBROUTINE RMCCFA( NR, INARR, MTI, MLI, MMI, MPI, MTO, MLO, MMO,
     1          FAC, NBN, NDRVS, NDRVM, ILP, IRP, ILT, IRT, XARR,
     2          NFDCM, A, N1, N2, IMF, KL, KU, KLE, NDCS, SVFDC, 
     3          NVI, IHNALP, IHNGAM, TVHI, CVI )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MTI( * ), MLI( * ), MMI( * ), MPI( * ),
     1        MTO( * ), MLO( * ), MMO( * ), NBN, NDRVS, NDRVM, ILP,
     2        IRP, ILT, IRT, NFDCM, N1, N2, IMF, KL, KU, KLE, NDCS,
     3        NVI, IHNALP( * ), IHNGAM( * )
      CHARACTER *(4) TVHI( * )
      DOUBLE PRECISION FAC, XARR( NR ), A( N1, N2 ),
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS ), CVI( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IFORMF, NRR, NH, IHI, IHO, ILNR, IRNR,
     1        IHD, IPARS( 4 ), IS, LI, MI, ICSI, IQSTI,
     2        LO, MO, ICSO, IQSTO, INDO, INDSHC, INDI
      DOUBLE PRECISION LOW, DPARS( 1 ), WORK( 2 ),
     1        QAB, SAB, TAB
      PARAMETER ( LOW = 1.0d-9 )
      EXTERNAL AMCCFT
      LOGICAL ONZIC
      CHARACTER *(4) REQVIT
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
         PRINT *,' Subroutine RMCCFA.'
         PRINT *,' IFORMF = ', IFORMF
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine RMCCFA.'
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
          INDI   = INDSHC( LI, MI, ICSI )
C         .
C         . This will give us K_{qs} and K_{qt} 
C         . (c.f. eq B.42 in thesis) - K_{qq} should be zero!
C         .
          IQSTI = 1
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
            INDO = INDSHC( LO, MO, ICSO )
C           .
C           . Case of poloidal harmonic
C           .
            IF ( MTO( IHO ).EQ.1 ) THEN
C              .
C              . Need to do K_{qt} interaction
C              .
               REQVIT = 'KCQT'
               CALL VLICER( NVI, IHNALP, IHNGAM, INDI, INDO, TVHI,
     1                      REQVIT, CVI, TAB, ONZIC )
C              .
               IF ( ONZIC ) THEN
C              .
               IHD        = 0
               ILNR       = ILP
               IRNR       = IRP
               IQSTO      = 3
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = TAB
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
               REQVIT = 'KCQS'
               CALL VLICER( NVI, IHNALP, IHNGAM, INDI, INDO, TVHI,
     1                      REQVIT, CVI, SAB, ONZIC )
C              .
               IF ( ONZIC ) THEN
C              .
               IHD        = 1
               ILNR       = ILT
               IRNR       = IRT
               IQSTO      = 2
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = SAB
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
C         . This will give us K_{sq}, K_{ss} and K_{st}
C         . (c.f. eq B.43 in thesis)
C         .
          IQSTI = 2
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
            INDO = INDSHC( LO, MO, ICSO )
C           .
C           . Case of poloidal harmonic
C           .
            IF ( MTO( IHO ).EQ.1 ) THEN
C              .
C              . Need to do K_{st} interaction
C              .
               REQVIT = 'KCST'
               CALL VLICER( NVI, IHNALP, IHNGAM, INDI, INDO, TVHI,
     1                      REQVIT, CVI, TAB, ONZIC )
C              .
               IF ( ONZIC ) THEN
C              .
               IHD        = 1
               ILNR       = ILP
               IRNR       = IRP
               IQSTO      = 3
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = TAB
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
               REQVIT = 'KCSQ'
               CALL VLICER( NVI, IHNALP, IHNGAM, INDI, INDO, TVHI,
     1                      REQVIT, CVI, QAB, ONZIC )
C              .
               IF ( ONZIC ) THEN
C              .
               IHD        = 1
               ILNR       = ILT
               IRNR       = IRT
               IQSTO      = 1
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = QAB
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
               REQVIT = 'KCSS'
               CALL VLICER( NVI, IHNALP, IHNGAM, INDI, INDO, TVHI,
     1                      REQVIT, CVI, SAB, ONZIC )
C              .
               IF ( ONZIC ) THEN
C              .
               IHD        = 2
               ILNR       = ILT
               IRNR       = IRT
               IQSTO      = 2
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = SAB
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
          INDI   = INDSHC( LI, MI, ICSI )
C         .
C         . This will give us K_{tq}, K_{ts} and K_{tt}
C         . (c.f. eq B.44 in thesis)
C         .
          IQSTI = 3
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
            INDO = INDSHC( LO, MO, ICSO )
C           .
C           . Case of poloidal harmonic
C           .
            IF ( MTO( IHO ).EQ.1 ) THEN
C              .
C              . Need to do K_{tt} interaction
C              .
               REQVIT = 'KCTT'
               CALL VLICER( NVI, IHNALP, IHNGAM, INDI, INDO, TVHI,
     1                      REQVIT, CVI, TAB, ONZIC )
C              .
               IF ( ONZIC ) THEN
C              .
               IHD        = 0
               ILNR       = ILP
               IRNR       = IRP
               IQSTO      = 3
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = TAB
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
               REQVIT = 'KCTQ'
               CALL VLICER( NVI, IHNALP, IHNGAM, INDI, INDO, TVHI,
     1                      REQVIT, CVI, QAB, ONZIC )
C              .
               IF ( ONZIC ) THEN
C              .
               IHD        = 0
               ILNR       = ILT
               IRNR       = IRT
               IQSTO      = 1
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = QAB
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
               REQVIT = 'KCTS'
               CALL VLICER( NVI, IHNALP, IHNGAM, INDI, INDO, TVHI,
     1                      REQVIT, CVI, SAB, ONZIC )
C              .
               IF ( ONZIC ) THEN
C              .
               IHD        = 1
               ILNR       = ILT
               IRNR       = IRT
               IQSTO      = 2
               IPARS( 3 ) = IQSTO
               IPARS( 4 ) = LO
               DPARS( 1 ) = SAB
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
