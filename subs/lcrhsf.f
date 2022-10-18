C*********************************************************************
C subroutine Landau Couette Right Hand Side Form *********************
C            -      -       -     -    -    -    *********************
C Steve Gibbons Sun May  7 18:34:46 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Let a_r be the real function (harmonic 1) and let a_i be the       C
C imaginary function (harmonic 2 ). Let ^{i} denote current time     C
C step and ^{i+1} the next time step.                                C
C                                                                    C
C LCRHSF adds to the real part of the RHS vector:                    C
C                                                                    C
C [ 1 + deltat.c.(lambda-x^2) + deltat.c.d^2/dx^2 ] a_r^{i}          C
C - 2.deltat.c.kappa.x a_i^{i}                                       C
C - c.deltat.[ a_r^{i}a_r^{i} + a_i^{i}a_i^{i} ] a_r^{i}             C
C -(1-c).deltat.[ a_r^{i+1}a_r^{i+1} + a_i^{i+1}a_i^{i+1}] a_r^{i+1} C
C                                                                    C
C LCRHSF adds to the imag part of the RHS vector:                    C
C                                                                    C
C [ 1 + deltat.c.(lambda-x^2) + deltat.c.d^2/dx^2 ] a_i^{i}          C
C + 2.deltat.c.kappa.x a_r^{i}                                       C
C - c.deltat.[ a_r^{i}a_r^{i} + a_i^{i}a_i^{i} ] a_i^{i}             C
C -(1-c).deltat.[ a_r^{i+1}a_r^{i+1} + a_i^{i+1}a_i^{i+1}] a_i^{i+1} C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VI        : Dim (*). Vector at time step i.                    C
C     VIP1      : Dim (*). Vector at time step i+1.                  C
C     RHS       : Dim (*). Right hand vector.                        C
C                                                                    C
C     DPARS     : Array dim (*).                                     C
C                 DPARS( 1 ) = LAMBDA                                C
C                 DPARS( 2 ) = KAPPA                                 C
C                                                                    C
C     DELTAT    : Time step.                                         C
C     CFAC      : Implicit/Explicit ratio.                           C
C                 CFAC must be greater than 0 and less than 1.       C
C                 If cfac was 1, this would be implicit.             C
C                 If cfac was 0, this would be explicit.             C
C                 cfac = 0.5 --> Crank-Nicolson.                     C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( 2*NBN+1, NR, 3, 1).                   C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LCRHSF( NR, NBN, VI, VIP1, RHS, DPARS,
     1                   DELTAT, CFAC, XARR, SVFDC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NBN
      DOUBLE PRECISION VI( * ), VIP1( * ), RHS( * ), DPARS( * ),
     1                 SVFDC( 2*NBN+1, NR, 3, 1 ), DELTAT, CFAC,
     2                 XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER ILNR, IRNR, IFORMF, INARR( 3 ), NDCS, IHD, IOP, NCFM,
     1        NDRVS, IS, N2, IR, INDFUN, INDRF, INDIF, IH
      DOUBLE PRECISION ZERO, FAC, DERV( 3 ), ARTSI, AR2TSI, X,
     1                 ARTSI1, AITSI, AI2TSI, AITSI1, LAMBDA, KAPPA
      PARAMETER ( ZERO = 0.0d0, IFORMF = 3, NDRVS = 2 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( CFAC.LE.0.0d0 .OR. CFAC.GE.1.0d0 ) THEN
        PRINT *,' Subroutine LCRHSF.'
        PRINT *,' CFAC = ', CFAC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      LAMBDA = DPARS( 1 )
      KAPPA  = DPARS( 2 )
C
      IS   = 1
      NCFM = 2*NBN+1
      NDCS = 1
      IHD  = 2
C
      INARR( 1 ) = IFORMF
      INARR( 2 ) = NR
      INARR( 3 ) = 2
C
      N2 = 2*NR
C
C Zero the RHS vector
C
      IOP = 0
      CALL VECOP( RHS, ZERO, N2, IOP )
C
      ILNR       = 2
      IRNR       = NR - 1
C
      DO IR = ILNR, IRNR
        X = XARR( IR )
C       .
C       . Calculate derivatives of current time-step vector
C       . First: real part
C       .
        IH    = 1
        INDRF = INDFUN( IR, IH, INARR )
        CALL ASVDR ( VI, IR, IS, IH, NBN, IHD, NCFM, NR, NDRVS,
     1               NDRVS, DERV, INARR, SVFDC, NDCS )
        ARTSI  = DERV( 1 )
        AR2TSI = DERV( 3 )
        ARTSI1 = VIP1( INDRF )
C       .
C       . artsi contains a_r^{i}
C       . ar2tsi contains d^2 a_r^{i} / dx^2
C       . artsi1 contains a_r^{i+1}
C       . Now - imaginary part.
C       .
        IH    = 2
        INDIF = INDFUN( IR, IH, INARR )
        CALL ASVDR ( VI, IR, IS, IH, NBN, IHD, NCFM, NR, NDRVS,
     1               NDRVS, DERV, INARR, SVFDC, NDCS )
        AITSI  = DERV( 1 )
        AI2TSI = DERV( 3 )
        AITSI1 = VIP1( INDIF )
C       .
C       . aitsi contains a_i^{i}
C       . ai2tsi contains d^2 a_i^{i} / dx^2
C       . aitsi1 contains a_i^{i+1}
C       . We can now evaluate the term at this node.
C       .
        FAC = 1.0d0 + DELTAT*CFAC*(LAMBDA - X*X)
        RHS( INDRF ) = RHS( INDRF ) + FAC*ARTSI
        RHS( INDIF ) = RHS( INDIF ) + FAC*AITSI
C
        RHS( INDRF ) = RHS( INDRF ) + DELTAT*CFAC*AR2TSI
        RHS( INDIF ) = RHS( INDIF ) + DELTAT*CFAC*AI2TSI
C
        FAC = 2.0d0*DELTAT*CFAC*KAPPA*X
        RHS( INDRF ) = RHS( INDRF ) - FAC*AITSI
        RHS( INDIF ) = RHS( INDIF ) + FAC*ARTSI
C
        FAC = (ARTSI*ARTSI + AITSI*AITSI)*CFAC*DELTAT
        RHS( INDRF ) = RHS( INDRF ) - FAC*ARTSI
        RHS( INDIF ) = RHS( INDIF ) - FAC*AITSI
        FAC = (ARTSI1*ARTSI1 + AITSI1*AITSI1)*(1.0d0-CFAC)*DELTAT
        RHS( INDRF ) = RHS( INDRF ) - FAC*ARTSI1
        RHS( INDIF ) = RHS( INDIF ) - FAC*AITSI1
C       .
      ENDDO
C
      RETURN
      END
C*********************************************************************
