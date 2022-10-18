C*********************************************************************
C subroutine Newton Raphson Convection Woodbury Matrices Form ********
C            -      -       -          -        -        -    ********
C Steve Gibbons Wed Mar 15 14:52:18 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C When solving the vorticity and heat equations using a Newton-      C
C Raphson type method: either for boundary forced stationary flows   C
C or for solutions steady in a drifitng frame, we will need to       C
C provide a column for solving for the drift rate (if applicable -   C
C this is when OTCDC = .TRUE.) and we will need to change rows to    C
C prevent arbitrary rotations.                                       C
C                                                                    C
C The t_1^0 harmonic needs such a relation if boundaries are         C
C stress free. In this case, OT10 is .TRUE. and IT10 contains the    C
C number of the T_1^0 harmonic.                                      C
C                                                                    C
C The t_1^{1c} and t_1^{1s} harmonics needs such a relation if       C
C boundaries are stress free and the Coriolis force is absent from   C
C the matrix. In this case, OT11C and OT11S are .TRUE. and IT11C and C
C IT11S respectively contain the numbers of the T_1^{1c} and         C
C T_1^{1s} harmonics.                                                C
C                                                                    C
C There are a maximum number of 4 such operations. The actual        C
C number is given by the number NTS.                                 C
C                                                                    C
C The appropriate terms are then dealt with in the matrix A, the     C
C right hand vector, RHS, and the matrices U and V for use in        C
C the routine BMWDFS.                                                C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. Must be equal to KL in this routine.      C
C                                                                    C
C     NTS       : Number of special row/columns (see above).         C
C     ITCDC     : Number of harmonic to be fixed in rotating frame.  C
C     IT10      : Number of T_1^{0c} harmonic.                       C
C     IT11C     : Number of T_1^{1c} harmonic.                       C
C     IT11S     : Number of T_1^{1s} harmonic.                       C
C                                                                    C
C Note that ITCDC, IT10, IT11C and IT11S are only referred to if     C
C OTCDC, OT10, OT11C and OT11S are respectively .TRUE.               C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARR     : Dim( * ). Format info. (See INDFUN)                C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dimensions ( N1, N2 ) format as above      C
C     RHS       : Right hand vector. Dim (N2).                       C
C     U         : Matrix dim ( N2, NTS ). See above for guidance.    C
C     V         : Matrix dim ( N2, NTS ). See above for guidance.    C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     DFTDV     : Drifting frame time derivative vector. Dim (N2).   C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     OTCDC     : .T. if and only if drift rate column to be formed. C
C     OT10      : .T. if and only if T_1^0 needs frame fixing.       C
C     OT11C     : .T. if and only if T_1^{1c} needs frame fixing.    C
C     OT11S     : .T. if and only if T_1^{1s} needs frame fixing.    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NRCWMF( N1, N2, KL, KU, KLE, NTS, ITCDC, IT10, IT11C,
     1                   IT11S, NR, INARR, OTCDC, OT10, OT11C, OT11S,
     2                   A, RHS, U, V, XARR, DFTDV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, NTS, ITCDC, IT10, IT11C, IT11S,
     1        NR, INARR( * )
      LOGICAL OTCDC, OT10, OT11C, OT11S
      DOUBLE PRECISION A( N1, N2 ), RHS( N2 ), U( N2, NTS ),
     1                 V( N2, NTS ), XARR( NR ), DFTDV( N2 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NNDS, NTS2, IOP, IRR, INDFUN, IND, IR, IRC, IGO, IH,
     1        INODE
      DOUBLE PRECISION ZERO, HI, DLOW, QUOT, X0
      PARAMETER ( NNDS = 6, ZERO = 0.0d0, DLOW = 1.0d-7 )
      INTEGER IWORK( NNDS )
      DOUBLE PRECISION WORK1( NNDS ), WORK2( NNDS ),
     1                 COEFM( NNDS, NNDS )
      LOGICAL ODO
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C 
      NTS2 = 0
      IF ( OTCDC ) NTS2 = NTS2 + 1
      IF ( OT10  ) NTS2 = NTS2 + 1
      IF ( OT11C ) NTS2 = NTS2 + 1
      IF ( OT11S ) NTS2 = NTS2 + 1
      IF ( NTS2.NE.NTS ) THEN
        PRINT *,' Subroutine NRCWMF.'
        PRINT *,' NTS = ', NTS,' but logical var.s suggest ', NTS2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Early escape??
C
      IF ( NTS.EQ.0 ) RETURN
C
C Clear the arrays U and V
C
      IOP = 0
      CALL MATOP( U, ZERO, N2, NTS, IOP )
      CALL MATOP( V, ZERO, N2, NTS, IOP )
C
      NTS2 = 0
      IRR  = NR/2
C
C Do time derivative column.
C This is a column operation so put vector in U
C matrix and location in V
C
      IF ( OTCDC ) THEN
        NTS2 = NTS2 + 1
        IND  = INDFUN( IRR, ITCDC, INARR )
        DO IR = 1, N2
          U( IR, NTS2 ) = DFTDV( IR )
        ENDDO
        U( IND, NTS2 ) = U( IND, NTS2 ) - 1.0d0
        V( IND, NTS2 ) = 1.0d0
        IRC = 2
        CALL BMRCOP( KL, KU, KLE, N2, IRC, IND, A, ZERO, ZERO )
        A( KLE + KU + 1, IND ) = 1.0d0
      ENDIF
C
C Sort out toroidal solid body rotations
C
      DO IGO = 1, 3
        ODO  = .FALSE.
        IF ( IGO.EQ.1 .AND. OT10  ) THEN
          ODO = .TRUE.
          IH  = IT10
        ENDIF
        IF ( IGO.EQ.2 .AND. OT11C ) THEN
          ODO = .TRUE.
          IH  = IT11C
        ENDIF
        IF ( IGO.EQ.3 .AND. OT11S ) THEN
          ODO = .TRUE.
          IH  = IT11S
        ENDIF
        IF ( ODO ) THEN
          NTS2 = NTS2 + 1
          IND  = INDFUN( IRR, IH, INARR )
          IRC = 1
          CALL BMRCOP( KL, KU, KLE, N2, IRC, IND, A, ZERO, ZERO )
          A( KLE + KU + 1, IND ) = 1.0d0
          U( IND, NTS2 )         = 1.0d0
          V( IND, NTS2 )         = (-1.0d0)
          RHS( IND )             = ZERO   
C         .
C         . Loop IR from 2, NR - 1 and fill in
C         . the 'obvious' values.
C         .
          DO IR = 2, NR - 1
            IND = INDFUN( IR, IH, INARR )
            HI = XARR( IR + 1 ) - XARR( IR - 1 )
            V( IND, NTS2 ) = V( IND, NTS2 ) + 0.5d0*HI/XARR( IR )
          ENDDO
C         .
C         . Now add components from inner boundary
C         .
          IF ( XARR( 1 ).GT.DLOW ) THEN
            X0 = XARR( 1 )
            DO INODE = 1, NNDS
              WORK1( INODE ) = XARR( INODE )
            ENDDO
            CALL GFDCFD( X0, WORK1, NNDS, COEFM, NNDS, IWORK, WORK2 )
C           .
C The coeff of f( x_i ) in the derivative of f at xarr( 1 )
C is now stored in COEFM( 2, i )
C           .
            QUOT = 1.0d0 - X0*COEFM( 2, 1 )
            IF ( ABS( QUOT ).LT.DLOW ) THEN
              PRINT *,' Subroutine NRCWMF.'
              PRINT *,' QUOT = ', QUOT
              PRINT *,' Program aborted.'
              STOP
            ENDIF
C           .
            HI = XARR( 2 ) - XARR( 1 )
            DO IR = 2, NNDS
              IND = INDFUN( IR, IH, INARR )
              V( IND, NTS2 ) = V( IND, NTS2 ) +
     1                    0.5d0*HI*COEFM( 2, IR )/(QUOT*X0)
            ENDDO
C           .
          ENDIF
C         .
C Now add components from outer boundary
C         .
          X0 = XARR( NR )
          HI = X0 - XARR( NR - 1 )
          DO INODE = 1, NNDS
            WORK1( INODE ) = XARR( NR + 1 - INODE )
          ENDDO
          CALL GFDCFD( X0, WORK1, NNDS, COEFM, NNDS, IWORK, WORK2 )
C         .
C The coeff of f( x_i ) in the derivative of f at xarr( nr )
C is now stored in COEFM( 2, nr + 1 - i )
C         .
          QUOT = 1.0d0 - X0*COEFM( 2, 1 )
          IF ( ABS( QUOT ).LT.DLOW ) THEN
            PRINT *,' Subroutine NRCWMF.'
            PRINT *,' QUOT = ', QUOT
            PRINT *,' Program aborted.'
            STOP
          ENDIF
C         .
          DO INODE = 2, NNDS
            IR = NR + 1 - INODE
            IND = INDFUN( IR, IH, INARR )
            V( IND, NTS2 ) = V( IND, NTS2 ) +
     1                  0.5d0*HI*COEFM( 2, INODE )/(QUOT*X0)
          ENDDO
C         .
        ENDIF
      ENDDO
C
      RETURN
      END
C*********************************************************************
