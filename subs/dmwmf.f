C*********************************************************************
C subroutine Diffusion Matrix Woodbury Matrices Form *****************
C            -         -      -        -        -    *****************
C Steve Gibbons Sat Mar  4 13:55:46 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C If a diffusion matrix, DMAT, is formed in advance by SFDDMF then   C
C this routine forms the matrices U and V which are used by the      C
C routine BMWDFS to solve (using the Woodbury formula) the linear    C
C system. The purpose of this procedure is to avoid a degeneracy     C
C to solid body rotations when stress free boundary conditions are   C
C employed.                                                          C
C                                                                    C
C In addition, the appropriate rows of RHS are put to zero.          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N2        : Second dimension of matrix DMAT.                   C
C                 Must be NH*NR                                      C
C     NR        : Total number of radial grid nodes.                 C
C                                                                    C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NRR. Must be consistent with NR.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). See below.                     C
C     MHOBC     : Dimension ( NDCS ). See below.                     C
C                                                                    C
C  MHIBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( ih ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                        where L = LARR( ih )                        C
C                                                                    C
C  Similarly, at the outer boundary:-                                C
C                                                                    C
C  MHOBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( ih ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                        where L = LARR( ih )                        C
C                                                                    C
C     NNDS      : Number of nodes allowed for derivatives.           C
C                 (Usually 2*NBN + 1 - but quite arbitrary provided  C
C                  that it is atleast 2).                            C
C                                                                    C
C     IWORK     : Work array - length ( NNDS ).                      C
C                                                                    C
C     NTS       : Number of toroidal singularities.                  C
C                 NTS = 0 --> no rows missing.                       C
C                 NTS = 1 --> row corresponding to t_1^0( NR/2 )     C
C                             must be replaced with a vector in V.   C
C                 NTS = 3 --> row corresponding to t_1^0( NR/2 )     C
C                             t_1^{1c}( NR/2 ) and t_1^{1s}( NR/2 )  C
C                             must be replaced with a vector in V.   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     U         : Dimension ( N2, NTS ).                             C
C                 Matrix indicating location of rows in matrix       C
C     V         : Dimension ( N2, NTS ).                             C
C                 Matrix containing additional rows in matrix        C
C                   (See BMWDFS for clearer guidance on U and V).    C
C                                                                    C
C     WORK1     : Work array dim ( NNDS )                            C
C     WORK2     : Work array dim ( NNDS )                            C
C     COEFM     : Work array dim ( NNDS, NNDS )                      C
C                                                                    C
C     RHS       : Right hand side vector. Dim ( NR*NH )              C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DMWMF( N2, NR, INARR, MHT, MHL, MHP, MHIBC, MHOBC,
     1                  NNDS, IWORK, NTS, XARR, U, V, WORK1,
     2                  WORK2, COEFM, RHS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N2, NR, INARR( * ), MHT( * ), MHL( * ), MHP( * ),
     1        MHIBC( * ), MHOBC( * ), NNDS, IWORK( NNDS ), NTS
      DOUBLE PRECISION XARR( NR ), U( N2, NTS ), V( N2, NTS ),
     1                 WORK1( NNDS ), WORK2( NNDS ), RHS( * ),
     2                 COEFM( NNDS, NNDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, IND, INDFUN, IH, NH, IS, NTS2, IR, IRAD, INODE
      DOUBLE PRECISION ZERO, ONE, HI, QUOT, X0, LOW
      PARAMETER ( ZERO = 0.0d0, ONE = 1.0d0, LOW = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Early escape??
C
      IF ( NTS.EQ.0 ) RETURN
C
      NH = INARR( 3 )
C
C Zero matrices U and V
C
      IOP = 0
      CALL MATOP( U, ZERO, N2, NTS, IOP )
      CALL MATOP( V, ZERO, N2, NTS, IOP )
C
C Now search for lines to replace
C
      NTS2 = 0
      IRAD = NR/2
      DO IH = 1, NH
        IS  = MHP( IH )
        IF ( MHT( IH ).EQ.2   .AND. MHL( IH ).EQ.1   .AND.
     1       MHIBC( IS ).EQ.6 .AND. MHOBC( IS ).EQ.6 ) THEN
          NTS2 = NTS2 + 1
          IND  = INDFUN( IRAD, IH, INARR )
          U( IND, NTS2 ) = ONE
          V( IND, NTS2 ) = (-1.0d0)
          RHS( IND ) = ZERO
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
          IF ( XARR( 1 ).GT.LOW ) THEN
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
            IF ( ABS( QUOT ).LT.LOW ) THEN
              PRINT *,' Subroutine DMWMF.'
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
          IF ( ABS( QUOT ).LT.LOW ) THEN
            PRINT *,' Subroutine DMWMF.'
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
      IF ( NTS.NE.NTS2 ) THEN
        PRINT *,' NTS = ', NTS,' NTS2 = ', NTS2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
