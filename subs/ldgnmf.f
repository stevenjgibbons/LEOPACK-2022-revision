C*********************************************************************
C subroutine Linear Dependence of Grid Node Matrix Form **************
C            -      -             -    -    -      -    **************
C Steve Gibbons Sat Oct 23 15:01:52 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Let f be a function of x and let f_j denote the value of f at      C
C the grid node j (x value is x_j given by XARR( j ) ... )           C
C                                                                    C
C If f has to satisfy a particular boundary condition then           C
C all the f_j may not be linearly independent.                       C
C                                                                    C
C For a group of n ( = NNDS ) grid nodes, from s = k+1 to k+n        C
C for some integer, k,                                               C
C                                                                    C
C  f_{k+i} = \sum_{s=1}^n DMAT( i, s ) f_{k+s},                      C
C                                                                    C
C and the routine LDGNMF returns the matrix DMAT.                    C
C                                                                    C
C The number of linearly dependent nodes to the left and right are   C
C given by NALF and NARF respectively.                               C
C                                                                    C
C The boundary conditions at the inner and outer boundaries are      C
C specified by the integers IIBC and IOBC which may take the         C
C following values                                                   C
C                                                                    C
C    IIBC              Inner Boundary Condition                      C
C    ====              ========================                      C
C                                                                    C
C      1       No condition to be applied                            C
C      2       Function must vanish                                  C
C      3       First derivative must vanish                          C
C      4       Function AND first derivative must vanish             C
C      5       Function AND second derivative must vanish            C
C      6       rdf/dr - f(r) = 0                                     C
C      7       r df/dr - l f(r) = 0                                  C
C                                                                    C
C    IOBC              Outer Boundary Condition                      C
C    ====              ========================                      C
C                                                                    C
C      1       No condition to be applied                            C
C      2       Function must vanish                                  C
C      3       First derivative must vanish                          C
C      4       Function AND first derivative must vanish             C
C      5       Function AND second derivative must vanish            C
C      6       rdf/dr - f(r) = 0                                     C
C      7       r df/dr + (l+1) f(r) = 0                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Total number of radial grid nodes.                 C
C     NNDS      : Number of nodes needed to take derivative.         C
C     NALF      : Number of nodes to left which are a linear         C
C                 combination of the other nodes.                    C
C     NARF      : Number of nodes to right which are a linear        C
C                 combination of the other nodes.                    C
C     L         : Spherical harmonic degree, l.                      C
C     IIBC      : Inner boundary flag - see above.                   C
C     IOBC      : Outer boundary flag - see above.                   C
C     NCFM      : Leading dimension of array DMAT etc.               C
C     IPCM      : Dimension ( NCFM ). Working array.                 C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension ( NR ).                         C
C                  XARR( i ) contains the value of r at the          C
C                   i^{th} grid node.                                C
C                                                                    C
C     DMAT      : Dimension ( NCFM, NCFM ). See above.               C
C     WMAT      : Dimension ( NCFM, NCFM ). Working array.           C
C     WORK1     : Dimension ( NCFM ). Working array.                 C
C     WORK2     : Dimension ( NCFM ). Working array.                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LDGNMF( NR, NNDS, NALF, NARF, L, IIBC, IOBC, NCFM,
     1                   XARR, DMAT, WMAT, WORK1, WORK2, IPCM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NNDS, NALF, NARF, L, IIBC, IOBC, NCFM,
     1        IPCM( NCFM )
      DOUBLE PRECISION XARR( NR ), DMAT( NCFM, NCFM),
     1                 WMAT( NCFM, NCFM), WORK1( NCFM ),
     2                 WORK2( NCFM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER I, IFN, ILN, NDS2, ID1, IHND
      DOUBLE PRECISION ZERO, TOL, X0, FAC, QUOT
      PARAMETER ( ZERO = 0.0d0, TOL = 1.0d-9 )
C ifn is first linearly independent node
C iln is last linearly independent node
C ihnd is the highest number derivative which
C will be needed to be calculated during this
C routine (by GFDCFD)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First check on values of NALF and NARF
C
      I = NALF*NARF
      IF ( I.NE.0 .OR. (NALF.LT.0) .OR. (NALF.GT.2) .OR.
     1                 (NARF.LT.0) .OR. (NARF.GT.2)       ) THEN
        PRINT *,' Subroutine LDGNMF.'
        PRINT *,' NALF = ', NALF
        PRINT *,' NARF = ', NARF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Ok - the number of affected nodes is valid
C     . so let's proceed. First, zero the matrix DMAT
C     .
      I = 0
      CALL MATOP( DMAT, ZERO, NCFM, NCFM, I )
C     .
C     . Now add the diagonal elements for the
C     . linearly independent variables
C     .
      IF ( NALF.EQ.0 ) IFN = 1
      IF ( NALF.EQ.1 ) IFN = 2
      IF ( NALF.EQ.2 ) IFN = 3
C     .
      IF ( NARF.EQ.0 ) ILN = NNDS
      IF ( NARF.EQ.1 ) ILN = NNDS - 1
      IF ( NARF.EQ.2 ) ILN = NNDS - 2
C     .
      DO I = IFN, ILN
        DMAT( I, I ) = 1.0d0
      ENDDO
C     .
C     . We can now return if the identity matrix is required
C     .
      IF ( NALF.EQ.0 .AND. NARF.EQ.0 ) RETURN
C     .
C     .
C     . OK - now decide whether we are doing inner or
C     . outer boundary
C     .
      IF ( NARF.EQ.0 ) THEN
C       .
C       . We are considering the inner boundary
C       . First, return if IIBC = 1 or 2
C       .
        IF ( IIBC.EQ.1 .OR. IIBC.EQ.2 ) RETURN
C       .
C       . So we need to calculate deriv.s at
C       . the inner boundary, XARR( 1 ).
C       . NDS2 is the number of nodes we can use at
C       . the inner boundary.
C       .
        IF ( IIBC.EQ.5 ) THEN
          IHND = 2
        ELSE
          IHND = 1
        ENDIF
C       .
        NDS2 = 1 + NNDS - NALF
C       .
        IF ( (NDS2-1).LT.IHND ) THEN
           PRINT *,' Subroutine LDGNMF '
           PRINT *,' NDS2                = ', NDS2
           PRINT *,' Required derivative = ', IHND
           PRINT *,' Program aborted.'
           STOP
        ENDIF
C       .
        DO I = 1, NDS2
          WORK1( I ) = XARR( I )
        ENDDO
        X0 = XARR( 1 )
        CALL GFDCFD ( X0, WORK1, NDS2, WMAT, NCFM, IPCM, WORK2 )
C       .
C       . WMAT( m + 1, i ) now contains the coefficient
C       . with which you multiply f_i to get the m^{th}
C       . derivative of f at x_0.
C       .
C       . Consider the case IIBC = 3
C       . We require the first derivative = 0
C       .
        IF ( IIBC.EQ.3 ) THEN
C          .
C          . Just quickly check that NALF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NALF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IIBC = ',IIBC,' and NALF = ',NALF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . Check for imminent division by zero
C          .
           IF ( ABS( WMAT( 2, 1 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT( 2, 1 ) = ',WMAT( 2, 1 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 2, NDS2
             DMAT( 1, I ) = (-1.0d0)*WMAT( 2, I )/WMAT( 2, 1 )
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . Consider the case IIBC = 4(5)
C       . We need both the function and the first
C       . (second) derivative to vanish. We do not need to
C       . to change anything to make the function
C       . vanish so just need do the first (second) deriv. cond.
C       .
        IF ( IIBC.EQ.4 .OR. IIBC.EQ.5 ) THEN
C          .
C          . Check for imminent division by zero
C          .
           IF ( IIBC.EQ.4 ) ID1 = 2
           IF ( IIBC.EQ.5 ) ID1 = 3
C          .
           IF ( ABS( WMAT( ID1, 2 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT(',ID1,', 2 ) = ',WMAT( ID1, 2 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 3, NDS2
             DMAT( NALF, NALF - 2 + I ) =
     1                    (-1.0d0)*WMAT( ID1, I )/WMAT( ID1, 2 )
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . We may treat the cases IIBC.EQ.6 and IIBC.EQ.7
C       . together as they are both conditions of the form
C       . rdf/dr + FAC f(r) = 0
C       .
        IF ( IIBC.EQ.6 .OR. IIBC.EQ.7 ) THEN
C          .
C          . First we make an early escape if r_{inner bnd} = 0
C          . this is equivlent to setting the function to zero
C          . which is what the current status of the matrix is
C          .
           IF ( ABS( X0 ).LT.TOL ) RETURN
C          .
C          . OK - so it's non-trivial!
C          . First, let's allocate the correct value of FAC
C          .
           IF ( IIBC.EQ.6 ) FAC = -1.0d0
           IF ( IIBC.EQ.7 ) FAC = -1.0d0*DBLE( L )
C          .
C          . Just quickly check that NALF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NALF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IIBC = ',IIBC,' and NALF = ',NALF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . Check for imminent division by zero
C          . We can safely divide by X0 now since the
C          . X0 = 0.0 case has been discounted above.
C          .
           QUOT = FAC/X0 + WMAT( 2, 1 )
           IF ( ABS( QUOT ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' QUOT = ', QUOT
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 2, NDS2
             DMAT( 1, I ) = (-1.0d0)*WMAT( 2, I )/QUOT
           ENDDO
           RETURN
C          .
        ENDIF
C       .
      ELSE
C       .
C       . We are considering the outer boundary
C       . First, return if IOBC = 1 or 2
C       .
        IF ( IOBC.EQ.1 .OR. IOBC.EQ.2 ) RETURN
C       .
C       . So we need to calculate deriv.s at
C       . the outer boundary, XARR( NR ).
C       . NDS2 is the number of nodes we can use at
C       . the outer boundary.
C       .
        IF ( IOBC.EQ.5 ) THEN
          IHND = 2
        ELSE
          IHND = 1
        ENDIF
C       .
        NDS2 = 1 + NNDS - NARF
C       .
        IF ( (NDS2-1).LT.IHND ) THEN
           PRINT *,' Subroutine LDGNMF '
           PRINT *,' NDS2                = ', NDS2
           PRINT *,' Required derivative = ', IHND
           PRINT *,' Program aborted.'
           STOP
        ENDIF
C       .
        DO I = 1, NDS2
          WORK1( I ) = XARR( NR - NDS2 + I )
        ENDDO
        X0 = XARR( NR )
        CALL GFDCFD ( X0, WORK1, NDS2, WMAT, NCFM, IPCM, WORK2 )
C       .
C       . WMAT( m + 1, i ) now contains the coefficient
C       . with which you multiply f_{nr-nds2+i} to get the m^{th}
C       . derivative of f at x_0.
C       .
C       . Consider the case IOBC = 3
C       . We require the first derivative = 0
C       .
        IF ( IOBC.EQ.3 ) THEN
C          .
C          . Just quickly check that NARF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NARF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IOBC = ',IOBC,' and NARF = ',NARF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . Check for imminent division by zero
C          .
           IF ( ABS( WMAT( 2, NDS2 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT( 2,',NDS2,') = ',WMAT( 2, NDS2 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          .
           DO I = 1, NDS2-1
             DMAT( NNDS, NNDS - NDS2 + I ) =
     1             (-1.0d0)*WMAT( 2, I )/WMAT( 2, NDS2 )
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . Consider the case IOBC = 4(5)
C       . We need both the function and the first
C       . (second) derivative to vanish. We do not need to
C       . to change anything to make the function
C       . vanish so just need do the first (second) deriv. cond.
C       .
        IF ( IOBC.EQ.4 .OR. IOBC.EQ.5 ) THEN
C          .
C          . Check for imminent division by zero
C          .
           IF ( IOBC.EQ.4 ) ID1 = 2
           IF ( IOBC.EQ.5 ) ID1 = 3
C          .
           IF ( ABS( WMAT( 2, NDS2-1 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT(',ID1,',',NDS2-1,') = ',WMAT( ID1,NDS2-1)
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 1, NDS2-2
             DMAT( NNDS + 1 - NARF, NNDS - NDS2 - NARF + 2 + I ) =
     1                    (-1.0d0)*WMAT( ID1, I )/WMAT( ID1, NDS2-1)
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . We may treat the cases IOBC.EQ.6 and IOBC.EQ.7
C       . together as they are both conditions of the form
C       . rdf/dr + FAC f(r) = 0
C       .
        IF ( IOBC.EQ.6 .OR. IOBC.EQ.7 ) THEN
C          .
C          . Just quickly check that NARF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NARF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IOBC = ',IOBC,' and NARF = ',NARF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . First, let's allocate the correct value of FAC
C          .
           IF ( IOBC.EQ.6 ) FAC = -1.0d0
           IF ( IOBC.EQ.7 ) FAC = DBLE( L + 1 )
C          .
C          . Check for imminent division by zero
C          . (Once again, we safely divide by X0)
C          .
           QUOT = FAC/X0 + WMAT( 2, NDS2 )
           IF ( ABS( QUOT ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' QUOT = ', QUOT
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          .
           DO I = 1, NDS2-1
             DMAT( NNDS, NNDS - NDS2 + I ) =
     1             (-1.0d0)*WMAT( 2, I )/QUOT
           ENDDO
           RETURN
C          .
        ENDIF
C       .
      ENDIF
C     .
      PRINT *,' Subroutine LDGNMF '
      PRINT *,' Problem with inputs.'
      PRINT *,' IIBC = ', IIBC
      PRINT *,' IOBC = ', IOBC
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
