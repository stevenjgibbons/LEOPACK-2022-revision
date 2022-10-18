C*********************************************************************
C subroutine Solid Body Rotation Radial Function Constraint **********
C            -     -    -        -      -        -          **********
C Steve Gibbons Mon Nov 15 12:17:50 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C If we are solving the vorticity equation with a T_1^0 velocity     C
C harmonic with stress free boundaries, then we need to impose       C
C a condition of no solid body rotation.                             C
C                                                                    C
C SBRRFC finds out whether we have to impose this additional         C
C condition. If not, then OSBR is returned .FALSE. and no further    C
C action is taken. Otherwise, OSBR is set to .TRUE. and SBRRFC forms C
C a vector (lenth NR*NH) with which the vorticity equation may       C
C carefully be solved (subroutine VESR).                             C
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
C     MHM       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     SBRVEC    : Solid body rotation vector. Dim ( NR*NH )          C
C                                                                    C
C     WORK1     : Work array dim ( NNDS )                            C
C     WORK2     : Work array dim ( NNDS )                            C
C     COEFM     : Work array dim ( NNDS, NNDS )                      C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     OSBR      : Returned .FALSE. if there is no problem            C
C                 Returned .TRUE. if VESR must be used to solve      C
C                 the problem.                                       C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SBRRFC( NR, INARR, MHT, MHL, MHM, MHP, MHIBC, MHOBC,
     1                   NNDS, XARR, SBRVEC, OSBR, WORK1, WORK2,
     2                   IWORK, COEFM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MHT( * ), MHL( * ), MHM( * ), MHP( * ),
     1        MHIBC( * ), MHOBC( * ), NNDS, IWORK( NNDS )
      DOUBLE PRECISION XARR( NR ), SBRVEC( * ), WORK1( NNDS ),
     1                 WORK2( NNDS ), COEFM( NNDS, NNDS )
      LOGICAL OSBR
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR, NH, IH, IS, ILEN, IOP, IND, INDFUN, INODE
      DOUBLE PRECISION ZERO, LOW, X0, QUOT, HI
      PARAMETER ( ZERO = 0.0d0, LOW = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( INARR( 2 ).NE.NR ) THEN
        PRINT *,' Subroutine SBRRFC.'
        PRINT *,' NR = ', NR,' INARR( 2 ) = ', INARR( 2 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      NH = INARR( 3 )
C
      IF ( NNDS.LT.2 .OR. NNDS.GT.NR ) THEN
        PRINT *,' Subroutine SBRRFC.'
        PRINT *,' NNDS = ', NNDS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO IH = 1, NH
        IF ( MHT( IH ).EQ.2 .AND. MHL( IH ).EQ.1 .AND.
     1       MHM( IH ).EQ.0 ) THEN
C          .
C          . This is the T_1^0 velocity harmonic
C          .
           IS = MHP( IH )
           IF ( MHIBC( IS ).EQ.6 .AND. MHOBC( IS ).EQ.6 ) THEN
             OSBR = .TRUE.
             GOTO 50
           ENDIF
C          .
        ENDIF
      ENDDO
      OSBR = .FALSE.
 50   CONTINUE
      IF ( .NOT. OSBR ) RETURN
C
C Zero the array SBRVEC
C
      ILEN = NR*NH
      IOP  = 0
      CALL VECOP( SBRVEC, ZERO, ILEN, IOP )
C
C Loop IR around 2, NR - 1 and fill in the
C "obvious values" ...
C
      DO IR = 2, NR - 1
        IND = INDFUN( IR, IH, INARR )
        HI = XARR( IR + 1 ) - XARR( IR - 1 )
        SBRVEC( IND ) = 0.5d0*HI/XARR( IR )
      ENDDO
C
C Now add components from inner boundary
C
      IF ( XARR( 1 ).GT.LOW ) THEN
        X0 = XARR( 1 )
        DO INODE = 1, NNDS
          WORK1( INODE ) = XARR( INODE )
        ENDDO
        CALL GFDCFD( X0, WORK1, NNDS, COEFM, NNDS, IWORK, WORK2 )
C
C The coeff of f( x_i ) in the derivative of f at xarr( 1 )
C is now stored in COEFM( 2, i )
C
        QUOT = 1.0d0 - X0*COEFM( 2, 1 )
        IF ( ABS( QUOT ).LT.LOW ) THEN
          PRINT *,' Subroutine SBRRFC.'
          PRINT *,' QUOT = ', QUOT
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
        HI = XARR( 2 ) - XARR( 1 )
        DO IR = 2, NNDS
          IND = INDFUN( IR, IH, INARR )
          SBRVEC( IND ) = SBRVEC( IND ) +
     1            0.5d0*HI*COEFM( 2, IR )/(QUOT*X0)
        ENDDO
C
      ENDIF
C
C Now add components from outer boundary
C
      X0 = XARR( NR )
      HI = X0 - XARR( NR - 1 )
      DO INODE = 1, NNDS
        WORK1( INODE ) = XARR( NR + 1 - INODE )
      ENDDO
      CALL GFDCFD( X0, WORK1, NNDS, COEFM, NNDS, IWORK, WORK2 )
C
C The coeff of f( x_i ) in the derivative of f at xarr( nr )
C is now stored in COEFM( 2, nr + 1 - i )
C
      QUOT = 1.0d0 - X0*COEFM( 2, 1 )
      IF ( ABS( QUOT ).LT.LOW ) THEN
        PRINT *,' Subroutine SBRRFC.'
        PRINT *,' QUOT = ', QUOT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO INODE = 2, NNDS
        IR = NR + 1 - INODE
        IND = INDFUN( IR, IH, INARR )
        SBRVEC( IND ) = SBRVEC( IND ) +
     1                 0.5d0*HI*COEFM( 2, INODE )/(QUOT*X0)
      ENDDO
C
      RETURN
      END
C*********************************************************************
