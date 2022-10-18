C*********************************************************************
C subroutine SPHERical SATellite 2 WORLD coordinates transform *******
C            -----     ---       - -----                       *******
C Steve Gibbons Mon Mar 12 10:47:53 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Has input values R, THETA, PHI and EARCM.                          C
C R, THETA, PHI are the spherical coordinates (REAL) and EARCM       C
C is a 3 by 3 double precision matrix which has been prepared by     C
C EARCMC and, probably, inverted by EARMIR.                          C
C                                                                    C
C Spherical coordinates are in radians.                              C
C                                                                    C
C Returns coordinates XWORLD and YWORLD                              C
C DEPTH is the X2 co-ordinate.                                       C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SPHER_SAT_2_WORLD( R, THETA, PHI, EARCM, XWORLD,
     1                              YWORLD, DEPTH )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      REAL             R, THETA, PHI, XWORLD, YWORLD, DEPTH
      DOUBLE PRECISION EARCM( 3, 3 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      REAL             X1, X2, X3
      DOUBLE PRECISION ALPHA, BETA, XVEC( 3 ), YVEC( 3 )
      CHARACTER *(1)   TRANS
      INTEGER M, N, LDA, INCX, INCY
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      TRANS = 'N'
      M     = 3
      N     = 3
      LDA   = 3
      INCX  = 1
      INCY  = 1
      ALPHA = 1.0d0
      BETA  = 0.0d0
C
      X1 = R*SIN( THETA )*COS( PHI )
      X2 = R*SIN( THETA )*SIN( PHI )
      X3 = R*COS( THETA )
C
      XVEC( 1 ) = DBLE( X1 )
      XVEC( 2 ) = DBLE( X2 )
      XVEC( 3 ) = DBLE( X3 )
C
      CALL DGEMV ( TRANS, M, N, ALPHA, EARCM, LDA, XVEC, INCX,
     $                   BETA, YVEC, INCY )
C
      XWORLD = REAL( YVEC( 1 ) )*(-1.0)
      DEPTH  = REAL( YVEC( 2 ) )
      YWORLD = REAL( YVEC( 3 ) )
C
      RETURN
      END
C*********************************************************************
