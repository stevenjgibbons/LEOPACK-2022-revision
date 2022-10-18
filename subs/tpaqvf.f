C*********************************************************************
C subroutine Three Point Amplitude Quadratic Value Find **************
C            -     -     -         -         -     -    **************
C Steve Gibbons Tue Oct 16 09:53:15 WEST 2001                        C
C____________________________________________________________________C
C                                                                    C
C We have the equation  f(x) = | m(x-c) |                            C
C                                                                    C
C Is it possible to find m and c from the values x1, f(x1),          C
C x2, f(x2) and x3, f(x3) ??                                         C
C                                                                    C
C Well if we square both sides of the equation, we get               C
C                                                                    C
C f^2 = m^2 ( x^2 - 2xc + c^2 )                                      C
C                                                                    C
C If we then solve the linear system                                 C
C                                                                    C
C  (f_1^2)  =  ( x_1^2   x_1    1 ) (   m^2   )                      C
C  (f_2^2)  =  ( x_2^2   x_2    1 ) ( -2m^2c  )                      C
C  (f_3^2)  =  ( x_3^2   x_3    1 ) ( m^2 c^2 )                      C
C                                                                    C
C We should be able to work out m and c?                             C
C                                                                    C
C  IERR = 1 means that m^2 is less than zero and we have a problem!! C
C                                                                    C
C  If IERR = 0, then RETC and RETM should contain the correct        C
C  values for M and C.                                               C
C                                                                    C
C*********************************************************************
      SUBROUTINE TPAQVF( X1, F1, X2, F2, X3, F3, RETC, RETM, IERR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          IERR
      DOUBLE PRECISION X1, F1, X2, F2, X3, F3, RETC, RETM
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER          NV, INFO, NRHS
      DOUBLE PRECISION DZERO, DLOW
      PARAMETER      ( NV = 3, DZERO = 0.0d0,
     1                 DLOW = 1.0d-7 )
      INTEGER          IWORK( NV )
      DOUBLE PRECISION COEFM( NV, NV ), WORK( NV )
      CHARACTER *(1)   TRANS
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      COEFM( 1, 1 )  = X1*X1
      COEFM( 1, 2 )  = -X1
      COEFM( 1, 3 )  = 1.0d0
C
      COEFM( 2, 1 )  = X2*X2
      COEFM( 2, 2 )  = -X2
      COEFM( 2, 3 )  = 1.0d0
C
      COEFM( 3, 1 )  = X3*X3
      COEFM( 3, 2 )  = -X3
      COEFM( 3, 3 )  = 1.0d0
C
      WORK( 1 )      = F1*F1
      WORK( 2 )      = F2*F2
      WORK( 3 )      = F3*F3
C
C Ok - this matrix is now ready for inversion -
C For this we use the LAPACK routines DGETRF and DGETRI
C First perform LU decomposition
C
      CALL DGETRF( NV, NV, COEFM, NV, IWORK, INFO )
C
C     . Check that LU decomposition has gone without
C     . problem.
C     .
C
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine TPAQVF.'
         PRINT *,' The LAPACK subroutine DGETRF has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in LU decomposition of COEFM matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Now perform solve on work
C     .
      TRANS = 'N'
      NRHS  = 1
C     .
      CALL DGETRS( TRANS, NV, NRHS, COEFM, NV, IWORK, WORK,
     1             NV, INFO )
C     .
C     . Check that solution has gone without problem.
C     .
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine TPAQVF.'
         PRINT *,' The LAPACK subroutine DGETRS has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in solving the COEFM matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
      IF ( WORK( 1 ).LT.DZERO )  THEN
        IERR = 1
      ELSE
        IERR = 0
        RETM = DSQRT( WORK( 1 ) )
        IF ( RETM.LT.DLOW ) THEN
          RETC = 0.0d0
        ELSE
          RETC = 0.5d0*WORK( 2 )/WORK( 1 )
        ENDIF
      ENDIF
C     .
      RETURN
      END
C*********************************************************************
