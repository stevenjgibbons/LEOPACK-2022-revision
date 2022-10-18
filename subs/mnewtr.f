C*********************************************************************
C subroutine Multi-dimensional NEWTon raphson Solve ******************
C            -                 ----           -     ******************
C Steve Gibbons Sat Oct  9 18:11:31 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Say we have NEQN equations and NEQN unknowns.                      C
C For IEQN = 1, NEQN; we seek a solution ( x_1 , x_2, ... , x_NEQN ) C
C to the equations                                                   C
C                                                                    C
C f_IEQN( x_1 , x_2, ... , x_NEQN ) = 0.                             C
C                                                                    C
C FUNC is a double precision function with the calling sequence      C
C                                                                    C
C FUNC( NEQN, IEQN, ID, XVEC, INTARR, DPRARR )                       C
C                                                                    C
C FUNC returns f_IEQN( x_1 , x_2, ... , x_NEQN ) when ID = 0.        C
C If ID = ICMP, FUNC returns the partial derivative of f_IEQN        C
C with respect to x_ICMP.                                            C
C                                                                    C
C XVEC is a vector of length NEQN with X(icmp) = the icmp^{th} var.  C
C                                                                    C
C INTPAR and DPRARR are unlimited arrays of integer and double       C
C precision variables respectively which are not referred to by      C
C MNEWTS other than to supply information to FUNC.                   C
C                                                                    C
C                                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FUNC      : Declared EXTERNALly. Has the form                  C
C                 FUNC( NEQN, IEQN, ID, XVEC, INTARR, DPRARR )       C
C                 (see above).                                       C
C                                                                    C
C     XVEC      : Vector containing initial estimates of the         C
C                 unknowns. On return, XVEC contains the solution    C
C                 provided that convergence is attained.             C
C                 Dim ( NEQN )                                       C
C                                                                    C
C     ERR       : User specified tolerance of error of the solution. C
C                 If the norm of the r.h.s. vector is less than      C
C                 ERR, convergence will be judged to have been       C
C                 achieved.                                          C
C                                                                    C
C     WORKV     : Work array dim ( NEQN )                            C
C     WORKA     : Work array dim ( NEQN, NEQN )                      C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NEQN      : Number of equations and unknowns.                  C
C     ITMX      : Maximum number of iterations permitted.            C
C                                                                    C
C     INFO      : Output variable. If solution converges             C
C                 succesfully then INFO returns the number of        C
C                 iterations required to converge - this is          C
C                 an integer between 1 and ITMX.                     C
C                                                                    C
C                 A failiure will result in a negative value         C
C                 of INFO.                                           C
C                                                                    C
C                 INFO = -1: Jacobian matrix was singular.           C
C                 The LAPACK routine DGETRF was unable to perform    C
C                 an LU decomposition.                               C
C                 The value of 'INFO' returned by DGETRF             C
C                 is returned from MNEWTR in IWORK( 1 )              C
C                                                                    C
C                 INFO = -2:                                         C
C                 The LAPACK routine DGETRS was unable to solve      C
C                 the matrix equation.                               C
C                 The value of 'INFO' returned by DGETRS             C
C                 is returned from MNEWTR in IWORK( 1 )              C
C                                                                    C
C                 INFO = -3: Maximum number of iterations exceeded.  C
C                                                                    C
C     IWORK     : Work array dim ( NEQN ).                           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MNEWTR( FUNC, XVEC, NEQN, ITMX, ERR, INFO, WORKV,
     1                   WORKA, IWORK, INTARR, DPRARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NEQN, ITMX, INFO, IWORK( NEQN ), INTARR( * )
      DOUBLE PRECISION FUNC, XVEC( NEQN ), WORKV( NEQN ), ERR,
     1                 WORKA( NEQN, NEQN ), DPRARR( * )
      EXTERNAL FUNC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION DPNORM
      INTEGER NOIT, ID, IEQN, IERR, NRHS
      CHARACTER *(1) TRANS
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NRHS = 1
      TRANS = 'N'
      NOIT = 0
C
 50   CONTINUE
      NOIT = NOIT + 1
C
C Check maximum iterations haven't been exceeded
C
      IF ( NOIT.GT.ITMX ) THEN
        INFO = -3
        RETURN
      ENDIF
C
C O.k. - we form the right hand side
C Store this in WORKV
C Also calculate the solution norm
C
      ID = 0
      DPNORM = 0.0d0
      DO IEQN = 1, NEQN
        WORKV( IEQN ) =
     1    (-1.0d0)*FUNC( NEQN, IEQN, ID, XVEC, INTARR, DPRARR )
        DPNORM = DPNORM + ABS( WORKV( IEQN ) )
      ENDDO
C
C Return if change in solution is sufficiently small
C
      IF ( DPNORM.LT.ERR ) THEN
        INFO = NOIT
        RETURN
      ENDIF
C
C Now form the Jacobian matrix
C WORKA( i, j ) must contain the derivative of F_i with
C respect to x_j
C
      DO ID = 1, NEQN
        DO IEQN = 1, NEQN
          WORKA( IEQN, ID ) = FUNC( NEQN, IEQN, ID, XVEC,
     1                              INTARR, DPRARR )
        ENDDO
      ENDDO
C
C Now perform an LU decomposition upon the matrix WORKA
C Use the LAPACK routine DGETRF
C
      CALL DGETRF( NEQN, NEQN, WORKA, NEQN, IWORK, IERR )
C
C Return if LU decomposition has been unsuccessful
C
      IF ( IERR.NE.0 ) THEN
        IWORK( 1 ) = IERR
        INFO = -1
        RETURN
      ENDIF
C
C Ok - the LU decomposition seems to have gone
C without a problem. So let's solve the equation
C WORKA . DX = WORKV
C Use LAPACK routine DGETRS
C
      CALL DGETRS( TRANS, NEQN, NRHS, WORKA, NEQN, IWORK,
     1             WORKV, NEQN, IERR )
C
C Return if solution has been unsuccessful
C
      IF ( IERR.NE.0 ) THEN
        IWORK( 1 ) = IERR
        INFO = -2
        RETURN
      ENDIF
C
C Update solution.
C
      DO ID = 1, NEQN
        XVEC( ID ) = XVEC( ID ) + WORKV( ID )
      ENDDO
C
      GOTO 50
C
      END
C*********************************************************************
