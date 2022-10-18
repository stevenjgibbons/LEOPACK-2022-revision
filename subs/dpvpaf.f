C*********************************************************************
C subroutine Double Precision Value Position in Array Find ***********
C            -      -         -     -           -     -    ***********
C Steve Gibbons Wed Jun 14 13:46:56 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C If an array A currently has N elements and may contain a maximum   C
C of NMAX elements. A double precision scalar, VAL, is compared to   C
C every element of A and if DABS( VAL - A[i] ).GT.TOL for all i=1,N  C
C then VAL is made the N+1^{th} element of the array, A is then      C
C reordered using SHELL1 and IPOS is returned -1.                    C
C N is returned as the new value (i.e. N+1).                         C
C                                                                    C
C If DABS( VAL - A[ipos] ).LE.TOL for some i then no change is made  C
C and IPOS is returned.                                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N         : Elements currently in A.                           C
C     NMAX      : Limit for N. Dimension of A.                       C
C     IPOS      : See above.                                         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VAL       : New value to be investigated.                      C
C     A         : Array with maximum dimension (NMAX).               C
C     TOL       : Distance criterion. (See above).                   C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPVPAF( N, NMAX, VAL, A, TOL, IPOS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N, NMAX, IPOS
      DOUBLE PRECISION VAL, A( * ), TOL
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IFLAG
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0, IFLAG = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( TOL.LT.ZERO ) THEN
        PRINT *,' Subroutine DPVPAF.'
        PRINT *,' TOL = ', TOL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C First treat first element case ...
C
      IF ( N.LE.0 ) THEN
        IF ( NMAX.GT.0 ) THEN
          A( 1 ) = VAL
          N      = 1
          IPOS   = -1
          RETURN
        ELSE
          PRINT *,' Subroutine DPVPAF.'
          PRINT *,' NMAX = ', NMAX
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDIF
C
C OK so N is above zero
C
      DO I = 1, N
        IF ( DABS( VAL - A( I ) ).LT.TOL ) THEN
          IPOS = I
          RETURN
        ENDIF
      ENDDO
C
C If we reach this point then VAL is distinct from
C all elements of A sp we need to add VAL to A.
C
      IF ( NMAX.LE.N ) THEN
        PRINT *,' NMAX = ', NMAX,' N = ', N
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      N      = N+1
      A( N ) = VAL
      IPOS   = -1
      CALL SHELL1( N, A, IFLAG )
C
      RETURN
      END
C*****************************************************************
