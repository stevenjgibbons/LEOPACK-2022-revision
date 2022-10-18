C*********************************************************************
C subroutine Adapted Matrix Landau Couette Flow Terms ****************
C            -       -      -      -       -    -     ****************
C Steve Gibbons Sat Apr 29 15:21:33 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Supplies coefficients to the routine AMLICA for terms with the     C
C linear stability matrix for the equation                           C
C                                                                    C
C  A_t = ( lambda + 2 i kappa x - x^2 ) A + d^2A/dx^2                C
C                                                                    C
C This routine must be declared EXTERNAL in the routine which calls  C
C AMLICA.                                                            C
C                                                                    C
C The solution vector has two radial functions ( inarr(3) = 2 )      C
C the first is a_r and the second is a_i.                            C
C                                                                    C
C At present there are four options for adding terms to this         C
C matrix.                                                            C
C                                                                    C
C OPTION    ROW   COLUMN       TERM                                  C
C ------    ---   ------       ----                                  C
C                                                                    C
C   1       a_r    a_r        lambda  - x^2  + d^/dx^2               C
C   2       a_i    a_i        lambda  - x^2  + d^/dx^2               C
C   3       a_r    a_i         - 2 kappa x                           C
C   4       a_i    a_r         + 2 kappa x                           C
C                                                                    C
C More options may be added later for the time-step matrix terms.    C
C                                                                    C
C  Set IPARS( 1 ) to the OPTION number above.                        C
C                                                                    C
C  IFLAG = 1 or 2 -->  CVEC( 1 ) = LAMBDA - RAD*RAD                  C
C                      CVEC( 2 ) = 0.0d0                             C
C                      CVEC( 3 ) = 1.0d0                             C
C                                                                    C
C  IFLAG = 3      -->  CVEC( 1 ) = (-2.0d0)*KAPPA*RAD                C
C                                                                    C
C  IFLAG = 4      -->  CVEC( 1 ) = 2.0d0*KAPPA*RAD                   C
C                                                                    C
C  DPARS( 1 ) = LAMBDA                                               C
C  DPARS( 2 ) = KAPPA                                                C
C                                                                    C
C  IHD is set to 2 for IPARS( 1 ) = 1 and IPARS( 1 ) = 2.            C
C  IHD is set to 0 for IPARS( 1 ) = 3 and IPARS( 1 ) = 4.            C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AMLCFT( CVEC, RAD, IPARS, DPARS, IHD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IPARS( * ), IHD
      DOUBLE PRECISION CVEC( * ), RAD, DPARS( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFLAG, ND
      DOUBLE PRECISION LAMBDA, KAPPA
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFLAG = IPARS( 1 )
      IF (  ( IFLAG.EQ.1 .AND. IHD.NE.2 ) .OR.
     1      ( IFLAG.EQ.2 .AND. IHD.NE.2 ) .OR.
     2      ( IFLAG.EQ.3 .AND. IHD.NE.0 ) .OR.
     3      ( IFLAG.EQ.4 .AND. IHD.NE.0 )      ) THEN
        PRINT *,' Subroutine AMLCFT.'
        PRINT *,' IFLAG = ', IFLAG,' IHD = ', IHD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      LAMBDA = DPARS( 1 )
      KAPPA  = DPARS( 2 )
C
C     . Check for valid value of IFLAG
C     .
      IF ( IFLAG.LT.1 .OR. IFLAG.GT.4 ) THEN
         PRINT *,' Subroutine AMLCFT.'
         PRINT *,' IFLAG = ', IFLAG
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Zero all coefficients up to IHD
C     .
      DO ND = 1, IHD + 1
        CVEC( ND ) = 0.0d0
      ENDDO
C     .
C     . Option 1 or 2:
C     .
      IF ( IFLAG.EQ.1 .OR. IFLAG.EQ.2 ) THEN
         CVEC( 1 ) = LAMBDA - RAD*RAD
         CVEC( 3 ) = 1.0d0
         RETURN
      ENDIF
C     .
C     . Option 3:
C     .
      IF ( IFLAG.EQ.3 ) THEN
         CVEC( 1 ) = (-2.0d0)*KAPPA*RAD
         RETURN
      ENDIF
C     .
C     . Option 4:
C     .
      IF ( IFLAG.EQ.4 ) THEN
         CVEC( 1 ) = 2.0d0*KAPPA*RAD
         RETURN
      ENDIF
C     .
      END
C*********************************************************************
