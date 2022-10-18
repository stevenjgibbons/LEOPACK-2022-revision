C*********************************************************************
C subroutine Double precision VECtor Zero ****************************
C Steve Gibbons Wed Feb 14 09:01:11 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Sets to zero a double precision vector, VEC, of length N.          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N		: Length of the vector.                              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VEC	: Vector - dimension ( N )                           C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DVECZ( VEC, N )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          N
      DOUBLE PRECISION VEC( N )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          I
      DOUBLE PRECISION DZERO
      PARAMETER      ( DZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO I = 1, N
        VEC( I ) = DZERO
      ENDDO
      RETURN
      END
C*********************************************************************


