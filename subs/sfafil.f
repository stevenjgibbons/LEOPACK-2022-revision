C*********************************************************************
C subroutine Scalar Function Array FILl ******************************
C            -      -        -     ---- ******************************
C Steve Gibbons Fri Sep 10 08:55:15 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Simply performs the task of looping around an array, SF, of        C
C dimensions ( NPHPTS, NTHPTS ) and filling it according to the      C
C EXTERNALly declared double precision function FUNC which has the   C
C calling sequence                                                   C
C                                                                    C
C  FUNC( THETA, PHI, INTPAR, DPPAR )                                 C
C                                                                    C
C  where THETA and PHI are the colatitude and longitude resp.        C
C  with angle measured in radians.                                   C
C  The arrays INTPAR and DPPAR are for the arbitrary use of          C
C  FUNC.                                                             C
C                                                                    C
C  SF( i, j ) contains the value of FUNC                             C
C  at theta = acos( x_j ) where x_j is the j^{th}                    C
C  element of GAUX (possibly returned by GAUWTS) (for NTHPTS)        C
C  and phi = 2 pi ( i - 1 ) / nphpts                                 C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NPHPTS    : Number of phi points.                              C
C     NTHPTS    : Number of theta points.                            C
C     INTPAR    : Array dim. ( * ) for use by FUNC                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX      : Array containing abscissae of the Gauss-Legendre   C
C                  NTHPTS-points quadrature formula.                 C
C                   (Dimension NTHPTS). Output from                  C
C                      subroutine GAUWTS (or maybe totally           C
C                        arbitrary if you do not wish to transform)  C
C                                                                    C
C     SF        : Scalar Function. This is an array containing a     C
C                  function over a set of theta and phi points.      C
C                    Dimensions are                                  C
C                      ( NPHPTS , NTHPTS ) See above.                C
C                                                                    C
C     FUNC      : Declared EXTERNAL in calling routine.              C
C                                                                    C
C     DPPAR     : Array dim. ( * ) for use by FUNC                   C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SFAFIL( NPHPTS, NTHPTS, INTPAR, GAUX,
     1                   SF, FUNC, DPPAR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHPTS, NTHPTS, INTPAR( * )
      DOUBLE PRECISION SF( NPHPTS, NTHPTS ), FUNC, DPPAR( * ),
     1                 GAUX( NTHPTS )
      EXTERNAL FUNC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHE, IPHI
      DOUBLE PRECISION PHI, PI, DPHI, COSTH, THETA
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DPHI = 2.0d0*PI/DBLE( NPHPTS )
      DO ITHE = 1, NTHPTS
       COSTH = GAUX( ITHE )
       THETA = ACOS( COSTH )
       DO IPHI = 1, NPHPTS
        PHI = DBLE( IPHI - 1 )*DPHI
        SF( IPHI, ITHE ) = FUNC( THETA, PHI, INTPAR, DPPAR )
       ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
