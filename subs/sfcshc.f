C*********************************************************************
C double precision function Scalar Function Constant Sph. Harm. Find *
C                           -      -        -        -    -     -    *
C Steve Gibbons Tue Nov  7 15:19:30 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C If the array SF contains a scalar function f( theta, phi )         C
C then the integral                                                  C
C                                                                    C
C \int_0^{2 pi} \int_0^{pi} f( theta, phi ) sin \theta d\theta d\phi C
C  = 4.\pi SFCSHC( GAUX, SF, NTHPTS, NPHPTS )                        C
C                                                                    C
C See argument list for format of SF and GAUX.                       C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SF	: Scalar Function. This is an array containing a     C
C                  function over a set of theta and phi points.      C
C                    Dimensions are                                  C
C                      ( NPHPTS , NTHPTS )                           C
C                                                                    C
C                 SF( i, j ) contains the value a scalar function    C
C                 at theta = acos( x_j ) where x_j is the j^{th}     C
C                 element of GAUX returned by GAUWTS and             C
C                 phi = 2 pi ( i - 1 ) / nphpts                      C
C                                                                    C
C     GAUW      : Gauss weights as calculated by GAUWTS. ( NTHPTS )  C
C     SFCSHC    : Coefficient of the l=0 spherical harmonic.         C
C                                                                    C
C  Integer                                                           C
C  -------							     C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION SFCSHC( GAUW, SF, NTHPTS, NPHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NTHPTS, NPHPTS
      DOUBLE PRECISION SFCSHC, SF( NPHPTS, NTHPTS ), GAUW( NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER ITHETA,I
      DOUBLE PRECISION ZERO, WEIGHT, FAC
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      SFCSHC = ZERO
      FAC    = 0.5d0/DBLE( NPHPTS )
C ............... begin looping around theta points ..................
      DO ITHETA = 1, NTHPTS
        WEIGHT = GAUW( ITHETA )*FAC
        DO I = 1, NPHPTS
           SFCSHC = SFCSHC + WEIGHT*SF( I, ITHETA )
        ENDDO
      ENDDO
C ................. end looping around theta points ..................
      RETURN
      END
C*********************************************************************
