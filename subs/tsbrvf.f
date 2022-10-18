C*********************************************************************
C subroutine Toroidal Solid Body Rotation Vector Form ****************
C            -        -     -    -        -      -    ****************
C Steve Gibbons Wed Oct 31 13:40:16 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Using the double precision array XARR, TSBRVF forms a vector which C
C must be added to VT2 in order to time-step with an imposed         C
C drifting frame of reference.                                       C
C                                                                    C
C Use the routine TSBRVA to add the vector to VT2.                   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes in main soln.          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim (NR). Radial grid nodes of main solution.      C
C     TSBRV     : Dim (NR). Toroidal Solid Body Rotation Vector.     C
C     SBRSCL    : Solid body rotation scale.                         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE TSBRVF( NR, XARR, TSBRV, SBRSCL )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NR
      DOUBLE PRECISION XARR( NR ), TSBRV( NR ), SBRSCL
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IR
      DOUBLE PRECISION DFAC, RAD
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DFAC = (-1.0d0)*DSQRT( 2.0d0 )*SBRSCL
      DO IR = 1, NR
        RAD         = XARR( IR )
        TSBRV( IR ) = DFAC*RAD
      ENDDO
C
      RETURN
      END
C*********************************************************************
