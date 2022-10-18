C*********************************************************************
C subroutine Equally Spaced Node Abscissa Allocation Subroutine. *****
C            -       -      -    -        -          -           *****
C Steve Gibbons Thu Oct 21 11:35:00 BST 1999                         C
C                                                                    C
C This short routine simply fills the array XARR with xvalues        C
C which are evenly spaced with                                       C
C                                                                    C
C  x_j = r_i + ( j - 1 )*h with h = ( r_o - r_i )/( nr - 1 )         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     RI        : Radius of inner boundary                           C
C     RO        : Radius of outer boundary                           C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ESNAAS( NR, XARR, RI, RO)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR
      DOUBLE PRECISION RI, RO, XARR( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR
      DOUBLE PRECISION H, TOL
      PARAMETER ( TOL = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NR.LT.4 ) THEN
         PRINT *,' Subroutine ESNAAS.'
         PRINT *,' NR = ', NR
         STOP
      ENDIF
C
      IF ( (RO-RI).LT.TOL ) THEN
         PRINT *,' Subroutine ESNAAS.'
         PRINT *,' RI = ', RI
         PRINT *,' RO = ', RO
         STOP
      ENDIF
C
      H = ( RO - RI )/DBLE( NR - 1 )
      DO IR = 1, NR
        XARR( IR ) = RI + DBLE( IR - 1 )*H
      ENDDO
C
      RETURN
      END
C*********************************************************************
