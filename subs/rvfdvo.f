C*********************************************************************
C Subroutine Radial Vector Function Dot Vector function: Optimised ***
C            -      -      -        -   -                -         ***
C Steve Gibbons Fri Jul 14 10:25:41 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C If RVFO is a radial vector function, and VFO a vector function     C
C then RVFDVO calculates the scalar product between the two          C
C functions at a given radial grid nodes, IR.                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IR        : Number of current radial grid node.                C
C     NR        : Total number of radial grid nodes.                 C
C     NPHPTS    : Number of phi points.                              C
C     NTHPTS    : Number of theta points.                            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RVFO      : Radial vector func. Dim ( 3, NPHPTS , NTHPTS, NR ) C
C     VFO       : Vector function. Dim ( 3, NPHPTS, NTHPTS )         C
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SF	: Scalar function. Dim ( NPHPTS, NTHPTS )            C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RVFDVO( RVFO, VFO, SF, IR, NR, NPHPTS, NTHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHPTS, NTHPTS, IR, NR
      DOUBLE PRECISION RVFO( 3, NPHPTS, NTHPTS, NR ),
     1                 VFO( 3, NPHPTS, NTHPTS ),
     2                 SF( NPHPTS, NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHETA, IPHI
      DOUBLE PRECISION DTEMP1, DTEMP2, DTEMP3
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C .... let's loop around icomp, then nthpts, then nphpts
      DO ITHETA = 1, NTHPTS
        DO IPHI = 1, NPHPTS
          DTEMP1 = RVFO( 1, IPHI, ITHETA, IR )*VFO( 1, IPHI, ITHETA )
          DTEMP2 = RVFO( 2, IPHI, ITHETA, IR )*VFO( 2, IPHI, ITHETA )
          DTEMP3 = RVFO( 3, IPHI, ITHETA, IR )*VFO( 3, IPHI, ITHETA )
          SF( IPHI, ITHETA ) = DTEMP1 + DTEMP2 + DTEMP3
        ENDDO
      ENDDO
C____________________________________________________________________C
      RETURN
      END
C*********************************************************************
