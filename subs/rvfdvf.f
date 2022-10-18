C*********************************************************************
C Subroutine Radial Vector Function Dot Vector Function **************
C            -      -      -        -   -      -        **************
C Steve Gibbons Tue May  9 17:48:22 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C If RVF is a radial vector function, and VF a vector function       C
C then RVFDVF calculates the scalar product between the two          C
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
C     NPHPTS	: Number of phi points.                              C
C     NTHPTS	: Number of theta points.                            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RVF	: Radial vector func. Dim ( NR, NPHPTS , NTHPTS, 3 ) C
C     VF 	: Vector function. Dim ( NPHPTS , NTHPTS, 3 )        C
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
      SUBROUTINE RVFDVF ( RVF, VF, SF, IR, NR, NPHPTS, NTHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHPTS, NTHPTS, IR, NR
      DOUBLE PRECISION RVF( NR, NPHPTS, NTHPTS, 3 ),
     1                 VF( NPHPTS , NTHPTS, 3 ),
     2                 SF( NPHPTS, NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHETA, IPHI, ICOMP, IOP
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C .... let's set SF to zero.
c
      IOP = 0
      CALL MATOP ( SF, ZERO, NPHPTS, NTHPTS, IOP )
c
C .... let's loop around icomp, then nthpts, then nphpts
      DO ICOMP = 1, 3
         DO ITHETA = 1, NTHPTS
            DO IPHI = 1, NPHPTS
               SF( IPHI, ITHETA ) = SF( IPHI, ITHETA ) +
     1    RVF( IR, IPHI, ITHETA, ICOMP )*VF( IPHI, ITHETA, ICOMP )
            ENDDO
         ENDDO
      ENDDO
C____________________________________________________________________C
      RETURN
      END
C*********************************************************************
