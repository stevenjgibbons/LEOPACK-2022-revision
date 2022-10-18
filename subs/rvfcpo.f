C*********************************************************************
C Subroutine Radial Vector Function Cross Product: Optimised *********
C            -      -      -        -     -                  *********
C Steve Gibbons Fri Jul 14 19:31:51 BST 2000                         C
C____________________________________________________________________C
C To have two input vector function arrays of dimension              C
C ( 3, NPHPTS , NTHPTS, NR ) of which the vector product is taken to C
C give a thrid vector VF3 with dimension ( 3, NPHPTS, NTHPTS, NR )   C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NPHPTS	: Number of phi points.                              C
C     NTHPTS	: Number of theta points.                            C
C     NR        : Number of radial grid points.                      C
C     ILNR      : Lowest radial grid point.                          C
C     IRNR      : Highest radial grid point.                         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RVF1	: First vector func. Dim ( 3, NPHPTS, NTHPTS, NR )   C
C     RVF2	: Second vector func. Dim ( 3, NPHPTS, NTHPTS, NR )  C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RVF3	: Third vector func. Dim ( 3, NPHPTS, NTHPTS, NR )   C
C____________________________________________________________________C
C                                                                    C
C      Let u = ( u1, u2, u3 ) and v = ( v1, v2, v3 )                 C
C      If w = ( w1, w2, w3 ) = u x v    then                         C
C                                                                    C
C                  w1 = u2 * v3 - u3 * v2                            C
C                  w2 = u3 * v1 - u1 * v3                            C
C                  w3 = u1 * v2 - u2 * v1                            C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RVFCPO( RVF1, RVF2, RVF3, NPHPTS, NTHPTS, NR,
     1                   ILNR, IRNR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHPTS, NTHPTS, NR, ILNR, IRNR
      DOUBLE PRECISION RVF1 ( 3, NPHPTS, NTHPTS, NR ),
     1                 RVF2 ( 3, NPHPTS, NTHPTS, NR ),
     2                 RVF3 ( 3, NPHPTS, NTHPTS, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHETA, IPHI, IOP, I3, IR
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0, I3 = 3 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C .... let's set RVF3 to zero.
c
      IOP = 0
      CALL QUADOP ( RVF3, ZERO, I3, NPHPTS, NTHPTS, NR, IOP )
c
C .... let's loop around ir, then nthpts, then nphpts
c
      DO IR = ILNR, IRNR
       DO ITHETA = 1, NTHPTS
        DO IPHI = 1, NPHPTS
c
         RVF3 ( 1, IPHI, ITHETA, IR ) =
     1    RVF1 ( 2, IPHI, ITHETA, IR )*RVF2 ( 3, IPHI, ITHETA, IR ) -
     2    RVF1 ( 3, IPHI, ITHETA, IR )*RVF2 ( 2, IPHI, ITHETA, IR )
c
         RVF3 ( 2, IPHI, ITHETA, IR ) =
     1    RVF1 ( 3, IPHI, ITHETA, IR )*RVF2 ( 1, IPHI, ITHETA, IR ) -
     2    RVF1 ( 1, IPHI, ITHETA, IR )*RVF2 ( 3, IPHI, ITHETA, IR )
c
         RVF3 ( 3, IPHI, ITHETA, IR ) =
     1    RVF1 ( 1, IPHI, ITHETA, IR )*RVF2 ( 2, IPHI, ITHETA, IR ) -
     2    RVF1 ( 2, IPHI, ITHETA, IR )*RVF2 ( 1, IPHI, ITHETA, IR )
c
        ENDDO
       ENDDO
      ENDDO
c 
C____________________________________________________________________C
      RETURN
      END
C*********************************************************************
