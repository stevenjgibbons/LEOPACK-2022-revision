C*********************************************************************
C Subroutine Radial Vector Function Cross Product 2 ******************
C            -      -      -        -     -       - ******************
C Steve Gibbons Fri May 12 09:19:56 BST 2000                         C
C____________________________________________________________________C
C To have two input vector function arrays of dimension              C
C ( NR, NPHPTS, NTHPTS, 3 ) of which B multiplied by the vector      C
C product is added to A multiplied by a thrid vector VF3 with        C
C dimension ( NR, NPHPTS, NTHPTS, 3 )                                C
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
C     RVF1	: First vector func. Dim ( NR, NPHPTS , NTHPTS, 3 )  C
C     RVF2	: Second vector func. Dim ( NR, NPHPTS , NTHPTS, 3 ) C
C     A         : Mulitplier of RVF3                                 C
C     B         : Mulitplier of RVF1 x RVF2                          C
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RVF3	: Third vector func. Dim ( NR, NPHPTS , NTHPTS, 3 )  C
C____________________________________________________________________C
C                                                                    C
C      Let u = ( u1, u2, u3 ) and v = ( v1, v2, v3 )                 C
C      If w = ( w1, w2, w3 ) = u x v    then                         C
C                                                                    C
C                  w1 = A*w1 + B*u2*v3 - B*u3*v2                     C
C                  w2 = A*w2 + B*u3*v1 - B*u1*v3                     C
C                  w3 = A*w3 + B*u1*v2 - B*u2*v1                     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RVFCP2 ( RVF1, RVF2, RVF3, NPHPTS, NTHPTS, NR,
     1                   ILNR, IRNR, A, B )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHPTS, NTHPTS, NR, ILNR, IRNR
      DOUBLE PRECISION RVF1 ( NR, NPHPTS, NTHPTS, 3 ),
     1                 RVF2 ( NR, NPHPTS, NTHPTS, 3 ),
     2                 RVF3 ( NR, NPHPTS, NTHPTS, 3 ), A, B
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHETA, IPHI, IOP, I3, IR
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0, I3 = 3 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C .... let's loop around ir, then nthpts, then nphpts
c
      DO IR = ILNR, IRNR
       DO ITHETA = 1, NTHPTS
        DO IPHI = 1, NPHPTS
c
       RVF3( IR, IPHI, ITHETA, 1 ) = A*RVF3( IR, IPHI, ITHETA, 1 ) +
     1  B*RVF1 ( IR, IPHI, ITHETA, 2 )*RVF2 ( IR, IPHI, ITHETA, 3 ) -
     2  B*RVF1 ( IR, IPHI, ITHETA, 3 )*RVF2 ( IR, IPHI, ITHETA, 2 )
c
       RVF3( IR, IPHI, ITHETA, 2 ) = A*RVF3( IR, IPHI, ITHETA, 2 ) +
     1  B*RVF1 ( IR, IPHI, ITHETA, 3 )*RVF2 ( IR, IPHI, ITHETA, 1 ) -
     2  B*RVF1 ( IR, IPHI, ITHETA, 1 )*RVF2 ( IR, IPHI, ITHETA, 3 )
c
       RVF3( IR, IPHI, ITHETA, 3 ) = A*RVF3( IR, IPHI, ITHETA, 3 ) +
     1  B*RVF1 ( IR, IPHI, ITHETA, 1 )*RVF2 ( IR, IPHI, ITHETA, 2 ) -
     2  B*RVF1 ( IR, IPHI, ITHETA, 2 )*RVF2 ( IR, IPHI, ITHETA, 1 )
c
        ENDDO
       ENDDO
      ENDDO
c 
C____________________________________________________________________C
      RETURN
      END
C*********************************************************************
