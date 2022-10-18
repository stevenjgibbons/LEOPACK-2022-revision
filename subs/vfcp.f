C*********************************************************************
C Subroutine Vector Function Cross Product ***************************
C            -      -        -     -       ***************************
C Steve Gibbons 9.7.97                                               C
C____________________________________________________________________C
C To have two input vector function arrays of dimension              C
C ( NPHPTS , NTHPTS, 3 ) of which the vector product is taken to     C
C give a thrid vector VF3 with dimension ( NPHPTS , NTHPTS, 3 )      C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NPHPTS	: Number of phi points.                              C
C     NTHPTS	: Number of theta points.                            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VF1	: First vector function. Dim ( NPHPTS , NTHPTS, 3 )  C
C     VF2	: Second vector function. Dim ( NPHPTS , NTHPTS, 3 ) C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VF3	: Third vector function. Dim ( NPHPTS , NTHPTS, 3 )  C
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
      SUBROUTINE VFCP ( VF1, VF2, VF3, NPHPTS, NTHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHPTS, NTHPTS
      DOUBLE PRECISION VF1 ( NPHPTS, NTHPTS, 3 ),
     1                 VF2 ( NPHPTS, NTHPTS, 3 ),
     2                 VF3 ( NPHPTS, NTHPTS, 3 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHETA, IPHI, IOP, I3
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0, I3 = 3 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C .... let's set VF3 to zero.
c
      IOP = 0
      CALL CUBEOP ( VF3, ZERO, NPHPTS, NTHPTS, I3, IOP )
c
C .... let's loop around icomp, then nthpts, then nphpts
c
      DO ITHETA = 1, NTHPTS
         DO IPHI = 1, NPHPTS
c
            VF3 ( IPHI, ITHETA, 1 ) =
     1       VF1 ( IPHI, ITHETA, 2 )*VF2 ( IPHI, ITHETA, 3 ) -
     2       VF1 ( IPHI, ITHETA, 3 )*VF2 ( IPHI, ITHETA, 2 )
c
            VF3 ( IPHI, ITHETA, 2 ) =
     1       VF1 ( IPHI, ITHETA, 3 )*VF2 ( IPHI, ITHETA, 1 ) -
     2       VF1 ( IPHI, ITHETA, 1 )*VF2 ( IPHI, ITHETA, 3 )
c
            VF3 ( IPHI, ITHETA, 3 ) =
     1       VF1 ( IPHI, ITHETA, 1 )*VF2 ( IPHI, ITHETA, 2 ) -
     2       VF1 ( IPHI, ITHETA, 2 )*VF2 ( IPHI, ITHETA, 1 )
c
         ENDDO
      ENDDO
c 
C____________________________________________________________________C
      RETURN
      END
C*********************************************************************
