C*********************************************************************
C subroutine Radial Vector Function Coriolis Force 2 *****************
C            -      -      -        -        -     - *****************
C Steve Gibbons Tue May  9 17:09:32 BST 2000                         C
C____________________________________________________________________C
C Calculates coriolis vector CORV = K x V of V.                      C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     NR        : Number of radial grid nodes.                       C
C     ILNR      : Lowest radial grid node.                           C
C     IRNR      : Highest radial grid node.                          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RVF       : Vector Function. An array of dimensions            C
C                ( NR, NPHPTS, NTHPTS, 3) which contain the R, THETA C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IR, IPHI, ITHETA, 2 ) is the Theta  C
C                  compontent of the vector at (ir, iphi, itheta).   C
C     RVF2      : Output Vector Function. Contains (k x v).          C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RVFCF2( NTHPTS, NPHPTS, NR, ILNR, IRNR, RVF,
     1                    RVF2, GAUX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHPTS, NPHPTS, NR, ILNR, IRNR
      DOUBLE PRECISION RVF( NR, NPHPTS, NTHPTS, 3),
     1                 RVF2( NR, NPHPTS, NTHPTS, 3)
      DOUBLE PRECISION GAUX ( NTHPTS )
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      INTEGER ITHETA, IPHI, IR
      DOUBLE PRECISION VR,VTHE,VPHI,COSTH,SINTH,THETA
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Loop around radial grid nodes
C
      DO IR = ILNR, IRNR
C      ............. loop around theta and phi points ................
       DO ITHETA = 1, NTHPTS
         COSTH = GAUX( ITHETA )
         THETA = ACOS( COSTH )
         SINTH = DSIN( THETA )
         DO IPHI = 1, NPHPTS
            VR   = RVF( IR, IPHI, ITHETA, 1 )
            VTHE = RVF( IR, IPHI, ITHETA, 2 )
            VPHI = RVF( IR, IPHI, ITHETA, 3 )
            RVF2( IR, IPHI, ITHETA, 1 ) = -SINTH*VPHI
            RVF2( IR, IPHI, ITHETA, 2 ) = -COSTH*VPHI
            RVF2( IR, IPHI, ITHETA, 3 ) = COSTH*VTHE + SINTH*VR
         ENDDO
       ENDDO
C ............. ended looping around theta, phi points ...............
      ENDDO
C
C End loop around radial grid nodes
C
      RETURN
      END
C*********************************************************************
