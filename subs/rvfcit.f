C*********************************************************************
C subroutine Radial Vector Function Coriolis and Inertial Terms ******
C            -      -      -        -            -        -     ******
C Steve Gibbons Fri Jul 14 12:48:26 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Let RVF1 contain a velocity function, V; let RVF2 contain the      C
C curl of V. GAUX( ith ) contains the cosine of the colatitude of    C
C the THETA point number ith.                                        C
C                                                                    C
C RVF3 is returned containing                                        C
C                                                                    C
C (-CG)*( k x V ) + CF*( V x curl V )                                C
C                                                                    C
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
C     RVF1      : Velocity Function. An array of dimensions          C
C               ( 3, NPHPTS, NTHPTS, NR ) which contain the R, THETA C
C                  and PHI components of a VECTOR at each point      C
C                ... i.e. RVF1 ( 2, IPHI, ITHETA, IR ) is the Theta  C
C                  compontent of the vector at (ir, iphi, itheta).   C
C     RVF2      : Curl of velocity.                                  C
C     RVF3      : Output Vector Function. Contains                   C
C                       (-CG)*( k x V ) + CF*( V x curl V )          C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RVFCIT( NTHPTS, NPHPTS, NR, ILNR, IRNR, RVF1,
     1                    RVF2, RVF3, GAUX, CG, CF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHPTS, NPHPTS, NR, ILNR, IRNR
      DOUBLE PRECISION RVF1( 3, NPHPTS, NTHPTS, NR ),
     1                 RVF2( 3, NPHPTS, NTHPTS, NR ),
     2                 RVF3( 3, NPHPTS, NTHPTS, NR )
      DOUBLE PRECISION GAUX ( NTHPTS ), CG, CF
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      INTEGER ITHETA, IPHI, IR
      DOUBLE PRECISION VR, VTHE, VPHI, COSTH, SINTH,
     1                 OUTRAD, OUTTHE, OUTPHI,
     2                 RV2RAD, RV2THE, Rv2PHI
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
         SINTH = DSQRT( 1.0d0 - COSTH*COSTH )
         DO IPHI = 1, NPHPTS
            VR   = RVF1( 1, IPHI, ITHETA, IR )
            VTHE = RVF1( 2, IPHI, ITHETA, IR )
            VPHI = RVF1( 3, IPHI, ITHETA, IR )
            OUTRAD = CG*SINTH*VPHI
            OUTTHE = CG*COSTH*VPHI
            OUTPHI = (-1.0d0)*CG*COSTH*VTHE - CG*SINTH*VR
c
            RV2RAD = RVF2( 1, IPHI, ITHETA, IR )*CF
            RV2THE = RVF2( 2, IPHI, ITHETA, IR )*CF
            RV2PHI = RVF2( 3, IPHI, ITHETA, IR )*CF
c
            RVF3 ( 1, IPHI, ITHETA, IR ) = OUTRAD +
     1           VTHE*RV2PHI - VPHI*RV2THE
c
            RVF3 ( 2, IPHI, ITHETA, IR ) = OUTTHE +
     1           VPHI*RV2RAD - VR*RV2PHI
c
            RVF3 ( 3, IPHI, ITHETA, IR ) = OUTPHI +
     1           VR*RV2THE   - VTHE*RV2RAD
c
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
