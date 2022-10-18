C*********************************************************************
C subroutine Vector Function CORiolis force **************************
C            -      -        ---            **************************
C Steve Gibbons 5.5.97                                               C
C____________________________________________________________________C
C Calculates coriolis vector CORV = K x V of V.                      C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VF        : Vector Function. An array of dimensions            C
C                  ( NPHPTS, NTHPTS, 3) which contain the R, THETA   C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IPHI, ITHETA, 2 ) is the Theta      C
C                  compontent of the vector at (iphi, itheta).       C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VFCOR ( NTHPTS, NPHPTS, VF, GAUX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHPTS, NPHPTS
      DOUBLE PRECISION VF ( NPHPTS, NTHPTS, 3)
      DOUBLE PRECISION GAUX ( NTHPTS )
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      INTEGER ITHETA, IPHI
      DOUBLE PRECISION VR,VTHE,VPHI,COSTH,SINTH,THETA
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C ............. loop around theta and phi points .....................
      DO ITHETA = 1, NTHPTS
         COSTH = GAUX( ITHETA )
         THETA = ACOS( COSTH )
         SINTH = DSIN( THETA )
         DO IPHI = 1, NPHPTS
            VR = VF ( IPHI, ITHETA, 1 )
            VTHE = VF ( IPHI, ITHETA, 2 )
            VPHI = VF ( IPHI, ITHETA, 3 )
            VF ( IPHI, ITHETA, 1 ) = -SINTH*VPHI
            VF ( IPHI, ITHETA, 2 ) = -COSTH*VPHI
            VF ( IPHI, ITHETA, 3 ) = COSTH*VTHE + SINTH*VR
         ENDDO
      ENDDO
C ............. ended looping around theta, phi points ...............
      RETURN
      END
C*********************************************************************
