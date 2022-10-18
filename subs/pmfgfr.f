C*********************************************************************
C subroutine Potential Magnetic Field Grid Formation Routine *********
C            -         -        -     -    -         -       *********
C Steve Gibbons Tue Mar 27 10:48:45 MET DST 2001                     C
C____________________________________________________________________C
C                                                                    C
C Takes an array of spherical harmonic coefficients, SHC, in order   C
C                                                                    C
C   SHC(   1   )   = g_1^0                                           C
C   SHC(   2   )   = g_1^1                                           C
C   SHC(   3   )   = h_1^1                                           C
C   SHC(   4   )   = g_2^0                                           C
C   SHC(   5   )   = g_2^1                                           C
C   SHC(   6   )   = h_2^1                                           C
C   SHC(   7   )   = g_2^2                                           C
C   SHC(   8   )   = h_2^2                                           C
C   SHC(   9   )   = g_3^0                                           C
C                                                                    C
C  ... and in general ...                                            C
C                                                                    C
C   SHC( INDSHC( L, M, 1) ) = g_L^M                                  C
C   SHC( INDSHC( L, M, 2) ) = h_L^M                                  C
C                                                                    C
C and calculates a 3D array DPOTMF( NRAD, NPHI, NTHE )               C
C                                                                    C
C where DPOTMF( irad, iphi, ithe ) = B_r( r, theta, phi )            C
C                                                                    C
C given that                                                         C
C                                                                    C
C      r = (RADEAR - RADCMB)(IRAD - 1)/(NRAD-1) + RADCMB             C
C                                                                    C
C  theta = pi*( ITHE - 1 )/(NTHE - 1)                                C
C                                                                    C
C    phi = 2*pi*( IPHI - 1 )/(NPHI - 1)                              C
C                                                                    C
C and B_r =  sum{l,m} [-g_l^m*(l+1)*(R/r)^{l+2} P_l^m*(cos m*phi)    C
C          + sum{l,m} [-h_l^m*(l+1)*(R/r)^{l+2} P_l^m*(sin m*phi)    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NRAD      : Number of radial grid nodes.                       C
C     NPHI      : Number of points in longitude.                     C
C     NTHE      : Number of points in latitude.                      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SHC       : Spherical harmonic coefficients (see above).       C
C     RADEAR    : Radius of the Earth.                               C
C     RADCMB    : Radius of the core-mantle boundary.                C
C                                                                    C
C     DPOTMF    : Array Dim ( NRAD, NPHI, NTHE ).                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PMFGFR( LH, NRAD, NPHI, NTHE, SHC, RADEAR,
     1                   RADCMB, DPOTMF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          LH, NRAD, NPHI, NTHE
      DOUBLE PRECISION SHC( * ), RADEAR, RADCMB,
     1                 DPOTMF( NRAD, NPHI, NTHE )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IRAD, ITHE, IPHI, IH, NH, L, M, ICS
      DOUBLE PRECISION DRAD, DPHI, DTHE, RAD1, RAD2, PHI1, PHI2,
     1                 THE1, THE2, ZERO, DELTA, PI, SHMPLG, RADFAC,
     2                 FAC, DLP2, DLP1, RAD, THE, PHI, COEF, COSTH,
     3                 DPLM, PHITRM
      EXTERNAL         SHMPLG
      PARAMETER ( ZERO = 0.0d0, DELTA = 1.0d-5,
     1            PI = 3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO ITHE = 1, NTHE
        DO IPHI = 1, NPHI
          DO IRAD = 1, NRAD
            DPOTMF( IRAD, IPHI, ITHE ) = ZERO
          ENDDO
        ENDDO
      ENDDO
C     .
      RAD1   = RADCMB
      RAD2   = RADEAR
      DRAD   = (RAD2 - RAD1)/DBLE( NRAD - 1 )
C     .
      THE1   = DELTA
      THE2   = PI - DELTA
      DTHE   = (THE2 - THE1)/DBLE( NTHE - 1 )
C     .
      PHI1   = ZERO
      PHI2   = 2.0d0*PI
      DPHI   = (PHI2 - PHI1)/DBLE( NPHI - 1 )
C     .
      NH     = LH*(LH + 2)
C     .
      DO ITHE = 1, NTHE
        THE   = THE1 + DBLE( ITHE - 1 )*DTHE
        COSTH = DCOS( THE )
        DO IH = 1, NH
          CALL LMFIND( IH, L, M, ICS )
          COEF  = SHC( IH )
          DPLM  = SHMPLG( L, M, COSTH )
          DLP1  = DBLE( L+1 )
          DLP2  = DBLE( L+2 )
C         .
          DO IPHI = 1, NPHI
            PHI    = PHI1 + DBLE( IPHI - 1 )*DPHI
            IF ( ICS.EQ.1 ) PHITRM = DCOS( PHI*M )
            IF ( ICS.EQ.2 ) PHITRM = DSIN( PHI*M )
C           .
            DO IRAD = 1, NRAD
              RAD    = RAD1 + DBLE( IRAD - 1 )*DRAD
              RADFAC = (RADEAR/RAD)**DLP2
C             .
              FAC = COEF*RADFAC*DLP1*DPLM*PHITRM*(-1.0d0)
              DPOTMF( IRAD, IPHI, ITHE ) = 
     1                        DPOTMF( IRAD, IPHI, ITHE ) + FAC
            ENDDO
C           . end loop irad = 1, nrad
          ENDDO
C         . end loop iphi = 1, nphi
        ENDDO
C       . end loop ih = 1, nh
      ENDDO
C     . end loop ithe = 1, nthe
      RETURN
      END
C*********************************************************************
