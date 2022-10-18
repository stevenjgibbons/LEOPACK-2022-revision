C*********************************************************************
C subroutine Spherical Harmonic Coefficients 2 Phi and Theta *********
C            -         -        -            - -       -     *********
C Steve Gibbons Wed Oct 10 12:42:10 WEST 2001                        C
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
C and calculates a 2D array DFUNC( NPHI, NTHE )                      C
C                                                                    C
C  theta = pi*( ITHE - 1 )/(NTHE - 1)                                C
C                                                                    C
C    phi = 2*pi*( IPHI - 1 )/(NPHI - 1)                              C
C                                                                    C
C       f =  sum{l,m} g_l^m P_l^m*(cos m*phi)                        C
C          + sum{l,m} h_l^m P_l^m*(sin m*phi)                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NPHI      : Number of points in longitude.                     C
C     NTHE      : Number of points in latitude.                      C
C     IAS       : Set to 1 to use only the m=0 harmonics.            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SHC       : Spherical harmonic coefficients (see above).       C
C     DFUNC     : Array Dim ( NPHI, NTHE ).                          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHC2PT( LH, NPHI, NTHE, SHC, DFUNC, IAS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          LH, NPHI, NTHE, IAS
      DOUBLE PRECISION SHC( * ), DFUNC( NPHI, NTHE )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ITHE, IPHI, IH, NH, L, M, ICS
      DOUBLE PRECISION DPHI, DTHE, PHI1, PHI2,
     1                 THE1, THE2, ZERO, DELTA, PI, SHMPLG,
     2                 FAC, THE, PHI, COEF, COSTH,
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
          DFUNC( IPHI, ITHE ) = ZERO
        ENDDO
      ENDDO
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
          IF ( IAS.EQ.1 .AND. M.NE.0 ) GOTO 50
          COEF  = SHC( IH )
          DPLM  = SHMPLG( L, M, COSTH )
C         .
          DO IPHI = 1, NPHI
C           .
            PHI    = PHI1 + DBLE( IPHI - 1 )*DPHI
            IF ( ICS.EQ.1 ) PHITRM = DCOS( PHI*M )
            IF ( ICS.EQ.2 ) PHITRM = DSIN( PHI*M )
            FAC = COEF*DPLM*PHITRM
            DFUNC( IPHI, ITHE ) = DFUNC( IPHI, ITHE ) + FAC
C           .
          ENDDO
C         . end loop iphi = 1, nphi
  50    CONTINUE
        ENDDO
C       . end loop ih = 1, nh
      ENDDO
C     . end loop ithe = 1, nthe
      RETURN
      END
C*********************************************************************
