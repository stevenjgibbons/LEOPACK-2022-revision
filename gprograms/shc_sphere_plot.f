C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C Wed Oct 10 13:26:19 WEST 2001                                      C
C                                                                    C
C Spherical Harmonic Coefficient Sphere Plot                         C
C                                                                    C
C*********************************************************************
      PROGRAM shc_sphere_plot
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      REAL    THETA1, THETA2, PHI1, PHI2
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMC / THETA1, THETA2, PHI1, PHI2
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRADMX, NTHMAX, NPHMAX, NLEVM, LHMAX, NHMAX
      PARAMETER ( NRADMX = 1, NTHMAX = 250, NPHMAX = 250,
     1            NLEVM = 20, LHMAX = 160, 
     2            NHMAX = LHMAX*(LHMAX+2) )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      REAL    FUNCRAD( NPHMAX, NTHMAX ),
     1        FUNCP1( NRADMX, NTHMAX ),
     2        FUNCP2( NRADMX, NTHMAX ),
     3        FUNCTHE( NRADMX, NPHMAX )
      REAL    CONTL( NLEVM )
C
      DOUBLE PRECISION SHC( NHMAX ),
     1                 DFUNC( NPHMAX, NTHMAX )
C
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C     .
C     . General
C     .
      INTEGER          NLEV, LUIN, LU, IOP, IDEV,
     1                 LH, I, ILEN, IFORM, IAS,
     2                 NTHP, NTHEP, NPHI, NPHIP
      CHARACTER *(1)   ESCPCH
      CHARACTER *(180) STEM, LINE, FNAME, STEM2
      REAL             PI, ZERO, PAPWIDTH, DTHETA,
     1                 PHIV1, PHIV2, RADV1, RADV2
      DOUBLE PRECISION ALPHAD, BETAD, GAMMAD, EARCM( 3, 3 )
C     .
C Variables for setting the colour.
C
      REAL    HUEPOS, HUENEG, CS, SCAL
C     .
      PARAMETER ( PI=3.1415926, ZERO = 0.0 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      ESCPCH = '*'
      LUIN   = 5
      PRINT *,' Enter ALPHAD, BETAD, GAMMAD.'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) ALPHAD, BETAD, GAMMAD
      PRINT *,' ALPHAD = ', ALPHAD
      PRINT *,' BETAD  = ', BETAD
      PRINT *,' GAMMAD = ', GAMMAD
C
      CALL EARCMC( ALPHAD, BETAD, GAMMAD, EARCM )
      CALL EARMIR( EARCM )
C
      DO I = 1, 180
        STEM(I:I) = ' '
        STEM2(I:I) = ' '
      ENDDO
C
C First read output file-name stem
C
      LU     = 15
      PRINT *,' Enter filename stem.'
      CALL LINERD( LUIN, STEM, ESCPCH )
C
C Now read in spherical harmonic coefficients
C
      PRINT *,' Enter SHC file name.'
      CALL LINERD( LUIN, LINE, ESCPCH )
C
      DO I = 1, 180
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 52
        ENDIF
      ENDDO
 52   CONTINUE
      FNAME = LINE(1:ILEN)
C      
      IFORM = 1
      PRINT *,' Data file to read = ', FNAME
      CALL SHSFRD( LH, LHMAX, LU, FNAME, SHC, IFORM )
C
C Enter huepos, hueneg, csat, scal, iw
C
      PRINT *,' Enter huepos, hueneg, csat, scal, iw'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) HUEPOS, HUENEG, CS, SCAL
C
C Enter number of contour levels
C
      PRINT *,' Enter NLEV, IDEV, PAPWIDTH, IAS '
      PRINT *,' IDEV: 1 --> gif file. (Landscape)'
      PRINT *,'       2 --> gif file. (Portrait)'
      PRINT *,'       3 --> ps file. (Landscape)'
      PRINT *,'       4 --> ps file. (Portrait)'
      PRINT *,'       5 --> colour ps file. (Landscape)'
      PRINT *,'       6 --> colour ps file. (Portrait)'
      PRINT *,' IAS:  1 --> ignore all harmonics with non-zero m'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NLEV, IDEV, PAPWIDTH, IAS
C
      IF ( NLEV.GT.NLEVM ) THEN
        PRINT *,' NLEV  = ', NLEV 
        PRINT *,' NLEVM = ', NLEVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Enter resolutions and positions
C
      PRINT *,' Enter NTHP, NPHI '
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NTHP, NPHI
C
C Fill data arrays
C
      IOP = 0
      CALL RMATOP( FUNCRAD, ZERO, NPHMAX, NTHMAX, IOP )
      CALL RMATOP( FUNCP1, ZERO, NRADMX, NTHMAX, IOP )
      CALL RMATOP( FUNCP2, ZERO, NRADMX, NTHMAX, IOP )
      CALL RMATOP( FUNCTHE, ZERO, NRADMX, NPHMAX, IOP )
C
      CALL SHC2PT( LH, NPHI, NTHP, SHC, DFUNC, IAS )
C
      CALL DBLE_2_REAL( NPHI, NTHP, FUNCRAD, DFUNC )
C
      PHI1   = 0.0
      PHI2   = 2.0*PI
      THETA1 = 0.0
      THETA2 = PI
C
      DO I = 1, 180
        IF ( STEM(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 53
        ENDIF
      ENDDO
 53   CONTINUE
C
      STEM2(1:ILEN)        = STEM(1:ILEN)
      STEM2(ILEN+1:ILEN+1) = ' '
      STEM2                = STEM2(1:ILEN+6)
C
      NRAD = 0
      DTHETA = 0.0
      RADV1  = 1.0
      RADV2  = 1.0
C      
      CALL CUT_OUT_SPHERE_PLOT( NLEV, NPHI, NTHP, FUNCRAD,
     1    RADV1, RADV2, NRAD, NTHEP, FUNCP1, PHIV1, NRAD, NTHEP,
     2    FUNCP2, PHIV2, NRAD, NPHIP, FUNCTHE, DTHETA, NLEVM, CONTL,
     3    PAPWIDTH, EARCM, HUEPOS, HUENEG, CS, SCAL, IDEV, STEM2 )
C
      STOP
      END
C*********************************************************************

      SUBROUTINE DBLE_2_REAL( N1, N2, REALARR, DBLARR )
      IMPLICIT NONE
      INTEGER N1, N2
      REAL             REALARR( N1, N2 )
      DOUBLE PRECISION DBLARR( N1, N2 )
      INTEGER I1, I2
C
      DO I2 = 1, N2
        DO I1 = 1, N1
          REALARR( I1, I2 ) = REAL( DBLARR( I1, I2 ) )
        ENDDO
      ENDDO
C
      RETURN
      END
