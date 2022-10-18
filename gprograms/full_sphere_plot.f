C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C Thu Nov  8 12:04:25 WET 2001                                       C
C                                                                    C
C Full Sphere Plot                                                   C
C                                                                    C
C Modification Steve Gibbons Thu Jan 15 12:11:26 MET 2004            C
C ICOMP may be set to the negative of its intended value in the      C
C input file in order to remove the mean value from the final        C
C array. The default is no change.                                   C
C                                                                    C
C*********************************************************************
      PROGRAM full_sphere_plot
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      REAL    THETA1, THETA2, PHI1, PHI2
      COMMON  / PARAMC / THETA1, THETA2, PHI1, PHI2
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRADMX, NTHMAX, NPHMAX, NLEVM, LHMAX, NHMAX,
     1        ISVMAX, NNDM
      PARAMETER ( NRADMX = 250, NTHMAX = 250, NPHMAX = 250,
     1            NLEVM = 20, LHMAX = 160, NHMAX = 3000,
     2            ISVMAX = NRADMX*NHMAX, NNDM = 6 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      REAL    FUNCRAD( NPHMAX, NTHMAX )
      REAL    CONTL( NLEVM )
C
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     2        IWORK( NNDM )
C
      DOUBLE PRECISION XARR( NRADMX ), VEC( ISVMAX ),
     1                 WORK1( NNDM ), WORK2( NNDM ),
     2                 COEFM( NNDM, NNDM )
C
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C     .
C     . General
C     .
      INTEGER          NLEV, LUIN, LU, IOP, IAS, NNDS,
     1                 IDEV, INARR( 3 ), NH, ICOMP, RMFLAG,
     2                 LH, NR, NR1, I, ILEN, NTHP, NPHI
      CHARACTER *(1)   ESCPCH
      CHARACTER *(180) STEM, LINE, FNAME, STEM2
      REAL             PI, ZERO, PAPWIDTH, RADV2
      DOUBLE PRECISION RAD
      DOUBLE PRECISION ALPHAD, BETAD, GAMMAD, EARCM( 3, 3 )
C     .
C Variables for setting the colour.
C
      REAL    HUEPOS, HUENEG, CS, SCAL
C     .
      PARAMETER ( PI=3.1415926, ZERO = 0.0, IAS = 0 )
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
C Read integers file:
C
      PRINT *,' Enter integers file name.'
      CALL LINERD( LUIN, LINE, ESCPCH )
C
      DO I = 1, 180
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 50
        ENDIF
      ENDDO
 50   CONTINUE
      FNAME = LINE(1:ILEN)
      CALL BIHFRD( NH, NHMAX, MHT, MHL, MHM, LU, FNAME )
      INARR( 3 ) = NH
      LH = 0
      DO I = 1, NH
        IF ( MHL( I ).GT.LH ) LH = MHL( I )
      ENDDO
      IF ( LH.GT.LHMAX ) THEN
        PRINT *,' LH = ', LH
        PRINT *,' LHMAX = ', LHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now read in vectors file
C
      PRINT *,' Enter vector file name.'
      CALL LINERD( LUIN, LINE, ESCPCH )
C
      DO I = 1, 180
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 51
        ENDIF
      ENDDO
 51   CONTINUE
      FNAME = LINE(1:ILEN)
C
      CALL SVFRD( INARR, LU, NRADMX, VEC, FNAME )
      NR1 = INARR( 2 )
C
C Now read in radial spacings filename
C
      PRINT *,' Enter radial spacings file name.'
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
      CALL XARRRD( NR, NRADMX, XARR, LU, FNAME )
C
      IF ( NR1.NE.NR ) THEN
        PRINT *,' Solution vector and radial node '
        PRINT *,' file claim differing numbers of '
        PRINT *,' grid nodes. Program aborted.'
        STOP
      ENDIF
C
C Enter huepos, hueneg, csat, scal
C
      PRINT *,' Enter huepos, hueneg, csat, scal'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) HUEPOS, HUENEG, CS, SCAL
C
C Enter number of contour levels
C
      PRINT *,' Enter NLEV, IDEV, NNDS, PAPWIDTH, ICOMP'
      PRINT *,' IDEV: 1 --> gif file. (Landscape)'
      PRINT *,'       2 --> gif file. (Portrait)'
      PRINT *,'       3 --> ps file. (Landscape)'
      PRINT *,'       4 --> ps file. (Portrait)'
      PRINT *,'       5 --> colour ps file. (Landscape)'
      PRINT *,'       6 --> colour ps file. (Portrait)'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NLEV, IDEV, NNDS, PAPWIDTH, ICOMP
C
C Modification SJG Thu Jan 15 12:11:26 MET 2004
      IF ( ICOMP.LT.0 ) THEN
        ICOMP  = -ICOMP
        RMFLAG = 1
      ELSE
        RMFLAG = 0
      ENDIF
C END Modification SJG Thu Jan 15 12:11:26 MET 2004
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
      PRINT *,' Enter NTHP, NPHI, RADV2'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NTHP, NPHI, RADV2
C
      IF ( NTHP.GT.NTHMAX .OR. NPHI.GT.NPHMAX ) THEN
        PRINT *,' NTHMAX = ', NTHMAX
        PRINT *,' NPHMAX = ', NPHMAX
        PRINT *,' NTHP   = ', NTHP
        PRINT *,' NPHI   = ', NPHI
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Fill data arrays
C
      IOP = 0
      CALL RMATOP( FUNCRAD, ZERO, NPHMAX, NTHMAX, IOP )
C
C Enter the data arrays
C First constant radius data
C
      PHI1   = 0.0
      PHI2   = 2.0*PI
      THETA1 = 0.0
      THETA2 = PI
C
      RAD   = DBLE( RADV2 )
      CALL CONSTANT_R_RECT_EVAL( ICOMP, INARR, NNDS, IWORK,
     1          MHT, MHL, MHM, LH, IAS, RAD, VEC, XARR, WORK1,
     2          WORK2, COEFM, FUNCRAD, NPHI, NTHP, THETA1, THETA2,
     3          PHI1, PHI2 )
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
      CALL ARR_RM_MEAN( NPHI, NTHP, FUNCRAD, RMFLAG )
C
      CALL SIMPLE_SPHERE_PLOT( NLEV, NPHI, NTHP, FUNCRAD,
     1    RADV2, NLEVM, CONTL, PAPWIDTH, EARCM, HUEPOS, HUENEG, CS,
     2    SCAL, IDEV, STEM2 )
C
      STOP
      END
C*********************************************************************

C*********************************************************************
C Steve Gibbons 15 January 2004
C If RMFLAG is not equal to 1, return with no action.
C If RMFLAG is equal to 1, then calculate the mean of the real
C array and subtract it from each element.
C We also return with no action if N1 or N2 is less than 1.
C--------------------------------------------------
      SUBROUTINE ARR_RM_MEAN( N1, N2, ARR, RMFLAG )
      IMPLICIT NONE
      INTEGER        RMFLAG, N1, N2
      REAL           ARR( N1, N2 )
C
      INTEGER        I1, I2
      REAL           RTOTEL, RMEANV
C
      IF ( RMFLAG.NE.1 ) RETURN
      IF ( N1.LT.1 )     RETURN
      IF ( N2.LT.1 )     RETURN
C
      RTOTEL = 1.0/FLOAT( N1*N2 )
      RMEANV = 0.0
      DO I2 = 1, N2
        DO I1 = 1, N1
          RMEANV = RMEANV + ARR( I1, I2 )*RTOTEL
        ENDDO
      ENDDO
C
C Have now calculated the mean - now subtract it ...
C
      DO I2 = 1, N2
        DO I1 = 1, N1
          ARR( I1, I2 ) = ARR( I1, I2 ) - RMEANV
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************

