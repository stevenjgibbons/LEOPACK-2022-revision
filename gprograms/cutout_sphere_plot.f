C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C Tue Oct  9 16:15:12 WEST 2001                                      C
C                                                                    C
C General Cut Out Sphere Plot                                        C
C                                                                    C
C*********************************************************************
      PROGRAM general_cutout_sphere_plot
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
      INTEGER NRADMX, NTHMAX, NPHMAX, NLEVM, LHMAX, NHMAX,
     1        ISVMAX, NNDM, NPMAX
      PARAMETER ( NRADMX = 250, NTHMAX = 250, NPHMAX = 250,
     1            NLEVM = 20, LHMAX = 160, NHMAX = 3000,
     2            ISVMAX = NRADMX*NHMAX, NNDM = 6,
     3            NPMAX = (LHMAX+1)*(LHMAX+2)/2 )
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
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     2        IWORK( NNDM )
C
      DOUBLE PRECISION XARR( NRADMX ), VEC( ISVMAX ),
     1                 WORK1( NNDM ), WORK2( NNDM ),
     2                 COEFM( NNDM, NNDM ), WORK3( NTHMAX )
C
      DOUBLE PRECISION PA( NPMAX, NTHMAX ), DPA( NPMAX, NTHMAX )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C     .
C     . General
C     .
      INTEGER          NLEV, IW, LUIN, LU, IOP, IAS, NNDS,
     1                 IDEV, INARR( 3 ), NH, ICOMP,
     2                 LH, NR, NR1, I, ILEN,
     3                 NTHP, NTHEP, NPHI, NPHIP
      CHARACTER *(1)   ESCPCH
      CHARACTER *(180) STEM, LINE, FNAME, STEM2
      REAL             PI, ZERO, DPHI, DTHE, DTHETA, PAPWIDTH,
     1                 PHIV1, PHIV2, RADV1, RADV2
      DOUBLE PRECISION PHI, VMISS, THE, RAD
      DOUBLE PRECISION ALPHAD, BETAD, GAMMAD, EARCM( 3, 3 )
C     .
C Variables for setting the colour.
C
      REAL    HUEPOS, HUENEG, CS, SCAL
C     .
      PARAMETER ( PI=3.1415926, ZERO = 0.0, VMISS = 9.9d9, IAS = 0 )
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
      PRINT *,' Enter huepos, hueneg, csat, scal, iw'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) HUEPOS, HUENEG, CS, SCAL, IW
C
      RFIRST = REAL( XARR( 1 ) )
      RADV1  = RFIRST
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
      IF ( NLEV.GT.NLEVM ) THEN
        PRINT *,' NLEV  = ', NLEV 
        PRINT *,' NLEVM = ', NLEVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Enter resolutions and positions
C
      PRINT *,' Enter NRAD, NTHP, NTHEP, NPHI, NPHIP, PHIV1, RADV2'
      CALL LINERD( LUIN, LINE, ESCPCH )
      READ ( LINE, * ) NRAD, NTHP, NTHEP, NPHI, NPHIP, PHIV1, RADV2
C
      RLAST = RADV2
C
      IF ( NRAD.GT.NRADMX .OR. NTHP.GT.NTHMAX .OR.
     1     NTHEP.GT.NTHP .OR. NPHI.GT.NPHMAX .OR.
     2     NPHIP.GE.NPHI ) THEN
        PRINT *,' NRADMX = ', NRADMX
        PRINT *,' NTHMAX = ', NTHMAX
        PRINT *,' NPHMAX = ', NPHMAX
        PRINT *,' NRAD   = ', NRAD
        PRINT *,' NTHP   = ', NTHP
        PRINT *,' NTHEP  = ', NTHEP
        PRINT *,' NPHI   = ', NPHI
        PRINT *,' NPHIP  = ', NPHIP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DTHE   = PI/FLOAT( NTHP - 1 )
      DTHETA = DTHE*FLOAT( NTHEP )
      DPHI   = 2.0*PI/FLOAT( NPHI - 1 )
      PHIV2  = PHIV1 + FLOAT( NPHIP )*DPHI
C
C Fill data arrays
C
      IOP = 0
      CALL RMATOP( FUNCRAD, ZERO, NPHMAX, NTHMAX, IOP )
      CALL RMATOP( FUNCP1, ZERO, NRADMX, NTHMAX, IOP )
      CALL RMATOP( FUNCP2, ZERO, NRADMX, NTHMAX, IOP )
      CALL RMATOP( FUNCTHE, ZERO, NRADMX, NPHMAX, IOP )
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
C Now enter the constant phi section 1
C
      NTHE   = NTHEP
      PHI    = DBLE( PHIV1 )
      RFIRST = RADV1
      RLAST  = RADV2
      TFIRST = 0.4999*PI
      TLAST  = 0.5*PI - DTHETA
      CALL MERID_SEC_POLAR_EVAL( ICOMP, INARR, NNDS, IWORK, MHT,
     1           MHL, MHM, LH, IAS, PHI, VEC, XARR, WORK1, WORK2,
     2           COEFM, WORK3, PA, DPA, VMISS, FUNCP1 )
C
C Now enter the constant phi section 2
C
      NTHE   = NTHEP
      PHI    = DBLE( PHIV2 )
      RFIRST = RADV1
      RLAST  = RADV2
      TFIRST = 0.4999*PI
      TLAST  = 0.5*PI - DTHETA
      CALL MERID_SEC_POLAR_EVAL( ICOMP, INARR, NNDS, IWORK, MHT,
     1           MHL, MHM, LH, IAS, PHI, VEC, XARR, WORK1, WORK2,
     2           COEFM, WORK3, PA, DPA, VMISS, FUNCP2 )
C
C Now enter the constant theta section
C
      IF ( NTHEP.LT.NTHP ) THEN
        NTHE   = NPHIP
        THE    = DBLE( DTHETA )
        RFIRST = RADV1
        RLAST  = RADV2
        TFIRST = PHIV1
        TLAST  = PHIV2
        CALL EQ_SEC_POLAR_EVAL( ICOMP, INARR, NNDS, IWORK, MHT,
     1           MHL, MHM, LH, IAS, THE, VEC, XARR, WORK1, WORK2,
     2           COEFM, VMISS, FUNCTHE )
      ENDIF
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
      CALL CUT_OUT_SPHERE_PLOT( NLEV, NPHI, NTHP, FUNCRAD,
     1    RADV1, RADV2, NRAD, NTHEP, FUNCP1, PHIV1, NRAD, NTHEP,
     2    FUNCP2, PHIV2, NRAD, NPHIP, FUNCTHE, DTHETA, NLEVM, CONTL,
     3    PAPWIDTH, EARCM, HUEPOS, HUENEG, CS, SCAL, IDEV, STEM2 )
C
      STOP
      END
C*********************************************************************
