C*********************************************************************
C subroutine ARRAY SPHERical PROJection COLOUR ***********************
C            ----- -----     ----       ------ ***********************
C Steve Gibbons Tue Mar 20 07:27:37 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Colours in between (undrawn) contours for an array in polar coords C
C It automatically calculates appropriate contours etc.              C
C                                                                    C
C If ICONTOUR = 1, then the contours are scaled between the          C
C minimum and maximum values of the array (adjusted by DELTA)        C
C                                                                    C
C If ICONTOUR = 2, then the contours are scaled between the          C
C imposed minimum and maximum values VALMIN and VALMAX.              C
C ARRAY_GENERAL_COLOUR doesn't tend to trust humans very much        C
C and works out what the min and max are. He will abort the program  C
C if the imposed limit is too small!                                 C
C                                                                    C
C Depending upon the                                                 C
C integer flag, IPFLAG, the spherical coordinates which describe     C
C the elements of A can be as follows:-                              C
C                                                                    C
C  IPFLAG = 1:                                                       C
C  -----------                                                       C
C                                                                    C
C    This is a meridian section with PHI (in radians) = COORD.       C
C                                                                    C
C    A( I, J ) is a function evaluated at RAD, THE, PHI              C
C                                                                    C
C    where RAD = C1V1 + REAL( I - 1 )*DRAD                           C
C               with DRAD = (C1V2 - C1V1)/REAL( IDIM - 1 )           C
C  and                                                               C
C          THE = C2V1 + REAL( J - 1 )*DTHE                           C
C               with DTHE = (C2V2 - C2V1)/REAL( JDIM - 1 )           C
C  THE is in radians.                                                C
C                                                                    C
C  IPFLAG = 2:                                                       C
C  -----------                                                       C
C                                                                    C
C    This is a constant theta section with THE (in radians) = COORD. C
C                                                                    C
C    A( I, J ) is a function evaluated at RAD, THE, PHI              C
C                                                                    C
C    where RAD = C1V1 + REAL( I - 1 )*DRAD                           C
C               with DRAD = (C1V2 - C1V1)/REAL( IDIM - 1 )           C
C  and                                                               C
C          PHI = C2V1 + REAL( J - 1 )*DPHI                           C
C               with DPHI = (C2V2 - C2V1)/REAL( JDIM - 1 )           C
C  PHI is in radians.                                                C
C                                                                    C
C  IPFLAG = 3:                                                       C
C  -----------                                                       C
C                                                                    C
C    This is a constant radius section with RAD = COORD.             C
C                                                                    C
C    A( I, J ) is a function evaluated at RAD, THE, PHI              C
C                                                                    C
C    where PHI = C1V1 + REAL( I - 1 )*DPHI                           C
C               with DRAD = (C1V2 - C1V1)/REAL( IDIM - 1 )           C
C  and                                                               C
C          THE = C2V1 + REAL( J - 1 )*DTHE                           C
C               with DTHE = (C2V2 - C2V1)/REAL( JDIM - 1 )           C
C  PHI and THE are in radians.                                       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NLEV      : Number of contour levels.                          C
C                  NLEV = -1 --> 17 contour levels with Andy's RBG   C
C                  colour scheme.                                    C
C     IDIM      : Number of radial grid nodes.                       C
C     JDIM      : Number of theta grid nodes.                        C
C                                                                    C
C     ICONTOUR  : 1 --> scale contours based on local values         C
C                 2 --> scale contours based on imposed values.      C
C                                                                    C
C     IPFLAG    : See above.                                         C
C                                                                    C
C     IDEP      : = 0 --> contour all regions regardless of depth    C
C                 = 1 --> only contour with depth.ge.0.0             C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Data array with dimensions (IDIM,JDIM)             C
C     C1V1      : See above.                                         C
C     C1V2      : See above.                                         C
C     C2V1      : See above.                                         C
C     C2V2      : See above.                                         C
C     CONTL     : Dim (NLEV). Work array (to store contour values)   C
C     HUEPOS    : Colour for positive values                         C
C     HUENEG    : Colour for negative values                         C
C                                                                    C
C huepos and hueneg can take values in the interval [0,360] degrees. C
C 0 --> blue, 120 --> red and 240 --> green.                         C
C                                                                    C
C     CS        : 0.0 for monochrome and up to 1.0 for colour.       C
C     SCAL      : In interval ( 0.0, 1.0] 1.0 gives darkest colours  C
C                                                                    C
C   HUEPOS, HUENEG, CS and SCAL are only used if NLEV.NE.-1          C
C                                                                    C
C     DELTA     : Difference between highest contour level and       C
C                 maximum array value, or minimum array value and    C
C                 lowest contour level.                              C
C                                                                    C
C     VALMIN    : Imposed minimum value.                             C
C     VALMAX    : Imposed maximum value.                             C
C                                                                    C
C (Note that valmin and valmax are both ignored if ICONTOUR = 1)     C
C                                                                    C
C     COORD     : See above.                                         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     EARCM     : Euler angle matrix Dim( 3, 3 ). (See EARCMC)       C
C                                                                    C
C*********************************************************************
      SUBROUTINE ARRAY_SPHER_PROJ__COLOUR( F, C1V1, C1V2, C2V1, C2V2,
     1                IDIM, JDIM, NLEV, CONTL, HUEPOS, HUENEG, CS,
     2                SCAL, DELTA, ICONTOUR, VALMIN, VALMAX,
     3                IPFLAG, EARCM, COORD, IDEP )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IDIM, JDIM, NLEV, ICONTOUR, IPFLAG, IDEP
      REAL    F( IDIM, JDIM ), C1V1, C1V2, C2V1, C2V2,
     1        CONTL( NLEV ), HUEPOS, HUENEG, CS, SCAL, DELTA,
     2        VALMIN, VALMAX, COORD
      DOUBLE PRECISION EARCM( 3, 3 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I1, J1, I2, J2, ILEV, NLEV2
      REAL    RRED( 17 ), RRED1,
     1        RBLUE( 17 ), RBLUE1,
     2        RGREEN( 17 ), RGREEN1
      REAL    FMIN, FMAX, ZERO, FLIMIT, C1, C2, FEXTRM
      PARAMETER ( ZERO = 0.0, FLIMIT = 1.0e8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NLEV.EQ.-1 ) THEN
        NLEV2 = 17
      ELSE
        NLEV2 = NLEV
      ENDIF
C
      IF ( ICONTOUR.NE.1 .AND. ICONTOUR.NE.2 ) THEN
        PRINT *,' Subroutine ARRAY_SPHER_PROJ_COLOUR.'
        PRINT *,' ICONTOUR = ', ICONTOUR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DELTA.LT.ZERO ) THEN
        PRINT *,' Subroutine ARRAY_SPHER_PROJ_COLOUR.'
        PRINT *,' DELTA = ', DELTA
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Find FMIN and FMAX
C
      FMAX = FLIMIT*(-1.0)
      FMIN = FLIMIT
      DO J1 = 1, JDIM
        DO I1 = 1, IDIM
          IF ( F( I1, J1 ).GT.FMAX .AND.
     1         F( I1, J1 ).LT.FLIMIT ) FMAX = F( I1, J1 )
          IF ( F( I1, J1 ).LT.FMIN ) FMIN = F( I1, J1 )
        ENDDO
      ENDDO
C
      IF ( ICONTOUR.EQ.2 ) THEN
        IF ( VALMAX.LT.FMAX .OR. VALMIN.GT.FMIN ) THEN
          PRINT *,' Subroutine ARRAY_SPHER_PROJ_COLOUR.'
          PRINT *,' ICONTOUR = 2.'
          PRINT *,' VALMIN   = ', VALMIN
          PRINT *,'   FMIN   = ',   FMIN
          PRINT *,' VALMAX   = ', VALMAX
          PRINT *,'   FMAX   = ',   FMAX
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        FMIN = VALMIN
        FMAX = VALMAX
      ENDIF
C
C Modify FMIN and FMAX
C
      IF ( ICONTOUR.EQ.1 ) THEN
        FMIN = FMIN - DELTA
        FMAX = FMAX + DELTA
      ENDIF
C
C Define contour levels
C
      CALL RCLDR( NLEV2, FMIN, FMAX, CONTL )
C
C Select the colour scheme
C
      IF ( NLEV.EQ.-1 ) THEN
C
        FEXTRM = ABS( FMIN )
        IF ( ABS( FMAX ).GT.FEXTRM ) FEXTRM = ABS( FMAX )
        FMIN = (-1.0)*FEXTRM
        FMAX = FEXTRM
        CALL RCLDR( NLEV2, FMIN, FMAX, CONTL )
C
        RRED(   2 ) = 0.619608
        RRED(   3 ) = 0.666667
        RRED(   4 ) = 0.749020
        RRED(   5 ) = 0.819608
        RRED(   6 ) = 0.894118
        RRED(   7 ) = 1.000000
        RRED(   8 ) = 1.000000
        RRED(   9 ) = 1.000000
        RRED(  10 ) = 0.749020
        RRED(  11 ) = 0.647059
        RRED(  12 ) = 0.564706
        RRED(  13 ) = 0.466667
        RRED(  14 ) = 0.384314
        RRED(  15 ) = 0.329412
        RRED(  16 ) = 0.200000
        RRED(  17 ) = 0.094118
C
        RGREEN(   2 ) =  0.039216
        RGREEN(   3 ) =  0.156863
        RGREEN(   4 ) =  0.247059
        RGREEN(   5 ) =  0.349020
        RGREEN(   6 ) =  0.447059
        RGREEN(   7 ) =  0.549020
        RGREEN(   8 ) =  0.647059
        RGREEN(   9 ) =  0.756863
        RGREEN(  10 ) =  0.749020
        RGREEN(  11 ) =  0.647059
        RGREEN(  12 ) =  0.564706
        RGREEN(  13 ) =  0.466667
        RGREEN(  14 ) =  0.384314
        RGREEN(  15 ) =  0.329412
        RGREEN(  16 ) =  0.200000
        RGREEN(  17 ) =  0.188235
C
        RBLUE(   2 ) =  0.000000
        RBLUE(   3 ) =  0.039216
        RBLUE(   4 ) =  0.039216
        RBLUE(   5 ) =  0.039216
        RBLUE(   6 ) =  0.039216
        RBLUE(   7 ) =  0.156863
        RBLUE(   8 ) =  0.313726
        RBLUE(   9 ) =  0.466667
        RBLUE(  10 ) =  1.000000
        RBLUE(  11 ) =  1.000000
        RBLUE(  12 ) =  1.000000
        RBLUE(  13 ) =  1.000000
        RBLUE(  14 ) =  0.909804
        RBLUE(  15 ) =  0.839216
        RBLUE(  16 ) =  0.776471
        RBLUE(  17 ) =  0.713726
C
        DO ILEV = 2, 17
          RRED1 = RRED( ILEV )
          RBLUE1 = RBLUE( ILEV )
          RGREEN1 = RGREEN( ILEV )
          CALL PGSCR( ILEV, RRED1, RGREEN1, RBLUE1 )
        ENDDO
      ELSE
        CALL HLS_CONT_SET( NLEV2, CONTL, HUEPOS, HUENEG, CS, SCAL )
      ENDIF
C
C Fill in contour levels
C
      I1 = 1
      I2 = IDIM
C
      J1 = 1
      J2 = JDIM
C
      DO ILEV = 2, NLEV2
C
C Sets the colour for the interval between contours
C ILEV - 1 and ILEV.
C
        CALL PGSCI( ILEV )
C
        C1  = CONTL( ILEV-1 )
        C2  = CONTL( ILEV )
C
        CALL PGCONF_SPHER_PROJ( F, C1V1, C1V2, C2V1, C2V2, IDIM,
     1                          JDIM, I1, I2, J1, J2, C1, C2,
     2                          IPFLAG, EARCM, COORD, IDEP )
C
      ENDDO
C
      RETURN
      END
C*********************************************************************
