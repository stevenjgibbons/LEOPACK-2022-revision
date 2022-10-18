C*********************************************************************
C subroutine ARRAY GENERAL COLOUR ************************************
C            ----- ------- ------ ************************************
C Steve Gibbons Mon Mar 19 13:44:14 MET 2001                         C
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
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NLEV      : Number of contour levels.                          C
C     IDIM      : Number of radial grid nodes.                       C
C     JDIM      : Number of theta grid nodes.                        C
C                                                                    C
C     ICONTOUR  : 1 --> scale contours based on local values         C
C                 2 --> scale contours based on imposed values.      C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Data array with dimensions (IDIM,JDIM)             C
C     XWA       : Corresponding x world coordinates (IDIM,JDIM)      C
C     YWA       : Corresponding y world coordinates (IDIM,JDIM)      C
C     CONTL     : Dim (NLEV). Work array (to store contour values)   C
C     HUEPOS    : Colour for positive values                         C
C     HUENEG    : Colour for negative values                         C
C                                                                    C
C huepos and hueneg can take values in the interval [0,360] degrees. C
C 0 --> blue, 120 --> red and 240 --> green.                         C
C                                                                    C
C     CS        : 0.0 for monochrome and up to 1.0 for colour.       C
C     SCAL      : In interval ( 0.0, 1.0] 1.0 gives darkest colours  C
C     DELTA     : Difference between highest contour level and       C
C                 maximum array value, or minimum array value and    C
C                 lowest contour level.                              C
C                                                                    C
C     VALMIN    : Imposed minimum value.                             C
C     VALMAX    : Imposed maximum value.                             C
C                                                                    C
C (Note that valmin and valmax are both ignored if ICONTOUR = 1)     C
C                                                                    C
C*********************************************************************
      SUBROUTINE ARRAY_GENERAL_COLOUR( F, XWA, YWA, IDIM, JDIM, NLEV,
     1                       CONTL, HUEPOS, HUENEG, CS, SCAL, DELTA,
     2                       ICONTOUR, VALMIN, VALMAX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IDIM, JDIM, NLEV, ICONTOUR
      REAL    F( IDIM, JDIM ), XWA( IDIM, JDIM ), YWA( IDIM, JDIM ),
     1        CONTL( NLEV ), HUEPOS, HUENEG, CS, SCAL, DELTA,
     2        VALMIN, VALMAX
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I1, J1, I2, J2, ILEV
      REAL    FMIN, FMAX, ZERO, FLIMIT, C1, C2
      PARAMETER ( ZERO = 0.0, FLIMIT = 1.0e8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( ICONTOUR.NE.1 .AND. ICONTOUR.NE.2 ) THEN
        PRINT *,' Subroutine ARRAY_POLAR_COLOUR_OPTION.'
        PRINT *,' ICONTOUR = ', ICONTOUR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DELTA.LT.ZERO ) THEN
        PRINT *,' Subroutine ARRAY_POLAR_COLOUR_OPTION.'
        PRINT *,' DELTA = ', DELTA
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Find FMIN and FMAX
C
      FMAX = F( 8, 8 )
      FMIN = FMAX
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
          PRINT *,' Subroutine ARRAY_POLAR_COLOUR_OPTION.'
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
      FMIN = FMIN - DELTA
      FMAX = FMAX + DELTA
C
C Define contour levels
C
      CALL RCLDR( NLEV, FMIN, FMAX, CONTL )
C
C Select the colour scheme
C
      CALL HLS_CONT_SET( NLEV, CONTL, HUEPOS, HUENEG, CS, SCAL )
C
C Fill in contour levels
C
      I1 = 1
      I2 = IDIM
C
      J1 = 1
      J2 = JDIM
C
      DO ILEV = 2, NLEV
C
C Sets the colour for the interval between contours
C ILEV - 1 and ILEV.
C
        CALL PGSCI( ILEV )
C
        C1  = CONTL( ILEV-1 )
        C2  = CONTL( ILEV )
C
        CALL PGCONF_GENERAL( F, XWA, YWA, IDIM, JDIM, I1, I2, J1,
     1                          J2, C1, C2)
      ENDDO
C
      RETURN
      END
C*********************************************************************
