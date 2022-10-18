C*********************************************************************
C subroutine ARRAY RECT DRAW (with OPTION) MODified ******************
C            ----- ---- ----       ------  ---      ******************
C Steve Gibbons Tue Feb 27 11:31:44 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Draws contour lines for an array in Cartesian coords.              C
C It automatically calculates appropriate contours etc.              C
C                                                                    C
C If ICONTOUR = 1, then the contours are scaled between the          C
C minimum and maximum values of the array (adjusted by DELTA)        C
C                                                                    C
C If ICONTOUR = 2, then the contours are scaled between the          C
C imposed minimum and maximum values VALMIN and VALMAX.              C
C ARRAY_RECT_DRAW_OPTION doesn't tend to trust humans very much      C
C and works out what the min and max are. He will abort the program  C
C if the imposed limit is too small!                                 C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NLEV      : Number of contour levels.                          C
C                 (NLEV.EQ.-1 is alias for 17 for colouring)         C
C     NPHI      : Number of phi grid nodes.                          C
C     NTHE      : Number of theta grid nodes.                        C
C     INW       : Width of lines. (Negative contours)                C
C     IPW       : Width of lines. (Positive contours)                C
C     INS       : Style of lines. (Negative contours)                C
C     IPS       : Style of lines. (Positive contours)                C
C     INDNEG    : Colour index for negative contours                 C
C     INDPOS    : Colour index for positive contours                 C
C                                                                    C
C     ICONTOUR  : 1 --> scale contours based on local values         C
C                 2 --> scale contours based on imposed values.      C
C                                                                    C
C  Alternatives for INS, IPS are:-                                   C
C                                                                    C
C       1: full line                                                 C
C       2: long dashes                                               C
C       3: dash-dot-dash-dot                                         C
C       4: dotted                                                    C
C       5: das-dot-dot-dot                                           C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     F         : Data array with dimensions (NPHI,NTHE)             C
C     CONTL     : Dim (NLEV). Work array (to store contour values)   C
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
      SUBROUTINE ARRAY_RECT_DRAW_OPTION_MOD( F, NPHI, NTHE, NLEV,
     1                    CONTL, DELTA, ICONTOUR, VALMIN, VALMAX,
     2                    INW, IPW, INS, IPS, INDNEG, INDPOS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHI, NTHE, NLEV, ICONTOUR,
     1        INW, IPW, INS, IPS, INDNEG, INDPOS
      REAL    F( NPHI, NTHE ), CONTL( NLEV ), DELTA, VALMIN, VALMAX
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHE, IPHI, NLEV2
      REAL    FMIN, FMAX, ZERO, FLIMIT, FEXTRM
      PARAMETER ( ZERO = 0.0, FLIMIT = 1.0e8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NLEV2 = NLEV
      IF ( NLEV.EQ.-1 ) NLEV2 = 17
C
      IF ( ICONTOUR.NE.1 .AND. ICONTOUR.NE.2 ) THEN
        PRINT *,' Subroutine ARRAY_RECT_DRAW_OPTION_MOD.'
        PRINT *,' ICONTOUR = ', ICONTOUR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DELTA.LT.ZERO ) THEN
        PRINT *,' Subroutine ARRAY_RECT_DRAW_OPTION_MOD .'
        PRINT *,' DELTA = ', DELTA
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Find FMIN and FMAX
C
      FMAX = F( 8, 8 )
      FMIN = FMAX
      DO ITHE = 1, NTHE
        DO IPHI = 1, NPHI
          IF ( F( IPHI, ITHE ).GT.FMAX .AND.
     1         F( IPHI, ITHE ).LT.FLIMIT ) FMAX = F( IPHI, ITHE )
          IF ( F( IPHI, ITHE ).LT.FMIN ) FMIN = F( IPHI, ITHE )
        ENDDO
      ENDDO
      PRINT *,'INFO: fmin = ', FMIN,' fmax = ', FMAX
C
      IF ( ICONTOUR.EQ.2 ) THEN
        IF ( VALMAX.LT.FMAX .OR. VALMIN.GT.FMIN ) THEN
          PRINT *,' Subroutine ARRAY_RECT_DRAW_OPTION_MOD.'
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
      IF ( NLEV.EQ.-1 ) THEN
        FEXTRM = ABS( FMIN )
        IF ( ABS( FMAX ).GT.FEXTRM ) FEXTRM = ABS( FMAX )
        FMIN = (-1.0)*FEXTRM
        FMAX = FEXTRM
        CALL RCLDR( NLEV2, FMIN, FMAX, CONTL )
      ENDIF
C
C Fill in contour intervals
C
      CALL CONT_DRAW_RECT_MOD( NPHI, NTHE, F, NLEV2, CONTL,
     1                         INW, IPW, INS, IPS, INDNEG, INDPOS )
C
      RETURN
      END
C*********************************************************************
