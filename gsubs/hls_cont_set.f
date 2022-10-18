C*********************************************************************
C subroutine Hue Lightness Saturation CONTour SET ********************
C            -   -         -          ----    -   ********************
C Steve Gibbons Fri Jan 19 10:02:14 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C You have already defined an array CONTL with NLEV real values.     C
C This gives you (NLEV-1) intervals. When the contour map is filled  C
C with colours, the interval between the contour lines (ILEV-1) and  C
C ILEV is filled with the colour corresponding to the integer ILEV.  C
C                                                                    C
C HLS_CONT_SET, defines the colour for that level.                   C
C The hue-light-saturation scheme has three parameters:              C
C   CH (hue), CL(lightness) and CS(saturation).                      C
C                                                                    C
C CH must be in the interval [0,360] (measured in degrees)           C
C                                                                    C
C  CH =   0  --> blue                                                C
C  CH = 120  --> red                                                 C
C  CH = 240  --> green                                               C
C                                                                    C
C CL = 1 --> brightest  :  CL must take values in the                C
C CL = 0 --> darkest    :  interval [0, 1].                          C
C                                                                    C
C CS = 1 --> colour     :  CS must take values in the                C
C CS = 0 --> monochrome :  interval [0, 1].                          C
C                                                                    C
C If CS = 0 then the CH value is irrelevant.                         C
C                                                                    C
C If CS = 0, then HLS_CONT_SET has the minimum value white and       C
C the maximum value black, regardless of values.                     C
C                                                                    C
C If CS is greater than zero, HLS_CONT_SET grades the brightness     C
C with negative values in the colour HUENEG and positive values      C
C in the colour HUEPOS.                                              C
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
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     CONTL     : Dim (NLEV). Resulting contour levels.              C
C     HUEPOS    : Colour for positive values                         C
C     HUENEG    : Colour for negative values                         C
C     CS        : 0.0 for monochrome and up to 1.0 for colour.       C
C     SCAL      : In interval ( 0.0, 1.0] 1.0 gives darkest colours  C
C                                                                    C
C*********************************************************************
      SUBROUTINE HLS_CONT_SET( NLEV, CONTL, HUEPOS, HUENEG, CS, SCAL)
      IMPLICIT NONE
      INTEGER NLEV
      REAL    CONTL( NLEV ), HUEPOS, HUENEG, CS, SCAL
C-------------------------------------------
      INTEGER     ILEV, NINTERVALS
      REAL        FMIN, FMAX, CH, CL, ZERO,
     1            RINT, DIFF, CMID, FAC
      PARAMETER ( ZERO = 0.0 )
C-------------------------------------------
      FMIN    =  CONTL(    1 )
      FMAX    =  CONTL( NLEV )
      IF ( FMIN.GE.ZERO .OR. FMAX.LE.ZERO .OR. CS.EQ.ZERO ) THEN
C       .
C       . So our contours are either all of the same sign
C       . are it is a monochrome plot.
C       .
        DIFF       = FMAX - FMIN
        NINTERVALS = NLEV - 1
C       .
      ELSE
C       .
C       . So our contours are of different signs
C       .
        DIFF = -FMIN
        IF ( FMAX.GT.DIFF ) DIFF = FMAX
        NINTERVALS = 0
        IF ( FMAX.GT.ABS( FMIN ) ) THEN
          DO ILEV = 1, NLEV
            IF ( CONTL( ILEV ).GT.ZERO )
     1                NINTERVALS = NINTERVALS + 1
          ENDDO
        ELSE
          DO ILEV = 1, NLEV
            IF ( CONTL( ILEV ).LT.ZERO )
     1                NINTERVALS = NINTERVALS + 1
          ENDDO
        ENDIF
C       .
      ENDIF
      RINT = 1.0/REAL( NINTERVALS )
      FAC  = RINT
C-------------------------------------
      DO ILEV = 2, NLEV
        CMID = 0.5*( CONTL(ILEV-1) + CONTL(ILEV) )
        IF ( CMID.LT.ZERO ) THEN
          CH = HUENEG
        ELSE
          CH = HUEPOS
        ENDIF
        IF ( FMAX*FMIN.GT.ZERO .OR. CS.EQ.ZERO ) THEN
           CL = (1.0 - SCAL*ABS(CMID-FMIN)/DIFF)
        ELSE
           CL = (1.0 - SCAL*ABS(CMID)/DIFF)
        ENDIF
C-------------------------------------
C Note that CL = 1  -->  brightest   |
C Note that CL = 0  -->  darkest     |
C-------------------------------------
c     print *,' ilev= ',ilev,' cmid = ',cmid,' ch= ',ch,' cl= ',cl
C-------------------------------------
        CALL PGSHLS( ILEV, CH, CL, CS )
      ENDDO
C-------------------------------------------
      RETURN
      END
C*********************************************************************
