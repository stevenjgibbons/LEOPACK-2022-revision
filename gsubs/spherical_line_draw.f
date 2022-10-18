C*********************************************************************
C subroutine SPHERICAL LINE DRAW *************************************
C            --------- ---- ---- *************************************
C Steve Gibbons Tue Mar 20 11:48:52 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Takes in the 3 by 3 double precision rotation matrix EARCM and     C
C a character string of undetermined length, which contains all      C
C other parameters.                                                  C
C                                                                    C
C*********************************************************************
      SUBROUTINE SPHERICAL_LINE_DRAW( EARCM, LINE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION EARCM( 3, 3 )
      CHARACTER *(*)   LINE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IW, NPTS, ISTYLE, IPTS
      REAL    RADV, THEV, PHIV, THE1, THE2, DTHE,
     1        XWORLD, YWORLD, DEPTH, RAD1, RAD2, DRAD,
     1        PHI1, PHI2, DPHI, PI
      CHARACTER *(3)   CODE
      PARAMETER ( PI = 3.1415927 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      CODE(1:3) = LINE(1:3)
C
C First do case of  code.eq.'CRP' -
C Constant radius and phi
C
      IF ( CODE.EQ.'CRP' ) THEN
C
C Line(4:) must have the variables
C RADV, PHIV, THE1, THE2, IW, NPTS, ISTYLE
C
        READ ( LINE(4:), * ) RADV, PHIV, THE1, THE2, IW, NPTS,
     1                        ISTYLE
C
C RADV   = radial value
C PHIV   = phi value
C THE1   = starting theta value
C THE2   = ending theta value
C IW     = line thickness
C NPTS   = number of points in theta
C ISTYLE = 0: Solid: show uansett
C ISTYLE = 1: Solid: do NOT show if we are "hidden"
C ISTYLE = 2: Dashed: do NOT show if we are "hidden"
C ISTYLE = 3: Solid: ONLY show if we are "hidden"
C ISTYLE = 4: Dashed: ONLY show if we are "hidden"
C
        CALL PGSLW( IW )
C
        DTHE = (THE2 - THE1)/REAL( NPTS )
        DO IPTS = 0, NPTS
          THEV = THE1 + REAL( IPTS )*DTHE
          CALL SPHER_SAT_2_WORLD( RADV, THEV, PHIV, EARCM,
     1                          XWORLD, YWORLD, DEPTH )
          IF ( IPTS.EQ.0 ) THEN
            CALL PGMOVE( XWORLD, YWORLD )
          ELSE
C ..........+                                  istyle = 0
            IF ( ISTYLE.EQ.0 ) THEN
              CALL PGDRAW( XWORLD, YWORLD )
            ENDIF
C ..........+
C ..........+                                  istyle = 1
            IF ( ISTYLE.EQ.1 ) THEN
              IF ( DEPTH.LT.0.0 ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 1
C ..........+                                  istyle = 2
            IF ( ISTYLE.EQ.2 ) THEN
              IF ( DEPTH.LT.0.0 .OR. IPTS/2*2.EQ.IPTS ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 2
C ..........+                                  istyle = 3
            IF ( ISTYLE.EQ.3 ) THEN
              IF ( DEPTH.GE.0.0 ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 3
C ..........+                                  istyle = 4
            IF ( ISTYLE.EQ.4 ) THEN
              IF ( DEPTH.GE.0.0 .OR. IPTS/2*2.EQ.IPTS ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 4
          ENDIF
        ENDDO
C
        RETURN
      ENDIF
C
C Now do case of  code.eq.'CRT' -
C Constant radius and theta
C
      IF ( CODE.EQ.'CRT' ) THEN
C
C Line(4:) must have the variables
C RADV, THEV, PHI1, PHI2, IW, NPTS, ISTYLE
C
        READ ( LINE(4:), * ) RADV, THEV, PHI1, PHI2, IW, NPTS,
     1                        ISTYLE
C
C RADV   = radial value
C THEV   = theta value
C PHI1   = starting phi value
C PHI2   = ending phi value
C IW     = line thickness
C NPTS   = number of points in theta
C ISTYLE = 0: Solid: show uansett
C ISTYLE = 1: Solid: do NOT show if we are "hidden"
C ISTYLE = 2: Dashed: do NOT show if we are "hidden"
C ISTYLE = 3: Solid: ONLY show if we are "hidden"
C ISTYLE = 4: Dashed: ONLY show if we are "hidden"
C
        CALL PGSLW( IW )
C
        DPHI = (PHI2 - PHI1)/REAL( NPTS )
        DO IPTS = 0, NPTS
          PHIV = PHI1 + REAL( IPTS )*DPHI
          CALL SPHER_SAT_2_WORLD( RADV, THEV, PHIV, EARCM,
     1                          XWORLD, YWORLD, DEPTH )
          IF ( IPTS.EQ.0 ) THEN
            CALL PGMOVE( XWORLD, YWORLD )
          ELSE
C ..........+                                  istyle = 0
            IF ( ISTYLE.EQ.0 ) THEN
              CALL PGDRAW( XWORLD, YWORLD )
            ENDIF
C ..........+
C ..........+                                  istyle = 1
            IF ( ISTYLE.EQ.1 ) THEN
              IF ( DEPTH.LT.0.0 ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 1
C ..........+                                  istyle = 2
            IF ( ISTYLE.EQ.2 ) THEN
              IF ( DEPTH.LT.0.0 .OR. IPTS/2*2.EQ.IPTS ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 2
C ..........+                                  istyle = 3
            IF ( ISTYLE.EQ.3 ) THEN
              IF ( DEPTH.GE.0.0 ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 3
C ..........+                                  istyle = 2
            IF ( ISTYLE.EQ.4 ) THEN
              IF ( DEPTH.GE.0.0 .OR. IPTS/2*2.EQ.IPTS ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 4
          ENDIF
        ENDDO
C
        RETURN
      ENDIF
C
C Now do case of  code.eq.'CPT' -
C Constant phi and theta
C
      IF ( CODE.EQ.'CPT' ) THEN
C
C Line(4:) must have the variables
C PHIV, THEV, RAD1, RAD2, IW, NPTS, ISTYLE
C
        READ ( LINE(4:), * ) PHIV, THEV, RAD1, RAD2, IW, NPTS,
     1                        ISTYLE
C
C PHIV   = phi value
C THEV   = theta value
C RAD1   = starting radial value
C RAD2   = ending radial value
C IW     = line thickness
C NPTS   = number of points in theta
C ISTYLE = 0: Solid: show uansett
C ISTYLE = 1: Solid: do NOT show if we are "hidden"
C ISTYLE = 2: Dashed: do NOT show if we are "hidden"
C ISTYLE = 3: Solid: ONLY show if we are "hidden"
C ISTYLE = 4: Dashed: ONLY show if we are "hidden"
C
        CALL PGSLW( IW )
C
        DRAD = (RAD2 - RAD1)/REAL( NPTS )
        DO IPTS = 0, NPTS
          RADV = RAD1 + REAL( IPTS )*DRAD
          CALL SPHER_SAT_2_WORLD( RADV, THEV, PHIV, EARCM,
     1                          XWORLD, YWORLD, DEPTH )
          IF ( IPTS.EQ.0 ) THEN
            CALL PGMOVE( XWORLD, YWORLD )
          ELSE
C ..........+                                  istyle = 0
            IF ( ISTYLE.EQ.0 ) THEN
              CALL PGDRAW( XWORLD, YWORLD )
            ENDIF
C ..........+
C ..........+                                  istyle = 1
            IF ( ISTYLE.EQ.1 ) THEN
              IF ( DEPTH.LT.0.0 ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 1
C ..........+                                  istyle = 2
            IF ( ISTYLE.EQ.2 ) THEN
              IF ( DEPTH.LT.0.0 .OR. IPTS/2*2.EQ.IPTS ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 2
C ..........+                                  istyle = 3
            IF ( ISTYLE.EQ.3 ) THEN
              IF ( DEPTH.GE.0.0 ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 3
C ..........+                                  istyle = 4
            IF ( ISTYLE.EQ.4 ) THEN
              IF ( DEPTH.GE.0.0 .OR. IPTS/2*2.EQ.IPTS ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 4
          ENDIF
        ENDDO
C
        RETURN
      ENDIF
C
C Now do case of  code.eq.'CRD' -
C Circle Ring Draw
C
      IF ( CODE.EQ.'CRD' ) THEN
C
C Line(4:) must have the variables
C RADV, IW, NPTS, ISTYLE
C ISTYLE = 0: Solid: show uansett
C ISTYLE = 1 --> solid
C ISTYLE = 2 --> dashed
C
        READ ( LINE(4:), * ) RADV, IW, NPTS, ISTYLE
C
        CALL PGSLW( IW )
C
        DTHE = 2.0*PI/REAL( NPTS )
        DO IPTS = 0, NPTS
          THEV = THE1 + REAL( IPTS )*DTHE
          XWORLD = RADV*COS( THEV )
          YWORLD = RADV*SIN( THEV )
          IF ( IPTS.EQ.0 ) THEN
            CALL PGMOVE( XWORLD, YWORLD )
          ELSE
C ..........+                                  istyle = 0
            IF ( ISTYLE.EQ.0 ) THEN
              CALL PGDRAW( XWORLD, YWORLD )
            ENDIF
C ..........+
C ..........+                                  istyle = 1
            IF ( ISTYLE.EQ.1 ) THEN
              CALL PGDRAW( XWORLD, YWORLD )
            ENDIF
C ..........+                                  istyle = 1
C ..........+                                  istyle = 2
            IF ( ISTYLE.EQ.2 ) THEN
              IF ( IPTS/2*2.EQ.IPTS ) THEN
                CALL PGMOVE( XWORLD, YWORLD )
              ELSE
                CALL PGDRAW( XWORLD, YWORLD )
              ENDIF
            ENDIF
C ..........+                                  istyle = 2
          ENDIF
        ENDDO
C
        RETURN
      ENDIF
C
      PRINT *,' Subroutine SPHERICAL_LINE_DRAW.'
      PRINT *,' CODE = ', CODE
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
