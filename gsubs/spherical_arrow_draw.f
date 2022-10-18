C*********************************************************************
C subroutine SPHERICAL ARROW DRAW ************************************
C            --------- ----- ---- ************************************
C Steve Gibbons Tue Mar 20 11:48:52 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Takes in the 3 by 3 double precision rotation matrix EARCM and     C
C a character string of undetermined length, which contains all      C
C other parameters.                                                  C
C                                                                    C
C*********************************************************************
      SUBROUTINE SPHERICAL_ARROW_DRAW( EARCM, LINE )
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
C Now do case of  code.eq.'CRD' -
C Constant phi arrow
C
      IF ( CODE.EQ.'CPA' ) THEN
C
C Line(4:) must have the variables
C RADV, IW, NPTS, ISTYLE
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
      PRINT *,' Subroutine SPHERICAL_LINE.'
      PRINT *,' CODE = ', CODE
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
