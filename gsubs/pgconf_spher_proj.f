C*********************************************************************
C subroutine PGCONF_GENERAL: Adapted from PGCONF *********************
C            --------------                      *********************
C Steve Gibbons Tue Mar 20 07:27:37 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C This subroutine is a very simple adaption from the PGPLOT routine  C
C PGCONF. However, I have made some changes. Depending upon the      C
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
C Set IDEP to 0 to ignore the depth, and to 1 to ignore values with  C
C  DEPTH .lt. 0.0                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C See argument list below.                                           C
C                                                                    C
C*********************************************************************
      SUBROUTINE PGCONF_SPHER_PROJ( A, C1V1, C1V2, C2V1, C2V2, IDIM,
     1                              JDIM, I1, I2, J1, J2, C1, C2,
     2                              IPFLAG, EARCM, COORD, IDEP )
      IMPLICIT NONE
C
      INTEGER          IDIM, JDIM, I1, I2, J1, J2, IPFLAG, IDEP
      REAL             A(IDIM,JDIM), C1V1, C1V2, C2V1, C2V2, C1, C2,
     1                 COORD
      DOUBLE PRECISION EARCM( 3, 3 )
C____________________________________________________________________C
C
C Shade the region between two contour levels of a function defined on
C the nodes of a rectangular grid. The routine uses the current fill
C attributes, hatching style (if appropriate), and color index.
C
C If you want to both shade between contours and draw the contour
C lines, call this routine first (once for each pair of levels) and
C then CALL PGCONT (or PGCONS) to draw the contour lines on top of the
C shading.
C
C Note 1: This routine is not very efficient: it generates a polygon
C fill command for each cell of the mesh that intersects the desired
C area, rather than consolidating adjacent cells into a single polygon.
C
C Note 2: If both contours intersect all four edges of a particular
C mesh cell, the program behaves badly and may consider some parts
C of the cell to lie in more than one contour range.
C
C Note 3: If a contour crosses all four edges of a cell, this
C routine may not generate the same contours as PGCONT or PGCONS
C (these two routines may not agree either). Such cases are always
C ambiguous and the routines use different approaches to resolving
C the ambiguity.
C
C Arguments:
C  A      (input)  : data array.
C  IDIM   (input)  : first dimension of A.
C  JDIM   (input)  : second dimension of A.
C  I1,I2  (input)  : range of first index to be contoured (inclusive).
C  J1,J2  (input)  : range of second index to be contoured (inclusive).
C  C1, C2 (input)  : contour levels; note that C1 must be less than C2.
C--
C 03-Oct-1996 - new routine [TJP].
C-----------------------------------------------------------------------
      INTEGER  I, J, IC, NPT, LEV
c     LOGICAL  PGNOTO
      REAL     DVAL(5), X(8), Y(8), DELTA, XX, YY, C, R,
     1         XVAL(5), YVAL( 5 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      REAL    FLIMIT, XWORLD, YWORLD, DEPTH, RDX, RDY,
     1        THE, PHI, RAD
      PARAMETER ( FLIMIT = 1.0e8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check arguments.
C
      IF ( IPFLAG.NE.1 .AND. IPFLAG.NE.2 .AND. IPFLAG.NE.3 ) THEN
        PRINT *,' IPFLAG = ', IPFLAG
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
c     IF (PGNOTO('PGCONF')) RETURN
      IF (I1.LT.1 .OR. I2.GT.IDIM .OR. I1.GE.I2 .OR.
     :    J1.LT.1 .OR. J2.GT.JDIM .OR. J1.GE.J2) RETURN
      IF (C1.GE.C2) RETURN
      CALL PGBBUF
C
      RDX = (C1V2 - C1V1)/REAL( IDIM - 1 )
      RDY = (C2V2 - C2V1)/REAL( JDIM - 1 )
C
      IF ( IPFLAG.EQ.1 ) PHI = COORD
      IF ( IPFLAG.EQ.2 ) THE = COORD
      IF ( IPFLAG.EQ.3 ) RAD = COORD
C
      DO 140 J=J1+1,J2
         DO 130 I=I1+1,I2
            DVAL(1) = A(I-1,J)
            DVAL(2) = A(I-1,J-1)
            DVAL(3) = A(I,J-1)
            DVAL(4) = A(I,J)
            DVAL(5) = DVAL(1)
C
            IF ( DVAL( 1 ).GT.FLIMIT ) GOTO 130
            IF ( DVAL( 2 ).GT.FLIMIT ) GOTO 130
            IF ( DVAL( 3 ).GT.FLIMIT ) GOTO 130
            IF ( DVAL( 4 ).GT.FLIMIT ) GOTO 130
C
            XVAL(1) = C1V1 + REAL( I - 2 )*RDX
            XVAL(2) = C1V1 + REAL( I - 2 )*RDX
            XVAL(3) = C1V1 + REAL( I - 1 )*RDX
            XVAL(4) = C1V1 + REAL( I - 1 )*RDX
            XVAL(5) = XVAL(1)
C
            YVAL(1) = C2V1 + REAL( J - 1 )*RDY
            YVAL(2) = C2V1 + REAL( J - 2 )*RDY
            YVAL(3) = C2V1 + REAL( J - 2 )*RDY
            YVAL(4) = C2V1 + REAL( J - 1 )*RDY
            YVAL(5) = YVAL(1)
C
C Loop around the corners of our square to see
C if any corner is out of sight - in which
C case we move onto the next square
C
            DO IC = 1, 4
C             .
              IF ( IPFLAG.EQ.1 ) THEN
                RAD = XVAL( IC )
                THE = YVAL( IC )
              ENDIF
C             .
              IF ( IPFLAG.EQ.2 ) THEN
                RAD = XVAL( IC )
                PHI = YVAL( IC )
              ENDIF
C             .
              IF ( IPFLAG.EQ.3 ) THEN
                PHI = XVAL( IC )
                THE = YVAL( IC )
              ENDIF
C             .
              CALL SPHER_SAT_2_WORLD( RAD, THE, PHI, EARCM, XWORLD,
     1                              YWORLD, DEPTH )
C             .
              IF ( DEPTH.LT.0.0  .AND. IDEP.EQ.1 ) GOTO 130
C             .
            ENDDO
C
            NPT = 0
            DO 120 IC=1,4
               IF (DVAL(IC).GE.C1 .AND. DVAL(IC).LT.C2) THEN
                  NPT  = NPT+1
C
                  IF ( IPFLAG.EQ.1 ) THEN
                    RAD = XVAL( IC )
                    THE = YVAL( IC )
                  ENDIF
C                 .
                  IF ( IPFLAG.EQ.2 ) THEN
                    RAD = XVAL( IC )
                    PHI = YVAL( IC )
                  ENDIF
C                 .
                  IF ( IPFLAG.EQ.3 ) THEN
                    PHI = XVAL( IC )
                    THE = YVAL( IC )
                  ENDIF
C                 .
                  CALL SPHER_SAT_2_WORLD( RAD, THE, PHI, EARCM,
     1                              XWORLD, YWORLD, DEPTH )
C                 .
                  X(NPT) = XWORLD
                  Y(NPT) = YWORLD
C                 .
               END IF
               R = DVAL(IC+1)-DVAL(IC)
               IF (R.EQ.0.0) GOTO 120
               DO 110 LEV=1,2
                  IF (R.GT.0.0) THEN
                     C = C1
                     IF (LEV.EQ.2) C = C2
                  ELSE
                     C = C2
                     IF (LEV.EQ.2) C = C1
                  END IF
                  DELTA = (C-DVAL(IC))/R
                  IF (DELTA.GT.0.0 .AND. DELTA.LT.1.0) THEN
                     IF (IC.EQ.1 .OR. IC.EQ.3) THEN
                        XX = XVAL( IC )
                        YY = YVAL( IC ) +
     :                       DELTA*( YVAL( IC+1 ) - YVAL( IC ) )
                     ELSE
                        XX = XVAL( IC ) +
     :                       DELTA*( XVAL( IC+1 ) - XVAL( IC ) )
                        YY = YVAL( IC )
                     END IF
                     NPT = NPT+1
C                    .
                     IF ( IPFLAG.EQ.1 ) THEN
                       RAD = XX
                       THE = YY
                     ENDIF
C                    .
                     IF ( IPFLAG.EQ.2 ) THEN
                       RAD = XX
                       PHI = YY
                     ENDIF
C                    .
                     IF ( IPFLAG.EQ.3 ) THEN
                       PHI = XX
                       THE = YY
                     ENDIF
C                    .
                     CALL SPHER_SAT_2_WORLD( RAD, THE, PHI, EARCM,
     1                              XWORLD, YWORLD, DEPTH )
C                    .
                     X(NPT) = XWORLD
                     Y(NPT) = YWORLD
C
                  END IF
 110           CONTINUE
 120        CONTINUE
            IF (NPT.GE.3) CALL PGPOLY(NPT, X, Y)
 130     CONTINUE
 140  CONTINUE
      CALL PGEBUF
      END
C*********************************************************************
