C*********************************************************************
C subroutine PGCONF_POLAR: Adapted from PGCONF ***********************
C            ------------                      ***********************
C Steve Gibbons Fri Jan 19 10:02:14 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C This subroutine is a very simple adaption from the PGPLOT routine. C
C The only difference is that the contour plot is projected onto a   C
C planar polar coordinate grid, instead of rectangular.              C
C                                                                    C
C However, it requires a common block PARAMP with two INTEGERs       C
C NRAD (number of equally spaced radial grid nodes) and              C
C NTHE (number of equally spaced theta grid nodes), and four REAL    C
C variables RFIRST, RLAST, TFIRST, TLAST                             C
C                                                                    C
C RFIRST is the inner boundary radius                                C
C RLAST is the outer boundary radius                                 C
C                                                                    C
C TFIRST is the first theta point, in radians.                       C
C TLAST is the last theta point, in radians.                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C See argument list below.                                           C
C                                                                    C
C*********************************************************************
      SUBROUTINE PGCONF_POLAR (A, IDIM, JDIM, I1, I2, J1, J2, C1, C2)
      IMPLICIT NONE
      INTEGER IDIM, JDIM, I1, I2, J1, J2
      REAL    A(IDIM,JDIM), C1, C2
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
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
      REAL     DVAL(5), X(8), Y(8), DELTA, XX, YY, C, R
      INTEGER  IDELT(6)
      DATA     IDELT/0,-1,-1,0,0,-1/
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      REAL R1, T1, XWORLD, YWORLD
C____________________________________________________________________C
C
C Check arguments.
C
c     IF (PGNOTO('PGCONF')) RETURN
      IF (I1.LT.1 .OR. I2.GT.IDIM .OR. I1.GE.I2 .OR.
     :    J1.LT.1 .OR. J2.GT.JDIM .OR. J1.GE.J2) RETURN
      IF (C1.GE.C2) RETURN
      CALL PGBBUF
C
      DO 140 J=J1+1,J2
         DO 130 I=I1+1,I2
            DVAL(1) = A(I-1,J)
            DVAL(2) = A(I-1,J-1)
            DVAL(3) = A(I,J-1)
            DVAL(4) = A(I,J)
            DVAL(5) = DVAL(1)
C
            NPT = 0
            DO 120 IC=1,4
               IF (DVAL(IC).GE.C1 .AND. DVAL(IC).LT.C2) THEN
                  NPT = NPT+1
                  XX = I+IDELT(IC+1)
                  YY = J+IDELT(IC)
C
                  R1 = RFIRST +
     1     (RLAST-RFIRST)*(XX - 1.0)/REAL( NRAD - 1 )
                  T1 = TFIRST +
     1     (TLAST-TFIRST)*(YY - 1.0)/REAL( NTHE - 1 )
C
                  XWORLD = R1*COS( T1 )
                  YWORLD = R1*SIN( T1 )
                  X(NPT) = XWORLD
                  Y(NPT) = YWORLD
C
c                 X(NPT) = TR(1) + TR(2)*XX + TR(3)*YY
c                 Y(NPT) = TR(4) + TR(5)*XX + TR(6)*YY
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
                        XX = I+IDELT(IC+1)
                        YY = REAL(J+IDELT(IC)) +
     :                       DELTA*REAL(IDELT(IC+1)-IDELT(IC))
                     ELSE
                        XX = REAL(I+IDELT(IC+1)) +
     :                       DELTA*REAL(IDELT(IC+2)-IDELT(IC+1))
                        YY = J+IDELT(IC)
                     END IF
                     NPT = NPT+1
C
                     R1 = RFIRST +
     1     (RLAST-RFIRST)*(XX - 1.0)/REAL( NRAD - 1 )
                     T1 = TFIRST +
     1     (TLAST-TFIRST)*(YY - 1.0)/REAL( NTHE - 1 )
C
                     XWORLD = R1*COS( T1 )
                     YWORLD = R1*SIN( T1 )
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
