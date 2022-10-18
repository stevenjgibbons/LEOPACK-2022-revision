C*********************************************************************
C function DL ********************************************************
C          -- ********************************************************
C Steve Gibbons 8.4.97						     C
C____________________________________________________________________C
C Spherical Harmonic Differential Operator Dl			     C
C 								     C
C Dl f = f'' + 2/r f' - l(l+1)/(rr).f				     C
C____________________________________________________________________C
C Input variables :-						     C
C ===============						     C
C  Double Precision						     C
C  ----------------						     C
C     R		: Radius					     C
C     F		: Value of function at r(i)			     C
C     DF	: Value of first derivative at r(i)		     C
C     DDF	: Value of second derivative at r(i)		     C
C  Integer						             C
C  -------						             C
C     L		: Spherical harmonic degree			     C
C____________________________________________________________________C
C Output :-						    	     C
C ======						  	     C
C  Double Precision						     C
C  ----------------						     C
C     DL	: as above					     C
C____________________________________________________________________C
C*********************************************************************
      FUNCTION DL(L,R,F,DF,DDF)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION DL,R,F,DF,DDF
      INTEGER L
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION FACT,TOL
      PARAMETER (TOL=1.0d-10)
C Note that both of these working variables are confirmed to
C be explicitly set sjg 15.3.97 .
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Avoid division by zero ....
      IF (ABS(R).LT.TOL) THEN
         WRITE (6,987)
         WRITE (6,988)
         WRITE (6,989)
 987     FORMAT ('Function DL has received too small a ')
 988     FORMAT ('value of R. This would result in a ')
 989     FORMAT ('division by zero error. Program aborted. ')
         STOP
      ENDIF
C
      FACT=DBLE(L*L+L)
      DL=DDF+2.0d0*DF/R-FACT*F/(R*R)
      RETURN
      END
C*********************************************************************
