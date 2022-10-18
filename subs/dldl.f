C*********************************************************************
C function DLDL ******************************************************
C          ---- ******************************************************
C Steve Gibbons 10.4.97						     C
C____________________________________________________________________C
C Spherical Harmonic Differential Operator Dl^2			     C
C 								     C
C Dl f = f'' + 2/r f' - l(l+1)/(rr).f				     C
C  and so from eqn (74) we see					     C
C          d^4 f    4 d^3 f    2l(l+1) d^2 f    (l+2)(l+1)l(l-1) f   C
C Dl^2 f = -----  + - ----- -  ------- ----- +  ----------------     C
C           dr^4    r dr^3       r^2   dr^2          r^4             C
C____________________________________________________________________C
C Input variables :-						     C
C ===============						     C
C  Double Precision						     C
C  ----------------						     C
C     R		: Radius					     C
C     F		: Value of function at r(i)			     C
C     DDF	: Value of second derivative at r(i)		     C
C     DDDF	: Value of third derivative at r(i)		     C
C     DDDDF	: Value of fourth derivative at r(i)		     C
C  Integer						             C
C  -------						             C
C     L		: Spherical harmonic degree			     C
C____________________________________________________________________C
C Output :-						    	     C
C ======						  	     C
C  Double Precision						     C
C  ----------------						     C
C     DLDL	: as above					     C
C____________________________________________________________________C
C*********************************************************************
      FUNCTION DLDL(L,R,F,DDF,DDDF,DDDDF)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION DLDL,R,F,DDF,DDDF,DDDDF
      INTEGER L
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION LFACT0,LFACT2,LFACT3,TOL,Q
      PARAMETER (TOL=1.0d-10)
C All local variables confirmed as being explicity set
C by program.
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Avoid division by zero ....
      IF (ABS(R).LT.TOL) THEN
         WRITE (6,987)
         WRITE (6,988)
         WRITE (6,989)
 987     FORMAT ('Function DLDL has received too small a ')
 988     FORMAT ('value of R. This would result in a ')
 989     FORMAT ('division by zero error. Program aborted. ')
         STOP
      ENDIF
C
      Q = DBLE(L)
      LFACT0 = (Q+2.0d0)*(Q+1.0d0)*Q*(Q-1.0d0)/(R*R*R*R)
      LFACT2 = -2.0d0*Q*(Q+1.0d0)/(R*R)
      LFACT3 = 4.0d0/R
      DLDL = DDDDF+LFACT3*DDDF+LFACT2*DDF+LFACT0*F
      RETURN
      END
C*********************************************************************
