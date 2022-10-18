C Derived from Numerical Recipes routine and is fully 
C described in that book. It is only called from NUMINT
C here and the only modifications are the arrays
C DPPAR and INTPAR to work with more general functions
C than the numerical recipes version will allow.
C Steve Gibbons - Wed Sep  8 10:52:05 BST 1999
C*********************************************************************
      SUBROUTINE TRAPZD( FUNC, A, B, S, N, INTPAR, DPPAR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N, INTPAR( * )
      DOUBLE PRECISION FUNC, A, B, S, DPPAR( * )
      EXTERNAL FUNC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER J, IT, TNM
      DOUBLE PRECISION SUM, DEL, X, FA, FB, FX
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( N.EQ.1 ) THEN
         FA = FUNC( A, INTPAR, DPPAR )
         FB = FUNC( B, INTPAR, DPPAR )
         S = 0.5d0*(B-A)*( FA + FB )
      ELSE
         IT = 2**(N-2)
         TNM = IT
         DEL = (B-A)/TNM
         X = A + 0.5d0*DEL
         SUM = 0.0d0
         DO J = 1, IT
            FX = FUNC( X, INTPAR, DPPAR )
            SUM = SUM + FX
            X = X + DEL
         ENDDO
         S = 0.5d0*( S + (B - A)*SUM/TNM )
      ENDIF
C
      RETURN
      END
C*********************************************************************
 
