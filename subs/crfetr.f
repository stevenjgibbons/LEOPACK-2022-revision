C*********************************************************************
C subroutine Chebyshev Radial Function Energy Trapezoidal Rule *******
C            -         -      -        -      -           -    *******
C Steve Gibbons Sat May  6 09:19:53 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C A direct adaption of the numerical recipes function TRAPZD which   C
C is specifically designed only to integrate the function CRFEIF.    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CRFETR( A, B, S, N, INARR, IPARR, IH, CCV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N, INARR( * ), IPARR, IH
      DOUBLE PRECISION A, B, S, CCV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER J, IT, TNM
      DOUBLE PRECISION SUM, DEL, X, FA, FB, FX, CRFEIF
      EXTERNAL CRFEIF
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( N.EQ.1 ) THEN
         FA = CRFEIF( INARR, IPARR, IH, A, CCV )
         FB = CRFEIF( INARR, IPARR, IH, B, CCV )
         S = 0.5d0*(B-A)*( FA + FB )
      ELSE
         IT = 2**(N-2)
         TNM = IT
         DEL = (B-A)/TNM
         X = A + 0.5d0*DEL
         SUM = 0.0d0
         DO J = 1, IT
            FX = CRFEIF( INARR, IPARR, IH, X, CCV )
            SUM = SUM + FX
            X = X + DEL
         ENDDO
         S = 0.5d0*( S + (B - A)*SUM/TNM )
      ENDIF
C
      RETURN
      END
C*********************************************************************
 
