C*********************************************************************
C d. p. function Dynamo Benchmark Project Initial State Components ***
C                -      -         -       -       -     -          ***
C Steve Gibbons Thu Feb  3 17:52:23 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Supplies values of components given for initial state for          C
C dynamo benchmark as suggested by U. Christensen.                   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ICOMP     : Selects the component to be calculated.            C
C                                                                    C
C                 (1) --> B_r (radial magnetic field)                C
C                 (2) --> B_theta (theta magnetic field)             C
C                 (3) --> B_phi (phi magnetic field)                 C
C                 (4) --> Theta scalar function.                     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Distance from centre of sphere.                    C
C     THE       : Co-latitude in radians.                            C
C     PHI       : Longitude in radians.                              C
C     RI        : Radius of inner boundary                           C
C     RO        : Radius of outer boundary                           C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DBPISC( ICOMP, RAD, THE, PHI, RI, RO )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION DBPISC, RAD, THE, PHI, RI, RO
      INTEGER ICOMP
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION PI, LOW, COSTH, RAD3, RI4, PIRMRI, TWOTHE,
     1                 X, X2, X4, X6, SINTH, SINTH4, FOURPH, BRAC,
     2                 A, FRAC
      PARAMETER ( PI = 3.14159265358979312D0, LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( RAD.LT.LOW ) THEN
        PRINT *,' Function DBPISC.'
        PRINT *,' RAD = ', RAD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ICOMP.NE.1 .AND. ICOMP.NE.2 .AND.
     1     ICOMP.NE.3 .AND. ICOMP.NE.4 ) THEN
        PRINT *,' Function DBPISC.'
        PRINT *,' ICOMP = ', ICOMP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Case of B_r
C
      IF ( ICOMP.EQ.1 ) THEN
        RI4    = RI*RI*RI*RI
        RAD3   = RAD*RAD*RAD
        COSTH  = COS( THE )
        BRAC   = 6.0d0*RAD - 8.0d0*RO + 2.0d0*RI4/RAD3
        DBPISC = -5.0d0*BRAC*COSTH/8.0d0
        RETURN
      ENDIF
C
C Case of B_theta
C
      IF ( ICOMP.EQ.2 ) THEN
        RI4    = RI*RI*RI*RI
        RAD3   = RAD*RAD*RAD
        SINTH  = SIN( THE )
        BRAC   = -9.0d0*RAD + 8.0d0*RO + RI4/RAD3
        DBPISC = -5.0d0*BRAC*SINTH/8.0d0
        RETURN
      ENDIF
C
C Case of B_phi
C
      IF ( ICOMP.EQ.3 ) THEN
        PIRMRI = PI*( RAD - RI )
        TWOTHE = 2.0d0*THE
        DBPISC = 5.0d0*SIN( PIRMRI )*SIN( TWOTHE )
        RETURN
      ENDIF
C
C Case of Theta
C
      IF ( ICOMP.EQ.4 ) THEN
        X      = 2.0d0*RAD - RI - RO
        X2     = X*X
        X4     = X2*X2
        X6     = X4*X2
        BRAC   = ( 1.0d0 - 3.0d0*X2 + 3.0d0*X4 - X6 )
        SINTH  = SIN( THE )
        SINTH4 = SINTH*SINTH*SINTH*SINTH
        FOURPH = 4.0d0*PHI
        A      = 0.1d0
        FRAC   = 210.0d0/SQRT( 17920.0d0*PI )
        DBPISC = A*FRAC*BRAC*SINTH4*COS( FOURPH )
        RETURN
      ENDIF
C
      END
C*********************************************************************
