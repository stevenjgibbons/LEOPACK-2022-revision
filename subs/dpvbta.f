C*********************************************************************
C subroutine Double Precision Vector Buoyancy Term Add ***************
C            -      -         -      -        -    -   ***************
C Steve Gibbons 20.6.99                                              C
C____________________________________________________________________C
C                                                                    C
C Add the buoyancy term to the right hand vector.                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NDIM      : Length of the solution vector.                     C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of harmonics (all types)                    C
C                                                                    C
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C     MHM       : Dimension ( NH )                                   C
C     MHC       : Dimension ( NH )                                   C
C                                                                    C
C  mht, mhl, mhm and mhc define the spherical harmonics present      C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RVEC      : DP vector of dimension ( NDIM )                    C
C     OVEC      : DP vector of dimension ( NDIM ) - old vector       C
C     FAC       : Multiplication factor for term                     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPVBTA ( NDIM, NR, NH, MHT, MHL, MHM, MHC,
     1                    RVEC, OVEC, FAC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM, NR, NH,
     1        MHT( NH ), MHM( NH ), MHC( NH ), MHL( NH )
      DOUBLE PRECISION RVEC( NDIM ), FAC, OVEC( NDIM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION POL, TOL
      PARAMETER ( TOL = 1.0d-8 )
      INTEGER LI, MI, ICSI, IIH, IOH, INDI, INDO, IRN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the length of the vector
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine DPVBTA. Bad dimensions '
         PRINT *,' Array length = ',NDIM
         PRINT *,' NR = ', NR,' NH = ',NH
         STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FAC ).LT.TOL ) RETURN
c
c     . Loop around the temperature harmonics
c
      DO IIH = 1, NH
         IF ( MHT( IIH ).EQ.3 ) THEN
         LI   = MHL( IIH )
         MI   = MHM( IIH )
         ICSI = MHC( IIH )
         DO IOH = 1, NH
           IF (    MHT( IOH ).EQ.1       .AND.
     1             MHL( IOH ).EQ.LI      .AND.
     2             MHM( IOH ).EQ.MI      .AND.
     3             MHC( IOH ).EQ.ICSI    ) THEN
c
c          . we have now found the appropriate
c          . poloidal harmonic. so add strat term.
c          . Loop around the radial grid nodes ...
c
            DO IRN = 1, NR
              INDI = ( IRN - 1 )*NH + IIH
              POL  = OVEC( INDI )
              INDO = ( IRN - 1 )*NH + IOH
              RVEC( INDO ) = RVEC( INDO ) + POL*FAC
            ENDDO
            GOTO 500
c
           ENDIF
         ENDDO
C
        ENDIF
 500    CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************
