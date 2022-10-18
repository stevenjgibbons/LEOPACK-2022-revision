C*********************************************************************
C subroutine Double Precision Vector Stratification Term Add *********
C            -      -         -      -              -    -   *********
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C Adds the term to the double precision vector resulting from        C
C the V . grad ( t_0 ( r ) ) term in the heat equation.              C
C                                                                    C
C The general form is t_0 ( r ) = CB1 * r + CB2 / r^2 and so in      C
C effect, we are adding  ( CB1 + CB2 / r^3 ) to the poloidal         C
C velocity term.                                                     C
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
C   MHT( I ) = 4 for a poloidal magnetic field vector                C
C   MHT( I ) = 5 for a poloidal magnetic field vector                C
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
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     RVEC      : DP vector of dimension ( NDIM )                    C
C     OLDVEC    : DP vector of dimension ( NDIM ) - old vector       C
C     FAC       : Multiplication factor for term                     C
C                                                                    C
C     CB1       : See above                                          C
C     CB2       : See above                                          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPVSTA ( NDIM, NR, NH, MHT, MHL, MHM, MHC, RI, RO,
     1                    RVEC, OLDVEC, FAC, CB1, CB2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM, NR, NH,
     1        MHT( NH ), MHM( NH ), MHC( NH ), MHL( NH )
      DOUBLE PRECISION RI, RO, RVEC( NDIM ), FAC, CB1, CB2,
     1                 OLDVEC( NDIM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION H, POL, STFAC, RAD, TOL
      PARAMETER ( TOL = 1.0d-8 )
      INTEGER LI, MI, ICSI, IIH, IOH, INDI, INDO, IRN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the length of the vector
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine DPVSTA. Bad dimensions '
         PRINT *,' Array length = ',NDIM
         PRINT *,' NR = ', NR,' NH = ',NH
         STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FAC ).LT.TOL ) RETURN
C
      H = ( RO - RI ) / DBLE( NR - 1 )
c
c     . Loop around the velocity harmonics
c
      DO IIH = 1, NH
        IF ( MHT( IIH ).EQ.1 ) THEN
         LI   = MHL( IIH )
         MI   = MHM( IIH )
         ICSI = MHC( IIH )
         DO IOH = 1, NH
           IF (    MHT( IOH ).EQ.3       .AND.
     1             MHL( IOH ).EQ.LI      .AND.
     2             MHM( IOH ).EQ.MI      .AND.
     3             MHC( IOH ).EQ.ICSI    ) THEN
c
c          . we have now found the appropriate
c          . temperature harmonic. so add strat term.
c          . Loop around the radial grid nodes ...
c
            DO IRN = 1, NR
              RAD = RI + H*DBLE( IRN - 1 )
              INDI = ( IRN - 1 )*NH + IIH
              POL  = OLDVEC( INDI )
              INDO = ( IRN - 1 )*NH + IOH
              STFAC =  DBLE(LI*LI + LI)*( CB1 + CB2/RAD**3 )
              RVEC( INDO ) = RVEC( INDO ) + POL*STFAC*FAC
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
