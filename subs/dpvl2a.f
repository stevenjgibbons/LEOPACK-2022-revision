C*********************************************************************
C subroutine Double Precision Vector Laplacian 2 (toroidal) Add ******
C            -      -         -      -         -            -   ******
C Steve Gibbons 17.6.99                                              C
C____________________________________________________________________C
C Adds the term to the double precision vector resulting from        C
C the curl of the Laplacian of the toroidal velocity.                C
C                                                                    C
C This results in a poloidal harmonic with a D_l operator.           C
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
C                                                                    C
C  mht and mhl define the spherical harmonics present                C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     RVEC      : DP vector of dimension ( NDIM )                    C
C     OLDVEC    : DP vector of dimension ( NDIM ) - old vector       C
C     FAC       : Multiplication factor for term                     C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     ORD       : *(2). Determines the order of accuracy. Options :  C
C             SS - Strictly second order                             C
C             SF - Strictly fourth order                             C
C             O5 - Optimum accuracy for bandwidth 5; this gives      C
C                  Fourth order accuracy for 1st and 2nd derivatives C
C                  and second order accuracy for 3rd and 4th der.s   C
C             O7 - Optimum accuracy for bandwidth 7; this gives      C
C                  Sixth order accuracy for 1st and 2nd derivatives  C
C                  and fourth order accuracy for 3rd and 4th der.s   C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPVL2A ( NDIM, NR, NH, MHT, MHL, RI, RO,
     1                    RVEC, OLDVEC, FAC, ORD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM, NR, NH,
     1        MHT( NH ), MHL( NH )
      DOUBLE PRECISION RI, RO, RVEC( NDIM ), FAC,
     1                 OLDVEC( NDIM )
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION H, LAPFAC, RAD, D0F, D1F, D2F, D3F, D4F, DL,
     1                 TOL
      PARAMETER ( TOL = 1.0d-8 )
      INTEGER LI, IIH, INDO, IRN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the length of the vector
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine DPVL2A. Bad dimensions '
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
         IF ( MHT( IIH ).EQ.2 ) THEN
         LI   = MHL( IIH )
         DO IRN = 1, NR
           RAD = RI + H*DBLE( IRN - 1 )
c
           CALL SVDERF ( NDIM, NH, IIH, NR, IRN, OLDVEC, ORD,
     1                    D0F, D1F, D2F, D3F, D4F, RI, RO )
c
           INDO = ( IRN - 1 )*NH + IIH
           LAPFAC = DL(LI,RAD,D0F,D1F,D2F)
           RVEC( INDO ) = RVEC( INDO ) + FAC*LAPFAC
         ENDDO
        ENDIF
      ENDDO
C
      RETURN
      END
C*********************************************************************
