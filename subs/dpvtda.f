C*********************************************************************
C subroutine Double Precision Vector Time Derivative Add *************
C            -      -         -      -    -          -   *************
C Steve Gibbons 22.3.99                                              C
C____________________________________________________________________C
C Performs the operation of the Time derivative matrix x a velocity  C
C vector ( without the need to actually form the matrix ) to form a  C
C right hand vector with time derivs.                                C
C                                                                    C
C The multiplication is of the form ( Y = alpha*A*X + beta*Y )       C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     RVEC      : DP vector of dimension ( NDIM )                    C
C     OLDVEC    : DP vector of dimension ( NDIM ) - old vector       C
C     FACT      : Multiplication factor for temperature term         C
C     FACV      : Multiplication factor for velocity term            C
C     ALPHA     : Multiplies the A matrix ( see above )              C
C     BETA      : Multiplies the Y vector ( see above )              C
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
      SUBROUTINE DPVTDA ( NDIM, NR, NH, MHT, MHL, RI, RO, RVEC,
     1                    OLDVEC, FACT, FACV, ORD, ALPHA, BETA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM, NR, NH,
     1        MHT( NH ), MHL( NH )
      DOUBLE PRECISION RI, RO, RVEC( NDIM ), FACT, FACV,
     1                 OLDVEC( NDIM ), ALPHA, BETA
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION H, RAD, DL, FAC, TOL,
     1                 D0F, D1F, D2F, D3F, D4F
      INTEGER IIH, INDI, L, IRN
      PARAMETER ( TOL = 1.0d-6 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the length of the vector
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine DPVTDA. Bad dimensions '
         PRINT *,' Array length = ',NDIM
         PRINT *,' NR = ', NR,' NH = ',NH
         STOP
      ENDIF
C
      FAC = 0.0d0
      H = ( RO - RI ) / DBLE( NR - 1 )
C
c     . Loop around the harmonics
c
      DO IIH = 1, NH
C First do poloidal velocity terms
        IF ( MHT( IIH ).EQ.1 ) THEN
c
           L = MHL( IIH )
           DO IRN = 1, NR
             INDI = ( IRN - 1 )*NH + IIH
             IF ( ABS( FACV ).GT.TOL ) THEN
               RAD = RI + H*DBLE( IRN - 1 )
               CALL SVDERF ( NDIM, NH, IIH, NR, IRN, OLDVEC, ORD,
     1                    D0F, D1F, D2F, D3F, D4F, RI, RO )
               FAC  = DL( L, RAD, D0F, D1F, D2F )*(-1.0d0)
               RVEC( INDI ) = BETA*RVEC( INDI ) + FAC*FACV*ALPHA
             ELSE
               RVEC( INDI ) = BETA*RVEC( INDI )
             ENDIF
           ENDDO
C
        ENDIF
C Now do toroidal velocity terms
        IF ( MHT( IIH ).EQ.2 ) THEN
c
           DO IRN = 1, NR
             INDI = ( IRN - 1 )*NH + IIH
             RVEC( INDI ) = BETA*RVEC( INDI ) +
     1                       FACV*ALPHA*OLDVEC( INDI )
           ENDDO
C
        ENDIF
C Now do temperature terms
        IF ( MHT( IIH ).EQ.3 ) THEN
c
           DO IRN = 1, NR
             INDI = ( IRN - 1 )*NH + IIH
             RVEC( INDI ) = BETA*RVEC( INDI ) +
     1                       FACT*ALPHA*OLDVEC( INDI )
           ENDDO
C
        ENDIF
      ENDDO
C
      RETURN
      END
C*********************************************************************
