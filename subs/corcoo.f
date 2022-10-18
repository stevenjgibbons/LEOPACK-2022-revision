C*********************************************************************
C subroutine CORiolis force COefficient find (Optimised) *************
C            ---            --                -          *************
C Steve Gibbons Tue Oct  9 10:28:38 WEST 2001                        C
C____________________________________________________________________C
C Takes an input of a single Q, S or T harmonic and returns the      C
C array QST with the appropriate coefficients for the coriolis force C
C of the single starting vector harmonic.                            C
C____________________________________________________________________C
C                                                                    C
C PRE-REQUISITES                                                     C
C --------------                                                     C
C     The Gaussian Weights and Legendre Polynomials must be found    C
C using routines GAUWTS and SCHNLA. Use the calling sequence ...     C
C                                                                    C
C >      X1 = -1.0d0                                                 C
C >      X2 =  1.0d0                                                 C
C >                                                                  C
C >      CALL GAUWTS ( X1, X2, GAUX, GAUW, NTHPTS, NTHMAX )          C
C >                                                                  C
C >      CALL SCHNLA ( PA, DPA, GAUX, LH, LHMAX, NTHPTS, NTHMAX)     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IQST      : IQST = 1 for starting with a Q harmonic.           C
C                 IQST = 2 for starting with an S harmonic.          C
C                 IQST = 3 for starting with a T harmonic.           C
C     L    	: Spherical harmonic degree, l.                      C
C     M    	: Spherical harmonic order, m.                       C
C     ICS  	: 1 for cos m phi, 2 for sin m phi ...               C
C     LH   	: Maximum spherical harmonic degree, l.              C
C     NPHPTS	: The number of phi points.                          C
C     NTHPTS	: The number of theta points.                        C
C     MMAX 	: Maximum spherical harmonic order, m.               C
C     M0        : Minimum non-zero wavenumber.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     QST       : Output array containing scaloidal/spheroidal       C
C                  toroidal decomposition of vector.                 C
C                  Has dimensions ( LH*(LH+2), 3).                   C
C               QST (l*l+2m,1) = q_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,1) = q_l^mc(r_i)                         C
C               QST (l*l+2m,2) = s_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,2) = s_l^mc(r_i)                         C
C               QST (l*l+2m,3) = t_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,3) = t_l^mc(r_i)                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Working Arrays   :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF1	: Array for Fourier transforming dim. (2*NPHPTS)     C
C     FTF2	: Array for Fourier transforming dim. (2*NPHPTS)     C
C     FTF3	: Array for Fourier transforming dim. (2*NPHPTS)     C
C     VF	: Vector function. dimensions ( NPHPTS, NTHPTS, 3)   C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CORCOO( IQST, L, M, ICS, QST, LH, NPHPTS, NTHPTS,
     1        FTF1, FTF2, FTF3, VF, GAUX, GAUW, PA, DPA, MMAX, M0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          IQST, L, M, ICS, LH, NPHPTS, NTHPTS, MMAX, M0
      DOUBLE PRECISION 
     1                 QST( LH*(LH+2), 3),
     2                 VF( NPHPTS, NTHPTS, 3), FTF1( 2*NPHPTS ),
     3                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS )
      DOUBLE PRECISION
     1                 GAUX( NTHPTS ), GAUW( NTHPTS ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION ZCOEF
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Let's evaluate this harmonic in space, stored in array VF
C
      CALL SHVECO( L, M, ICS, IQST, VF, GAUX, PA, DPA,
     1             NTHPTS, NPHPTS, LH, M0 )
C
C Now let's take the coriolis vector on it ..............
      CALL VFCOR ( NTHPTS, NPHPTS, VF, GAUX )
C
C Let's convert back to QST values ........................
C
      CALL VF2QSO( QST, VF, GAUX, GAUW, PA, DPA, FTF1, FTF2, FTF3,
     1             ZCOEF, LH, NTHPTS, NPHPTS, MMAX, M0 )
C 
      RETURN
      END
C*********************************************************************
