C*********************************************************************
C subroutine QST coefficient FiND ************************************
C            ---             - -- ************************************
C Steve Gibbons 30.4.97						     C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     PARAMS    : Parameter list for FUNC. Dimension ( * )           C
C                                                                    C
C ( External Function )						     C
C .....................						     C
C                                                                    C
C     FUNC	: Function of theta and phi. Must have the form      C
C                 FUNC( THETA, PHI, INFO, PARAMS, ICOMP )            C
C                 where PARAMS and INFO are double precision and     C
C                 integer arrays of dimension ( * ) respectively.    C
C                 Both PARAMS and INFO are just for                  C
C                 generality of the function.                        C
C                                                                    C
C                 ICOMP = 1 for radial component                     C
C                 ICOMP = 2 for theta  component                     C
C                 ICOMP = 3 for phi    component                     C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C (both NTHPTS and NPHPTS must be powers of 2. This is checked for). C
C     LH        : Level of harmonics.                                C
C     INFO	: Information flag for FUNC                          C
C     IDF       : Integer display flag. Set to 0 if no display is    C
C                 Required. Otherwise set to an integer value which  C
C                 is the logical unit number of output file.         C
C                 i.e. IDF = 6 for standard output.                  C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     QST       : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector c.f. eqn (38).            C
C                  Has dimensions (  LH*(LH+2) , 3).                 C
C               QST (l*l+2m,1) = q_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,1) = q_l^mc(r_i)                         C
C               QST (l*l+2m,2) = s_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,2) = s_l^mc(r_i)                         C
C               QST (l*l+2m,3) = t_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,3) = t_l^mc(r_i)                         C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VF        : Vector Function. An array of dimensions            C
C                  ( NPHPTS, NTHPTS, 3) which contain the R, THETA   C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IPHI, ITHETA, 2 ) is the Theta      C
C                  compontent of the vector at (iphi, itheta).       C
C     PA        : Schmidt Normalised Legendre Functions Dimension.   C
C                  {  ( LH + 1 )*( LH + 2 )/2 , NTHPTS }             C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Corresponding weights. As above.                   C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE QSTFND ( QST, VF, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1                 FTF3, FUNC, INFO, PARAMS, LH, NTHPTS, NPHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS, NPHPTS, INFO( * ), IDF
      DOUBLE PRECISION QST( LH*( LH + 2), 3 ),
     1                 VF( NPHPTS , NTHPTS, 3 ),
     2                 GAUX( NTHPTS ), GAUW(NTHPTS),
     3                 PARAMS( * )
      DOUBLE PRECISION FTF1( 2*NPHPTS ), FTF2( 2*NPHPTS ),
     1                 FTF3( 2*NPHPTS ), 
     2                 PA ( ( LH + 1 )*( LH + 2 )/2, NTHPTS ),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2, NTHPTS )
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IPHI, ITHETA, ICOMP
      DOUBLE PRECISION PHI, THETA, TOL, X1, X2
      PARAMETER (TOL = 1.0d-10)
      CHARACTER *(3) TYPE
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Calculate Gauss points and weights
C
      X1 = -1.0d0
      X2 = 1.0d0
      CALL GAUWTS ( X1, X2, GAUX, GAUW, NTHPTS, NTHPTS )
      CALL SCHNLA ( PA, DPA, GAUX, LH, LH, NTHPTS, NTHPTS)
C
C ................ fill VF grid with data from FUNC ..................
      DO IPHI = 1, NPHPTS
         PHI = 2.0d0*PI*DBLE( IPHI - 1 )/NPHPTS
         DO ITHETA = 1, NTHPTS
            THETA = ACOS( GAUX(ITHETA) )
            DO ICOMP = 1, 3
              VF( IPHI, ITHETA, ICOMP ) =
     1          FUNC( THETA, PHI, INFO, PARAMS, ICOMP )
            ENDDO
         ENDDO
      ENDDO
C ............... transform grid data into Spherical Harmonic coeffs .
C
      CALL VF2QSA ( QST, VF, GAUX, GAUW, PA, DPA, FTF1,
     1              FTF2, FTF3, LH, LH, NTHPTS, NTHPTS,
     2              NPHPTS, NPHPTS )
C
      RETURN
      END
C*********************************************************************
