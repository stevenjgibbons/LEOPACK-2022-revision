C*********************************************************************
C subroutine Spherical Harmonic Coefficient FiND *********************
C            -         -        -           - -- *********************
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
C                 FUNC( THETA, PHI, INFO, PARAMS )                   C
C                 where PARAMS and INFO are double precision and     C
C                 integer arrays of dimension ( * ) respectively.    C
C                 Both PARAMS and INFO are just for                  C
C                 generality of the function.                        C
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
C     SHC       : Spherical Harmonic Coefficients.                   C
C                  Has the usual ordering.                           C
C                   SHC( l*l ) = P_l^0                               C
C                   SHC( l*l + 2*m - 1 ) = P_l^mc ( m non zero )     C
C                   SHC( l*l + 2*m ) = P_l^ms ( m non zero )         C
C                  Dimension is ( LH*( LH + 2) )                     C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SF        : Scalar Function. This is an array containing a     C
C                  function over a set of theta and phi points.      C
C                    Dimensions are                                  C
C                      ( NTHPTS , NPHPTS )                           C
C     PA        : Schmidt Normalised Legendre Functions Dimension.   C
C                  {  ( LH + 1 )*( LH + 2 )/2 , NTHPTS }             C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C     FTF       : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Corresponding weights. As above.                   C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHCFND( SHC, SF, GAUX, GAUW, PA, DPA, FTF, FUNC, IDF,
     1                   INFO, PARAMS, LH, NTHPTS, NPHPTS, ZCOEF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS, NPHPTS, INFO( * ), IDF
      DOUBLE PRECISION SHC( LH*( LH + 2) ),
     1                 SF( NTHPTS , NPHPTS ),
     2                 GAUX( NTHPTS ), GAUW(NTHPTS),
     3                 PARAMS( * )
      DOUBLE PRECISION FTF( 2*NPHPTS ), ZCOEF,
     1                 PA ( ( LH + 1 )*( LH + 2 )/2, NTHPTS ),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2, NTHPTS )
      DOUBLE PRECISION FUNC
      EXTERNAL FUNC
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IPHI, ITHETA, IHARM, L, M, ICS
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
C ................ fill SF grid with data from FUNC ..................
      DO IPHI = 1, NPHPTS
         PHI = 2.0d0*PI*DBLE( IPHI - 1 )/NPHPTS
         DO ITHETA = 1, NTHPTS
            THETA = ACOS( GAUX(ITHETA) )
            SF( ITHETA, IPHI ) = FUNC( THETA, PHI, INFO, PARAMS )
         ENDDO
      ENDDO
C ............... transform grid data into Spherical Harmonic coeffs .
C
      CALL FORSSA ( SHC, SF, GAUW, PA, FTF, LH, LH,
     1              NTHPTS, NTHPTS, NPHPTS, NPHPTS, ZCOEF )
C
      IF ( IDF.EQ.0 ) RETURN
C
C     Output coefficients ...
C
      IF ( ABS( ZCOEF ).GT.TOL ) THEN
         TYPE = 'COS'
         WRITE ( IDF , 789 ) 0, 0, TYPE, ZCOEF
      ENDIF
      DO IHARM = 1, LH*(LH + 2)
         CALL LMFIND( IHARM , L, M, ICS )
         IF ( ICS.EQ.1 ) TYPE = 'COS'
         IF ( ICS.EQ.2 ) TYPE = 'SIN'
         IF ( ABS( SHC(IHARM)  ).GT.TOL ) THEN
           WRITE ( IDF , 789 ) L, M, TYPE, SHC(IHARM)
         ENDIF
 789     FORMAT('HARM_',I2,'^',I2,' ',A3,' = ',D14.5 )
      ENDDO
C
      RETURN
      END
C*********************************************************************
