C*********************************************************************
C subroutine Radial QST Coriolis Force *******************************
C            -      --- -        -     *******************************
C Steve Gibbons Sat Oct  2 17:51:02 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Given a vector, v, in RQST format, RQSTCF will return (k \times v) C
C in RQST format. v is supplied in the array RQST1, and (k \times v) C
C is output RQST2.                                                   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NR        : Number of radial grid nodes.                       C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum mode number to be used.                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C     NTHPTS    : Number of points in theta.                         C
C     NPHPTS    : Number of points in phi.                           C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Gauss weights computed bu GAUWTS. Dim ( NTHPTS )   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C       PA and DPA must be formed in advance by a call to SCHNLA.    C
C                                                                    C
C  RQST1 and ZCF1 store the qst decomposition of input vector.       C
C  RQST2 and ZCF2 store the qst decomposition of the output vector.  C
C                                                                    C
C   ZCF. ( j ) contains the coeff. to Q_0^0( r_j )                   C
C                                                                    C
C     RQST1     : Dim. (  LH*(LH+2) ,3, NR )                         C
C     ZCF1      : Dim. (  NR )                                       C
C     RQST2     : Dim. (  LH*(LH+2) ,3, NR )                         C
C     ZCF2      : Dim. (  NR )                                       C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VF1       : Storage for vector function. Dim.(NPHPTS,NTHPTS,3) C
C     FTF1      : Array for fourier transforming. Dim. ( 2*NPHPTS )  C
C     FTF2      : Array for fourier transforming. Dim. ( 2*NPHPTS )  C
C     FTF3      : Array for fourier transforming. Dim. ( 2*NPHPTS )  C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RQSTCF( NR, LH, MMAX, ILNR, IRNR, NTHPTS, NPHPTS,
     1                   GAUX, GAUW, PA, DPA, RQST1, ZCF1, RQST2,
     2                   ZCF2, VF1, FTF1, FTF2, FTF3 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, LH, MMAX, ILNR, IRNR, NTHPTS, NPHPTS
      DOUBLE PRECISION RQST1( LH*(LH+2) ,3 ,NR ), ZCF1( NR ),
     1                 RQST2( LH*(LH+2) ,3 ,NR ), ZCF2( NR )
      DOUBLE PRECISION FTF1( 2*NPHPTS ),
     1                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
      DOUBLE PRECISION VF1( NPHPTS, NTHPTS, 3), GAUX( NTHPTS ),
     1                 GAUW( NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Loop around points IR = ILNR, IRNR
C
      DO IR = ILNR, IRNR
C       .
C       . Put the vector function of RQST1 into VF1 for node IR
C       .
        CALL RQSTVF( RQST1, VF1, GAUX, PA, DPA, FTF1, FTF2, FTF3,
     1               LH, NTHPTS, NPHPTS, MMAX, ZCF1, NR, IR )
C       .
C       . Replace VF1 with the appropriate Coriolis term.
C       .
        CALL VFCOR( NTHPTS, NPHPTS, VF1, GAUX )
C       .
C       . Transform back into RQST3
C       .
        CALL VFRQST( RQST2, VF1, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1               FTF3, ZCF2, LH, NTHPTS, NPHPTS, MMAX, NR, IR )
C       .
C
      ENDDO
C
      RETURN
      END
C*********************************************************************
