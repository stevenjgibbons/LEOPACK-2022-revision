C*********************************************************************
C subroutine INVerse Scalar Spherical transform (Optimised) **********
C            ---     -      -                    -          **********
C Steve Gibbons Mon Oct  8 12:51:00 WEST 2001                        C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SHC       : Spherical Harmonic Coefficients.                   C
C                  Has the usual ordering.                           C
C                   SHC( l*l ) = P_l^0                               C
C                   SHC( l*l + 2*m - 1 ) = P_l^mc ( m non zero )     C
C                   SHC( l*l + 2*m ) = P_l^ms ( m non zero )         C
C                  Dimension is ( LH*( LH + 2) )                     C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C                  where X_j is the cosine of theta_j.               C
C                  PA is given by the routine SCHNLA and if SF       C
C                  is to be transformed back into spectral space,    C
C                  the GAUX array passed into SCHNLA must be output  C
C                  from GAUWTS.                                      C
C                                                                    C
C     ZCOEF  	: Scalar value of monopole term.                     C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     M0        : Smallest non-zero wavenumber.                      C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum spherical harmonic order, m.               C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SF        : Scalar Function. This is an array containing a     C
C                  function over a set of theta and phi points.      C
C                    Dimensions are                                  C
C                      ( NPHPTS , NTHPTS )                           C
C                 SF( i, j ) contains the value a scalar function    C
C                 at theta = acos( x_j ) where x_j is the j^{th}     C
C                 element of GAUX as supplied to SCHNLA and          C
C                 phi = 2 pi ( i - 1 ) / nphpts                      C
C                                                                    C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF       : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE INVSSO( SHC, SF, PA, FTF, ZCOEF, LH, MMAX,
     1                   NTHPTS, NPHPTS, M0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          LH, MMAX, NTHPTS, NPHPTS, M0
      DOUBLE PRECISION SHC( LH*( LH + 2) ), SF( NPHPTS, NTHPTS ),
     2                 FTF( 2*NPHPTS )
      DOUBLE PRECISION ZCOEF, PA( (LH+1)*(LH+2)/2, NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHETA, IPHI, L, M, INDPOL, IZERO, IND1, MM,
     1        INDFTC, INDFTS, ISIGN, LENV
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0 ,IZERO = 0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C   
      ISIGN = -1
C
C ......... first loop around theta points .........................
      DO ITHETA = 1, NTHPTS
C     ...................... Now nullify the transform vector ......
         LENV = 2*NPHPTS
         CALL VECOP ( FTF, ZERO, LENV, IZERO )
C        ......... Now loop around the harmonics for a given theta .
C        ... First do a summation over M ...........................
C        ............................................... case M = 0
         INDFTC = 1
         DO L = 1, LH
            INDPOL = L*(L+1)/2+1
            IND1 = L*L
            FTF ( INDFTC ) = FTF ( INDFTC ) +
     1                       PA( INDPOL, ITHETA )*SHC( IND1 )
         ENDDO
C        ............................................... case M > 0
         DO MM = 1, MMAX/M0
            M  = MM*M0
            INDFTC = 2*MM + 1
            INDFTS = 2*MM + 2
            DO L = M, LH
               INDPOL = L*(L+1)/2+M+1
               IND1 = L*L + 2*M - 1
               FTF ( INDFTC ) = FTF ( INDFTC ) +
     1                          PA( INDPOL, ITHETA )*SHC( IND1 )
               IND1 = IND1 + 1
               FTF ( INDFTS ) = FTF ( INDFTS ) +
     1                          PA( INDPOL, ITHETA )*SHC( IND1 )
            ENDDO
         ENDDO
C        ... ended summation over M ................................
C        ....... Ended loop around the harmonics for a given theta .
C .... now I think we are ready to do the Fourier transform (inverse)
         CALL FFTRLV( FTF, NPHPTS, ISIGN )
C .... now put the coefficients back into SF
         DO IPHI = 1, NPHPTS
            SF( IPHI, ITHETA ) = FTF( 2*IPHI - 1 )
            SF( IPHI, ITHETA ) = SF( IPHI, ITHETA ) + ZCOEF
         ENDDO
C
      ENDDO
C ......... ended loop around theta points .........................
C
      RETURN
      END
C*********************************************************************

