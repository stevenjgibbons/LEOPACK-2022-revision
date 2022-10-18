C*********************************************************************
C subroutine FORward Scalar Spherical transform (Optimised) **********
C            ---     -      -                    -                   C
C Steve Gibbons Mon Oct  8 12:06:34 WEST 2001                        C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SF	: Scalar Function. This is an array containing a     C
C                  function over a set of theta and phi points.      C
C                    Dimensions are                                  C
C                      ( NPHPTS , NTHPTS )                           C
C                                                                    C
C                 SF( i, j ) contains the value a scalar function    C
C                 at theta = acos( x_j ) where x_j is the j^{th}     C
C                 element of GAUX returned by GAUWTS and             C
C                 phi = 2 pi ( i - 1 ) / nphpts                      C
C                                                                    C
C     GAUW      : Gauss weights as calculated by GAUWTS. ( NTHPTS )  C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     ZCOEF     : Coefficient of the l=0 spherical harmonic.         C
C                                                                    C
C  Integer                                                           C
C  -------							     C
C                                                                    C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum spherical harmonic order, m.               C
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SHC	: Spherical Harmonic Coefficients.                   C
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
C     FTF       : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FORSSO( SHC, SF, GAUW, PA, FTF, LH, MMAX,
     1                    NTHPTS, NPHPTS, ZCOEF, M0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          LH, NTHPTS, NPHPTS, MMAX, M0
      DOUBLE PRECISION SHC( LH*( LH + 2) ),
     1                 SF( NPHPTS, NTHPTS ), GAUW( NTHPTS )
      DOUBLE PRECISION FTF( 2*NPHPTS ), ZCOEF,
     1                 PA ( (LH+1)*(LH+2)/2, NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER          IOP, ITHETA, ILEV, INDCOS, INDSIN, LENV,
     1                 ISIGN, I, M, IP, MM, ISHCOS, ISHSIN
      DOUBLE PRECISION ZERO,WEIGHT
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C    
      ISIGN = 1
C ... zero the array of spherical harm coeffs --> IOP = 0
      IOP = 0
      LENV = LH*(LH+2)
      CALL VECOP ( SHC, ZERO, LENV, IOP )
      ZCOEF = 0.0d0
C ............... begin looping around theta points ..................
      DO ITHETA = 1, NTHPTS
         LENV = 2*NPHPTS
         CALL VECOP( FTF, ZERO, LENV, IOP )
         WEIGHT = GAUW( ITHETA )
C ................ now (considering equation we have w_j,
C all of the P_l^m(cos(theta_j) ) so now need the Fourier transform
C of the SF for (theta_j) ......
         DO I = 1, NPHPTS
            FTF ( 2*I - 1 ) = SF( I, ITHETA )
            ZCOEF = ZCOEF + WEIGHT * SF( I, ITHETA )/DBLE(NPHPTS)
         ENDDO
C ................ o.k. now let's Fourier Transform .................
C Array ftf currently contains value of function at phi point
C (i) in element ftf( 2i - 1 ) ...
C
         CALL FFTRLV( FTF, NPHPTS, ISIGN )
C
C ftf now contains c_m in ftf(2m+1) and s_m in ftf(2m+2)
C ................ let's loop around the harmonics ..................
         DO ILEV = 1, LH         
C ................................... first let's do the M = 0 case .
            IP     = ILEV*(ILEV+1)/2+1
            ISHCOS = ILEV*ILEV
            INDCOS = 1
            SHC ( ISHCOS ) = SHC( ISHCOS ) + 
     1                         PA( IP, ITHETA )*WEIGHT*FTF( INDCOS )
C ............ now let's loop around M from 1 to ILEV ..............
            DO MM = 1, MIN( ILEV, MMAX )/M0
               M  = MM*M0
               IP = ILEV*(ILEV+1)/2+M+1
               ISHSIN = ILEV*ILEV + 2*M
               ISHCOS = ISHSIN          - 1
               INDCOS = 1 + 2*MM
               INDSIN = 1 + INDCOS
               SHC(ISHCOS) = SHC(ISHCOS) + PA(IP, ITHETA)*
     1                        WEIGHT*FTF(INDCOS)
               SHC(ISHSIN) = SHC(ISHSIN) + PA(IP, ITHETA)*
     1                        WEIGHT*FTF(INDSIN)
            ENDDO
         ENDDO
C .................. end loop around the harmonics ..................
      ENDDO
C ................. end looping around theta points ..................
C Now just loop around the SHC's to multiply by the schmidt normalisation
C and scale factor of fft.
      ZCOEF = ZCOEF / 2.0d0
      DO ILEV = 1, LH
         WEIGHT = DBLE( 2*ILEV + 1 )/4.0d0
         DO M = 0 , 2*MIN( ILEV, MMAX )
            SHC( ILEV*ILEV + M ) = SHC( ILEV*ILEV + M )*WEIGHT
         ENDDO
      ENDDO
      RETURN
      END
C*********************************************************************
