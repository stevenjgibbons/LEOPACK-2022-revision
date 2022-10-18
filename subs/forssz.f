C*********************************************************************
C subroutine FORward Scalar Spherical transform (aZimuthal sym.) *****
C            ---     -      -                     -                  C
C Steve Gibbons Sat Jul 22 12:02:15 BST 2000                         C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SFZ       : Scalar Function. This is an array containing a     C
C                  function over a set of theta and phi points.      C
C                    Dimensions are                                  C
C                      ( NPHPV, NTHPTS )                             C
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
C     NTHPTS    : Number of theta points.                            C
C     NPHPR     : Total number of phi points needed to perform       C
C                 Fourier transform.                                 C
C     NPHPV     : Number of phi points required to fully recreate    C
C                 the behaviour of the solution in a subspace of the C
C                 total physical space.                              C
C     M0        : Smallest non-zero wavenumber in solution.          C
C                 We check that NPHPV*M0 = NPHPR.                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum spherical harmonic order, m.               C
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
C                  Has dimensions ( 2*NPHPR )                        C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FORSSZ ( SHC, SFZ, GAUW, PA, FTF, LH, MMAX, M0,
     1                    NTHPTS, NPHPV, NPHPR, ZCOEF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          LH, NTHPTS, NPHPV, NPHPR, MMAX, M0
      DOUBLE PRECISION SHC( LH*( LH + 2) ),
     1                 SFZ( NPHPV, NTHPTS ), GAUW(NTHPTS)
      DOUBLE PRECISION FTF( 2*NPHPR ), ZCOEF,
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IOP,ITHETA,ILEV,INDCOS,INDSIN,LENV,ISIGN,I,M, IP
      DOUBLE PRECISION ZERO,WEIGHT
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C    
C Check the validity of arguments .....
C (Don't bother checking NPHPR is a power of 2
C as this will be done by FFTRLV)
C
      IF ( NTHPTS.LE.LH ) THEN
         PRINT *,' Subroutine FORSSZ.'
         PRINT *,' NTHPTS must be atleast ', LH + 1
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPR.LE.(2*MMAX+1) ) THEN
         PRINT *,' Subroutine FORSSZ.'
         PRINT *,' NPHPR must be atleast ',2*MMAX+1
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPV*M0.NE.NPHPR ) THEN
         PRINT *,' Subroutine FORSSZ.'
         PRINT *,' M0    = ', M0
         PRINT *,' NTHPV = ', NPHPV
         PRINT *,' NTHPR = ', NPHPR
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C_____________________________________________________________________
C
      ISIGN = 1
C ... zero the array of spherical harm coeffs --> IOP = 0
      IOP = 0
      LENV = LH*(LH+2)
      CALL VECOP ( SHC, ZERO, LENV, IOP )
      ZCOEF = 0.0d0
 
C ............... begin looping around theta points ..................
      DO ITHETA = 1, NTHPTS
         LENV = 2*NPHPR
         CALL VECOP ( FTF, ZERO, LENV, IOP )
         WEIGHT = GAUW( ITHETA )
C ................ now (considering equation we have w_j,
C all of the P_l^m(cos(theta_j) ) so now need the Fourier transform
C of the SF for (theta_j) ......
         DO I = 1, NPHPV
            FTF ( 2*I - 1 ) = SFZ( I, ITHETA )
            ZCOEF = ZCOEF + WEIGHT*SFZ( I, ITHETA )/DBLE(NPHPV)
         ENDDO
         CALL FTAASC( NPHPR, NPHPV, FTF )
C ................ o.k. now let's Fourier Transform .................
C Array ftf currently contains value of function at phi point
C (i) in element ftf( 2i - 1 ) ...
C
         CALL FFTRLV ( FTF, NPHPR, ISIGN )
C
C ftf now contains c_m in ftf(2m+1) and s_m in ftf(2m+2)
C ................ let's loop around the harmonics ..................
         DO ILEV = 1, LH         
            INDCOS = ILEV*ILEV
C ................................... first let's do the M = 0 case .
            M = 0
            IP = ILEV*(ILEV+1)/2+M+1
            SHC ( INDCOS ) = SHC( INDCOS ) + 
     1                       PA( IP, ITHETA )*WEIGHT*FTF( 1 )
            INDCOS = INDCOS + 1
            INDSIN = INDCOS + 1
C ............ now let's loop around M from 1 to ILEV ..............
            DO M = 1, MIN( ILEV, MMAX )
               IF ( MOD( M, M0 ).NE.0 ) GOTO 67
               IP = ILEV*(ILEV+1)/2+M+1
               SHC(INDCOS) = SHC(INDCOS) + PA(IP, ITHETA)*
     1                        WEIGHT*FTF(2*M+1)
               SHC(INDSIN) = SHC(INDSIN) + PA(IP, ITHETA)*
     1                        WEIGHT*FTF(2*M+2)
 67            CONTINUE
               INDCOS = INDCOS + 2
               INDSIN = INDSIN + 2
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
            SHC( ILEV*ILEV + M ) = SHC ( ILEV*ILEV + M ) * WEIGHT
         ENDDO
      ENDDO
      RETURN
      END
C*********************************************************************
