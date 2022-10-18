C*********************************************************************
C subroutine GRadient INVerse transform: (aZimuthal symmetry) ********
C            --       ---                  -                  ********
C Steve Gibbons Sat Jul 22 13:43:45 BST 2000                         C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SHC       : Spherical Harmonic Coefficients.                   C
C                  Has the usual ordering.                           C
C                   SHC( l*l ) = P_l^0                               C
C                   SHC( l*l + 2*m - 1 ) = P_l^mc ( m non zero )     C
C                   SHC( l*l + 2*m ) = P_l^ms ( m non zero )         C
C                  Dimension is ( LH*( LH + 2) )                     C
C     DSHC      : Radial derivatives of the above.                   C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     RAD       : Radial value at which function is to be evaluated. C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2,NTHPTS)       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C     DZCOEF    : Derivative of l=0 harmonic radial function.        C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
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
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GRFNZ     : GRadient FuNction                                  C
C                 has dimensions ( 3, NPHPV, NTHPTS )                C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Have dimensions ( 2*NPHPR )                       C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE GRINVZ( SHC, DSHC, GAUX, RAD, PA, DPA, FTF1, FTF2,
     1                   FTF3, GRFNZ, LH, NTHPTS, NPHPV, NPHPR, M0,
     2                   MMAX, DZCOEF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS, NPHPV, NPHPR, M0, MMAX
      DOUBLE PRECISION RAD, SHC( LH*( LH + 2) ),
     1                 DSHC( LH*( LH + 2) ),
     2                 GRFNZ( 3, NPHPV, NTHPTS ),
     3                 GAUX( NTHPTS ), DZCOEF
      DOUBLE PRECISION FTF1( 2*NPHPR ),
     1                 FTF2( 2*NPHPR ), FTF3( 2*NPHPR ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2, NTHPTS ),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2, NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, ISIGN, ITHREE, ITHETA, LENV, IPHI,
     1        IP, L, M, IND1, IND2, IND
      DOUBLE PRECISION ZERO,X,SINE,THETA,PSI, DPSI, PSI2, DPSI2, TOL
      PARAMETER ( ZERO = 0.0d0, ITHREE=3, TOL = 1.0d-6 )
C____________________________________________________________________C
C Functions used :-
c     INTEGER INDSHC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the validity of arguments .....
C No need to check NPHPR is power of 2 - this is done
C by FFTRLV ...
C ......... need to have MMAX < NPHPR/2
C
      IF ( NTHPTS.LE.LH ) THEN
         PRINT *,' Subroutine GRINVZ.'
         PRINT *,' NTHPTS must be atleast ', LH + 1
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPR.LE.(2*MMAX+1) ) THEN
         PRINT *,' Subroutine GRINVZ.'
         PRINT *,' NPHPR must be atleast ',(2*MMAX+1)
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPV*M0.NE.NPHPR ) THEN
         PRINT *,' Subroutine GRINVZ.'
         PRINT *,' M0    = ', M0
         PRINT *,' NTHPV = ', NPHPV
         PRINT *,' NTHPR = ', NPHPR
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
C Pre-empt a division by zero with RAD
C
      IF ( ABS( RAD ).LT.TOL ) THEN
         PRINT *,' Subroutine GRINVZ.'
         PRINT *,' RAD = ', RAD,' Division by zero imminent.'
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IOP = 0
      ISIGN = -1
C ............. set GRFNZ array to zero
      CALL CUBEOP ( GRFNZ, ZERO, ITHREE, NPHPV, NTHPTS, IOP )
C
C ............. start to loop around theta points ...................
      DO ITHETA = 1, NTHPTS
         X = GAUX( ITHETA )
         THETA = ACOS( X )
         SINE = DSIN( THETA )
C
C ............. set FTF1, FTF2, FTF3 all to zero
         LENV = 2*NPHPR
         CALL VECOP ( FTF1, ZERO, LENV, IOP )
         CALL VECOP ( FTF2, ZERO, LENV, IOP )
         CALL VECOP ( FTF3, ZERO, LENV, IOP )
C
C ............. start to loop around Harmonics ......................
C   						 ... case of M = 0 ..
         M = 0
         DO L = 1, LH
            IP = L*(L+1)/2+M+1
c           ICS = 1
c           IND = INDSHC( L, M, ICS )
            IND = L*L
            PSI = SHC( IND )
            DPSI = DSHC( IND )
            FTF1( 1 ) = FTF1( 1 ) + DPSI*PA( IP, ITHETA )
            FTF2( 1 ) = FTF2( 1 ) + PSI*DPA( IP, ITHETA )/RAD
         ENDDO
C					     ... end of case M = 0 ..
         DO M = 1, MIN( LH, MMAX )
            IF ( MOD( M, M0 ).NE.0 ) GOTO 67
            IND1 = 2*M + 1
            IND2 = 2*M + 2
            DO L = M, LH
               IP = L*(L+1)/2+M+1
c              ICS = 1
c              IND = INDSHC( L, M, ICS )
               IND = L*L + 2*M - 1
               PSI = SHC( IND )
               DPSI = DSHC( IND )
c              ICS = 2
c              IND = INDSHC( L, M, ICS )
               IND = L*L + 2*M
               PSI2 = SHC( IND )
               DPSI2 = DSHC( IND )
C
C ........( a little note ... here PSI = PSI_l^(mc) and PSI2 is
C ......... of course PSI_l^(ms) .............................)
C
C              ................................. radial component
               FTF1( IND1 ) = FTF1( IND1 ) + DPSI*PA(IP,ITHETA)
               FTF1( IND2 ) = FTF1( IND2 ) + DPSI2*PA(IP,ITHETA)
C              ................................. theta component
               FTF2( IND1 ) = FTF2( IND1 ) +
     1          DPA( IP, ITHETA )*PSI/RAD
               FTF2( IND2 ) = FTF2( IND2 ) +
     1          DPA( IP, ITHETA )*PSI2/RAD
C              ................................. phi component
               FTF3( IND1 ) = FTF3( IND1 ) +
     1          PSI2*PA( IP, ITHETA )*M/(RAD*SINE)
               FTF3( IND2 ) = FTF3( IND2 ) -
     1          PSI*PA( IP, ITHETA )*M/(RAD*SINE)
C 
            ENDDO
 67         CONTINUE
         ENDDO
C ............. ended looping around Harmonics ......................
C 
C ............... now perform Fourier Transforms on FTF1, FTF2, FTF3
         CALL FFTRLV ( FTF1, NPHPR, ISIGN )
         CALL FFTRLV ( FTF2, NPHPR, ISIGN )
         CALL FFTRLV ( FTF3, NPHPR, ISIGN )
C ...................................................................
         DO IPHI = 1, NPHPV
            GRFNZ( 1, IPHI, ITHETA ) = FTF1( 2*IPHI - 1 ) +
     1                                   DZCOEF
            GRFNZ( 2, IPHI, ITHETA ) = FTF2( 2*IPHI - 1 )
            GRFNZ( 3, IPHI, ITHETA ) = FTF3( 2*IPHI - 1 )
         ENDDO
C
      ENDDO
C ............. ended looping around theta points ...................
C
      RETURN
      END
C*********************************************************************
