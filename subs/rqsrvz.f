C*********************************************************************
C subroutine Radial QSt to Radial Vector function: aZimuthal sym. ****
C            -      ---    -      -                 -             ****
C Steve Gibbons Sat Jul 22 16:18:28 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
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
C     NR        : Number of radial grid nodes.                       C
C     ILNR      : Lowest radial node to be acted upon.               C
C     IRNR      : Highest radial node to be acted upon.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C                                                                    C
C     RQSTO     : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector.                          C
C                  Has dimensions ( NR, LH*(LH+2), 3 ).              C
C              RQSTO (I,l*l+2m,1) = q_l^ms(r_i)                      C
C   RQSTO (I,l*l+2m-1+delta_m0,1) = q_l^mc(r_i)                      C
C              RQSTO (I,l*l+2m,2) = s_l^ms(r_i)                      C
C   RQSTO (I,l*l+2m-1+delta_m0,2) = s_l^mc(r_i)                      C
C              RQSTO (I,l*l+2m,3) = t_l^ms(r_i)                      C
C   RQSTO (I,l*l+2m-1+delta_m0,3) = t_l^mc(r_i)                      C
C                                                                    C
C     ZCOEFA    : Coefficients of the monopole term in the R compon. C
C                 This is not necessarily zero in general.           C
C                 Dimension ( NR )                                   C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RVFZ      : Vector Function. An array of dimensions            C
C                ( 3, NPHPV, NTHPTS, NR) which contain the R, THETA  C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. RVFZ( 2, IPHI, ITHETA, NR ) is the Theta C
C                  compontent of the vector at (ir, iphi, itheta).   C
C                                                                    C
C                 RVFZ is completely zeroed on entry to RQSRVZ.      C
C____________________________________________________________________C
C                                                                    C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPR )                        C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RQSRVZ( RQSTO, RVFZ, GAUX, PA, DPA, FTF1, FTF2,
     1                   FTF3, LH, NTHPTS, NPHPR, NPHPV, M0, MMAX,
     2                   ZCOEFA, NR, ILNR, IRNR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS, NPHPV, NPHPR, M0, MMAX, NR, ILNR, IRNR
      DOUBLE PRECISION RQSTO( NR, LH*(LH+2), 3 ), ZCOEFA( NR ),
     1                 RVFZ( 3, NPHPV, NTHPTS, NR), GAUX( NTHPTS )
      DOUBLE PRECISION FTF1( 2*NPHPR ),
     1                 FTF2( 2*NPHPR ), FTF3( 2*NPHPR ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, ISIGN, ITHREE, ITHETA, LENV, IPHI, IR, 
     1        IP, IND1, IND2,  L, M, INDFTC, INDFTS
      DOUBLE PRECISION ZERO,X,SINE,TERM1, TERM2
      PARAMETER ( ZERO = 0.0d0 , ITHREE=3 )
C____________________________________________________________________C
C Functions used :-
c     INTEGER INDSHC
      DOUBLE PRECISION SQRLL1
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
         PRINT *,' Subroutine RQSRVZ.'
         PRINT *,' NTHPTS must be atleast ', LH + 1
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPR.LE.(2*MMAX+1) ) THEN
         PRINT *,' Subroutine RQSRVZ.'
         PRINT *,' NPHPR must be atleast ',(2*MMAX+1)
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPV*M0.NE.NPHPR ) THEN
         PRINT *,' Subroutine RQSRVZ.'
         PRINT *,' M0    = ', M0
         PRINT *,' NTHPV = ', NPHPV
         PRINT *,' NTHPR = ', NPHPR
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IOP = 0
C ............. set RVFZ array to zero
      CALL QUADOP( RVFZ, ZERO, ITHREE, NPHPV, NTHPTS, NR, IOP)
C
C     .................. start loop around radial grid nodes
      DO IR = ILNR, IRNR
C
C ............. start to loop around theta points ...................
        DO ITHETA = 1, NTHPTS
C
         X = GAUX( ITHETA )
         SINE = DSIN( ACOS( X ) )
C
C ............. set FTF1, FTF2, FTF3 all to zero
         LENV = 2*NPHPR
C
         CALL VECOP ( FTF1, ZERO, LENV, IOP )
         CALL VECOP ( FTF2, ZERO, LENV, IOP )
         CALL VECOP ( FTF3, ZERO, LENV, IOP )
C
C ............. start to loop around Harmonics ......................
         
         M = 0
         DO L = 1, LH
            IP = L*(L+1)/2+M+1
            FTF1( 1 ) = FTF1( 1 ) 
     1                   + PA( IP, ITHETA )*RQSTO( IR, L*L, 1 )
            FTF2( 1 ) = FTF2( 1 ) 
     1            + DPA( IP, ITHETA )*RQSTO( IR, L*L, 2 )/SQRLL1(L)
            FTF3( 1 ) = FTF3( 1 ) 
     1            + DPA( IP, ITHETA )*RQSTO( IR, L*L, 3 )/SQRLL1(L)
         ENDDO
C ................. end of case M = 0
         DO M = 1, MIN( LH, MMAX )
            IF ( MOD( M, M0 ).NE.0 ) GOTO 67
            INDFTC = 2*M + 1
            INDFTS = 2*M + 2
            DO L = M, LH
               IP = L*(L+1)/2+M+1
c              ICS = 1
c              IND1 = INDSHC( L, M, ICS )
               IND1 = L*L + 2*M - 1
c              ICS = 2
c              IND2 = INDSHC( L, M, ICS )
               IND2 = L*L + 2*M
               TERM1 = M*PA( IP , ITHETA )/( SINE*SQRLL1( L ) )
               TERM2 = DPA( IP , ITHETA )/SQRLL1( L )
C
               FTF1( INDFTC ) = FTF1( INDFTC ) +
     1          PA( IP, ITHETA )*RQSTO( IR, IND1, 1 )
               FTF1( INDFTS ) = FTF1( INDFTS ) +
     1          PA( IP, ITHETA )*RQSTO( IR, IND2, 1 )
               FTF2( INDFTC ) = FTF2( INDFTC ) +
     1          TERM2*RQSTO( IR, IND1, 2 ) - TERM1*RQSTO( IR, IND2, 3 )
               FTF2( INDFTS ) = FTF2( INDFTS ) +
     1          TERM2*RQSTO( IR, IND2, 2 ) + TERM1*RQSTO( IR, IND1, 3 )
               FTF3( INDFTC ) = FTF3( INDFTC ) +
     1          TERM2*RQSTO( IR, IND1, 3 ) + TERM1*RQSTO( IR, IND2, 2 )
               FTF3( INDFTS ) = FTF3( INDFTS ) +
     1          TERM2*RQSTO( IR, IND2, 3 ) - TERM1*RQSTO( IR, IND1, 2 )
C
            ENDDO
 67         CONTINUE
         ENDDO
C ............. ended looping around Harmonics ......................
C
C ............... now perform Fourier Transforms on FTF1, FTF2, FTF3
C
         ISIGN = -1
         CALL FFTRLV ( FTF1, NPHPR, ISIGN )
         CALL FFTRLV ( FTF2, NPHPR, ISIGN )
         CALL FFTRLV ( FTF3, NPHPR, ISIGN )
C
C ...................................................................
         DO IPHI = 1, NPHPV
          RVFZ( 1, IPHI, ITHETA, IR) = FTF1( 2*IPHI - 1 ) + ZCOEFA(IR)
          RVFZ( 2, IPHI, ITHETA, IR) = FTF2( 2*IPHI - 1 )
          RVFZ( 3, IPHI, ITHETA, IR) = FTF3( 2*IPHI - 1 )
         ENDDO
C
        ENDDO
C ............. ended looping around theta points ...................
      ENDDO
C     . Ended loop around radial grid nodes
C
      RETURN
      END
C*********************************************************************

