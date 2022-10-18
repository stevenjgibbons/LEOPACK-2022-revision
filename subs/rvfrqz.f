C*********************************************************************
C subroutine Radial Vector Function to Radial Qst (aZimuthal sym.) ***
C            -      -      -           -      -     -              ***
C Steve Gibbons Sat Jul 22 13:03:26 BST 2000                         C
C____________________________________________________________________C
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
C     GAUW      : Corresponding weights. As above.                   C
C     RVFOZ     : Vector Function. An array of dimensions            C
C                ( 3, NPHPV, NTHPTS, NR) which contain the R, THETA  C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. RVFO( 2, IPHI, ITHETA, NR ) is the Theta C
C                  compontent of the vector at (ir, iphi, itheta).   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
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
C____________________________________________________________________C
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
      SUBROUTINE RVFRQZ ( RQSTO, RVFOZ, GAUX, GAUW, PA, DPA, FTF1,
     1                    FTF2, FTF3, ZCOEFA, LH, NTHPTS, NPHPV,
     2                    NPHPR, M0, MMAX, NR, ILNR, IRNR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS, NPHPV, NPHPR, M0, MMAX, NR, ILNR, IRNR
      DOUBLE PRECISION RQSTO( NR, LH*(LH+2), 3 ), ZCOEFA( NR ),
     1                 RVFOZ( 3, NPHPV, NTHPTS, NR ),
     2                 GAUX( NTHPTS ), GAUW( NTHPTS )
      DOUBLE PRECISION FTF1( 2*NPHPR ),
     1                 FTF2( 2*NPHPR ), FTF3( 2*NPHPR ), 
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, ISIGN, ITHREE, ITHETA, IR, 
     1        INDCOS, INDSIN, IP, IND1, L, M, IPHI, I1, I2
      DOUBLE PRECISION ZERO,X,SINE,TERM, WEIGHT, W1, W2
      PARAMETER ( ZERO = 0.0d0 , ITHREE=3 )
C____________________________________________________________________C
C Functions used :-
      DOUBLE PRECISION SQRLL1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the validity of arguments .....
C No need to check NPHPTS is power of 2 - this is done
C by FFTRLV ...
C ......... need to have MMAX < NPHPTS/2
C
      IF ( NTHPTS.LE.LH ) THEN
         PRINT *,' Subroutine RVFRQZ.'
         PRINT *,' NTHPTS must be atleast ', LH + 1
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPR.LE.(2*MMAX+1) ) THEN
         PRINT *,' Subroutine RVFRQZ.'
         PRINT *,' NPHPR must be atleast ',(2*MMAX+1)
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( ILNR.LT.1 .OR. IRNR.GT.NR .OR. ILNR.GT.IRNR ) THEN
         PRINT *,' Subroutine RVFRQZ.'
         PRINT *,' ILNR = ', ILNR,' IRNR = ', IRNR,' NR = ', NR
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPV*M0.NE.NPHPR ) THEN
         PRINT *,' Subroutine RVFRQZ.'
         PRINT *,' M0    = ', M0
         PRINT *,' NTHPV = ', NPHPV
         PRINT *,' NTHPR = ', NPHPR
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C_____________________________________________________________________
C
C
C ............ set array RQSTO at node IR to zero ....................
      IOP = 0
      IND1 = LH*(LH+2)
      DO I2 = 1, ITHREE
        DO I1 = 1, IND1
          DO IR = ILNR, IRNR
           RQSTO( IR, I1, I2 ) = ZERO
          ENDDO
        ENDDO
      ENDDO
C     .
C     . Loop around radial grid points
C     .
      DO IR = ILNR, IRNR
       ZCOEFA( IR ) = 0.0d0
C ....................................................................
C ...... now start to loop around theta points .......................
       DO ITHETA = 1, NTHPTS
C
C******** Evaluate lengendre functions at given theta ****************
C
         X = GAUX( ITHETA )
         SINE = DSIN ( ACOS ( X ) )
         WEIGHT = GAUW( ITHETA )
C        .
C        . Zero the arrays ftf1, ftf2 and ftf3
C        .
         IOP = 0
         IND1 = 2*NPHPR
         CALL VECOP( FTF1, ZERO, IND1, IOP )
         CALL VECOP( FTF2, ZERO, IND1, IOP )
         CALL VECOP( FTF3, ZERO, IND1, IOP )
C
C .................. firstly enter the RVF information into FTF1, FTF2,
C .................. FTF3 for r, theta, phi components respectively of
C .................. VF.
C
         DO IPHI = 1, NPHPV
            FTF1( 2*IPHI - 1 ) = RVFOZ( 1, IPHI, ITHETA, IR )
            FTF2( 2*IPHI - 1 ) = RVFOZ( 2, IPHI, ITHETA, IR )
            FTF3( 2*IPHI - 1 ) = RVFOZ( 3, IPHI, ITHETA, IR )
            ZCOEFA( IR ) = ZCOEFA( IR ) +
     1                   WEIGHT*RVFOZ( 1, IPHI, ITHETA, IR )
         ENDDO
C
         CALL FTAASC( NPHPR, NPHPV, FTF1 )
         CALL FTAASC( NPHPR, NPHPV, FTF2 )
         CALL FTAASC( NPHPR, NPHPV, FTF3 )
C
C .................. now perform Forward Discreet Fourier Transforms
C .................. on FTF1, FTF2, FTF3
         ISIGN = 1
         CALL FFTRLV ( FTF1, NPHPR, ISIGN )
         CALL FFTRLV ( FTF2, NPHPR, ISIGN )
         CALL FFTRLV ( FTF3, NPHPR, ISIGN )
C ...................................................................
C .................. Now let's loop around the Harmonics. ...........
C ...................................................................
         DO L = 1, LH
            W1 = WEIGHT*(2.0d0*L + 1.0d0)/4.0d0
            W2 = W1/SQRLL1( L )
            INDCOS = L*L
C .................. Now let's loop around the order, but considering
C .................. the Case M = 0 separately. 
C ___ CASE M = 0 ___ ************************************************
            M = 0
            IP = L*(L+1)/2+M+1
            RQSTO( IR, INDCOS, 1 ) = RQSTO( IR, INDCOS, 1 ) +
     1                          W1*PA( IP, ITHETA )*FTF1( 1 )
            RQSTO( IR, INDCOS, 2 ) = RQSTO( IR, INDCOS, 2 ) +
     1                          W2*DPA( IP, ITHETA )*FTF2( 1 )
            RQSTO( IR, INDCOS, 3 ) = RQSTO( IR, INDCOS, 3 ) +
     1                          W2*DPA( IP, ITHETA )*FTF3( 1 )
C
            INDCOS = INDCOS - 1
C
C ___ Looping around from 1 to L ___ ********************************
C
C
            DO M = 1, MIN( MMAX, L )
               IP = IP + 1
               TERM = W2*DBLE(M)*PA( IP, ITHETA )/SINE
               INDCOS = INDCOS + 2
               INDSIN = INDCOS + 1
               IF ( MOD( M, M0 ).NE.0 ) GOTO 67
C                                                   ... eqn 194
C                                                       Qlm cos
C***.Equation for Qlm COS ***************************
               RQSTO( IR, INDCOS, 1 ) = RQSTO( IR, INDCOS, 1 ) +
     1          W1*PA( IP, ITHETA )*FTF1( 2*M + 1)
C ..............contribution from B(rad)cos .........
C                                                   ... eqn 195
C                                                       Qlm sin
C***.Equation for Qlm SIN ***************************
               RQSTO( IR, INDSIN, 1 ) = RQSTO( IR, INDSIN, 1 ) +
     1          W1*PA( IP, ITHETA )*FTF1( 2*M + 2)
C ..............contribution from B(rad)sin .........
C                                                   ... eqn 203
C                                                       Slm cos
C***.Equation for Slm COS ***************************
               RQSTO( IR, INDCOS, 2 ) = RQSTO( IR, INDCOS, 2 ) +
     1          W2*DPA( IP, ITHETA )*FTF2( 2*M + 1 ) -
     2          TERM*FTF3( 2*M + 2 )
C ..............contribution from B(theta)cos .......
C ..............contribution from B(phi)sin .........
C                                                   ... eqn 204
C                                                       Slm sin
C***.Equation for Slm SIN ***************************
               RQSTO( IR, INDSIN, 2 ) = RQSTO( IR, INDSIN, 2 ) +
     1          W2*DPA( IP, ITHETA )*FTF2( 2*M + 2 ) +
     2          TERM*FTF3( 2*M + 1 )
C ..............contribution from B(theta)sin .......
C ..............contribution from B(phi)cos .........
C                                                   ... eqn 206
C                                                       Tlm cos
C***.Equation for Tlm COS ***************************
               RQSTO( IR, INDCOS, 3 ) = RQSTO( IR, INDCOS, 3 ) +
     1          W2*DPA( IP, ITHETA )*FTF3( 2*M + 1 ) +
     2          TERM*FTF2( 2*M + 2 )
C ..............contribution from B(phi)cos .........
C ..............contribution from B(theta)sin .......
C                                                   ... eqn 207
C                                                       Tlm sin
C***.Equation for Tlm SIN ***************************
               RQSTO( IR,  INDSIN, 3 ) = RQSTO( IR, INDSIN, 3 ) -
     1          TERM*FTF2( 2*M + 1 ) +
     2          W2*DPA( IP, ITHETA )*FTF3( 2*M + 2 )
C ..............contribution from B(theta)cos .......
C ..............contribution from B(phi)sin .........
C
 67            CONTINUE
            ENDDO
C
C
C Done looping around from 1 to L ___ *******************************
C
C
         ENDDO
C ...................................................................
C .................. Ended looping around the harmonics .............
C ...................................................................
       ENDDO
C ...... ended looping around theta points ...........................
C ....................................................................
       ZCOEFA( IR ) = 0.5d0*ZCOEFA( IR )/DBLE(NPHPV)
      ENDDO
C
C Ended looping around radial grid nodes ...
C
      RETURN
      END
C*********************************************************************
