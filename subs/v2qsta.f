C*********************************************************************
C subroutine Vector function TO QST coefficients (Array version) *****
C            -               -- ---               -              *****
C Steve Gibbons 26.4.97                                              C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     LH        : Level of harmonics.                                C
C     IRAD      : Current radial grid node.                          C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Corresponding weights. As above.                   C
C     VF        : Vector Function. An array of dimensions            C
C                  ( NPHPTS, NTHPTS, 3) which contain the R, THETA   C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IPHI, ITHETA, 2 ) is the Theta      C
C                  compontent of the vector at (iphi, itheta).       C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     QST       : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector c.f. eqn (38).            C
C                  Has dimensions (  LH*(LH+2) , 3, NR ).            C
C               QST (l*l+2m,1) = q_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,1) = q_l^mc(r_i)                         C
C                                (as in equation (39)                C
C               QST (l*l+2m,2) = s_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,2) = s_l^mc(r_i)                         C
C                                (as in equation (40),(41)           C
C               QST (l*l+2m,3) = t_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,3) = t_l^mc(r_i)                         C
C                                (as in equation (40),(41)           C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE V2QSTA ( QST, VF, GAUX, GAUW, PA, DPA, FTF1,
     1                    FTF2, FTF3, LH, NTHPTS, NPHPTS, IRAD, NR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS, NPHPTS, IRAD, NR
      DOUBLE PRECISION QST(  LH*(LH+2) , 3, NR ),
     1                 VF( NPHPTS, NTHPTS, 3),
     2                 GAUX( NTHPTS ),
     3                 GAUW( NTHPTS )
      DOUBLE PRECISION FTF1( 2*NPHPTS ),
     1                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ), 
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      LOGICAL POWT
      INTEGER IOP, ISIGN, ITHREE, ITHETA, 
     1        INDCOS, INDSIN, IP, IND1, L, M, IPHI
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
C First that NPHPTS is a power of 2.
c     CALL POWTWO ( NPHPTS, POWT )
c     IF ( .NOT. POWT ) THEN
c        PRINT *,' Subroutine V2QSTA.'
c        PRINT *,' NPHPTS is not a power of 2.'
c        PRINT *,' Program stopped.'
c        STOP
c     ENDIF
C ......... need to have LH < NPHPTS/2
      IF ( 2*LH.GT.NPHPTS ) THEN
         PRINT *,' Subroutine V2QSTA.'
         PRINT *,' NPHPTS must be at least 2*LH.'
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C .........
C
C ............ set array QST to zero .................................
      IOP = 0
      ISIGN = 1
      IND1 = LH*(LH+2)
      CALL MATOP ( QST, ZERO, IND1, ITHREE, IOP )
C ....................................................................
C ...... now start to loop around theta points .......................
      DO ITHETA = 1, NTHPTS
C
C
C******** Evaluate lengendre functions at given theta ****************
         X = GAUX( ITHETA )
         SINE = DSIN ( ACOS ( X ) )
         WEIGHT = GAUW( ITHETA )
C
C .................. firstly enter the VF information into FTF1, FTF2,
C .................. FTF3 for r, theta, phi components respectively of
C .................. VF.
C
         DO IPHI = 1, NPHPTS
            FTF1( 2*IPHI - 1 ) =
     1            2.0d0*VF( IPHI, ITHETA, 1 )/DBLE( NPHPTS )
            FTF1( 2*IPHI ) = ZERO
            FTF2( 2*IPHI - 1 ) =
     1            2.0d0*VF( IPHI, ITHETA, 2 )/DBLE( NPHPTS )
            FTF2( 2*IPHI ) = ZERO
            FTF3( 2*IPHI - 1 ) =
     1            2.0d0*VF( IPHI, ITHETA, 3 )/DBLE( NPHPTS )
            FTF3( 2*IPHI ) = ZERO
         ENDDO
C .................. now perform Forward Discreet Fourier Transforms
C .................. on FTF1, FTF2, FTF3
         CALL FFTNEW ( FTF1, NPHPTS, ISIGN )
         CALL FFTNEW ( FTF2, NPHPTS, ISIGN )
         CALL FFTNEW ( FTF3, NPHPTS, ISIGN )
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
            QST( INDCOS , 1, IRAD ) = QST( INDCOS , 1, IRAD) +
     1                          W1*PA( IP ,ITHETA )*FTF1( 1 )
            QST( INDCOS , 2, IRAD ) = QST( INDCOS , 2, IRAD) +
     1                          W2*DPA( IP ,ITHETA )*FTF2( 1 )
            QST( INDCOS , 3, IRAD ) = QST( INDCOS , 3, IRAD) +
     1                          W2*DPA( IP ,ITHETA )*FTF3( 1 )
C
            INDCOS = INDCOS - 1
C
C ___ Looping around from 1 to L ___ ********************************
C
C
            DO M = 1, L
               IP = IP + 1
               TERM = W2*DBLE(M)*PA( IP ,ITHETA )/SINE
               INDCOS = INDCOS + 2
               INDSIN = INDCOS + 1
C                                                   ... eqn 194
C                                                       Qlm cos
C***.Equation for Qlm COS ***************************
               QST( INDCOS , 1, IRAD ) = QST( INDCOS , 1, IRAD) +
C ..............contribution from B(rad)cos .........
     1          W1*PA( IP ,ITHETA )*FTF1( 2*M + 1)
C                                                   ... eqn 195
C                                                       Qlm sin
C***.Equation for Qlm SIN ***************************
               QST( INDSIN , 1, IRAD ) = QST( INDSIN , 1, IRAD) +
C ..............contribution from B(rad)sin .........
     1          W1*PA( IP ,ITHETA )*FTF1( 2*M + 2)
C                                                   ... eqn 203
C                                                       Slm cos
C***.Equation for Slm COS ***************************
               QST( INDCOS , 2, IRAD ) = QST( INDCOS , 2, IRAD) +
C ..............contribution from B(theta)cos .......
     1          W2*DPA( IP ,ITHETA )*FTF2( 2*M + 1 ) -
C ..............contribution from B(phi)sin .........
     2          TERM*FTF3( 2*M + 2 )
C                                                   ... eqn 204
C                                                       Slm sin
C***.Equation for Slm SIN ***************************
               QST( INDSIN , 2, IRAD ) = QST( INDSIN , 2, IRAD) +
C ..............contribution from B(theta)sin .......
     1          W2*DPA( IP ,ITHETA )*FTF2( 2*M + 2 ) +
C ..............contribution from B(phi)cos .........
     2          TERM*FTF3( 2*M + 1 )
C                                                   ... eqn 206
C                                                       Tlm cos
C***.Equation for Tlm COS ***************************
               QST( INDCOS , 3, IRAD ) = QST( INDCOS , 3, IRAD) +
C ..............contribution from B(phi)cos .........
     1          W2*DPA( IP ,ITHETA )*FTF3( 2*M + 1 ) +
C ..............contribution from B(theta)sin .......
     2          TERM*FTF2( 2*M + 2 )
C                                                   ... eqn 207
C                                                       Tlm sin
C***.Equation for Tlm SIN ***************************
               QST( INDSIN , 3, IRAD ) = QST( INDSIN , 3, IRAD) -
C ..............contribution from B(theta)cos .......
     1          TERM*FTF2( 2*M + 1 ) +
C ..............contribution from B(phi)sin .........
     2          W2*DPA( IP ,ITHETA )*FTF3( 2*M + 2 )
C
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

      RETURN
      END
C*********************************************************************
