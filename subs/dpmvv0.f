C*********************************************************************
C subroutine Double Precision Matrix curl ( Velocity . grad Vel_0 ) **
C            -      -         -             -               -   -   **
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C                                                                    C
C This routine adds to the convection matrix the terms associated    C
C with the non-linear v . grad v_0.                                  C
C                                                                    C
C This routine assumes that the harmonic sets which constitute       C
C the old solution vector and the new are non-identical.             C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     LVEC      : Length of vector containing old solution.          C
C     NR        : Number of radial grid nodes.                       C
C     NOH       : Number of harmonics in old solution.               C
C     NNH       : Number of harmonics in new solution.               C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C                                                                    C
C     IMF       : Matrix format flag.                                C
C                                                                    C
C          imf = 1; Matrix is in LAPACK banded format                C
C                   ie element a_{i,j} is stored in                  C
C                   A( kle + ku + 1 + i - j , j )                    C
C                                                                    C
C          imf = 2; Matrix is banded but with element a_{i,j}        C
C                   stored in A( kl + 1 + j - i , i ).               C
C                                                                    C
C          imf = 3; Matrix is square - ie a_{i,j} is stored          C
C                   in A( i, j ).                                    C
C                                                                    C
C     MOHT       : Dimension ( NOH )                                 C
C     MOHL       : Dimension ( NOH )                                 C
C     MOHM       : Dimension ( NOH )                                 C
C     MOHC       : Dimension ( NOH )                                 C
C                                                                    C
C  moht, mohl, mohm and mohc define the spherical harmonics present  C
C  in the old solution vector. For spherical harmonic number I;      C
C                                                                    C
C   MOHT( I ) = 1 for a poloidal velocity vector                     C
C   MOHT( I ) = 2 for a toroidal velocity vector                     C
C   MOHT( I ) = 3 for a temperature / codensity term                 C
C                                                                    C
C   MOHL( I ) = spherical harmonic degree, l                         C
C                                                                    C
C   MOHM( I ) = spherical harmonic order, m                          C
C                                                                    C
C   MOHC( I ) = 1 for a cosine dependence in phi and                 C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C     MNHT       : Dimension ( NNH )                                 C
C     MNHL       : Dimension ( NNH )                                 C
C     MNHM       : Dimension ( NNH )                                 C
C     MNHC       : Dimension ( NNH )                                 C
C                                                                    C
C  mnht, mnhl, mnhm and mnhc define the spherical harmonics present  C
C  in the new solution vector. For spherical harmonic number I;      C
C                                                                    C
C   MNHT( I ) = 1 for a poloidal velocity vector                     C
C   MNHT( I ) = 2 for a toroidal velocity vector                     C
C   MNHT( I ) = 3 for a temperature / codensity term                 C
C                                                                    C
C   MNHL( I ) = spherical harmonic degree, l                         C
C                                                                    C
C   MNHM( I ) = spherical harmonic order, m                          C
C                                                                    C
C   MNHC( I ) = 1 for a cosine dependence in phi and                 C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C     LH        : Level of harmonics.                                C
C     NPHPTS    : Number of phi points.                              C
C     NTHPTS    : Number of theta points.                            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     AMAT      : Matrix. Dimensions ( N1, N2 )                      C
C              Will generally be banded due to the nature of the     C
C              numerical scheme. KL, KU and KLE parameterise this.   C
C                                                                    C
C     FACVV0    : Multiplication factor for curl ( v . grad ) v_0    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2,NTHPTS)       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C     OLDVEC    : Vector containing (v_0, theta_0 ). Dimen. ( N2 )   C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     ORD       : *(2). Determines the order of accuracy. Options :  C
C             SS - Strictly second order                             C
C             SF - Strictly fourth order                             C
C             O5 - Optimum accuracy for bandwidth 5; this gives      C
C                  Fourth order accuracy for 1st and 2nd derivatives C
C                  and second order accuracy for 3rd and 4th der.s   C
C             O7 - Optimum accuracy for bandwidth 7; this gives      C
C                  Sixth order accuracy for 1st and 2nd derivatives  C
C                  and fourth order accuracy for 3rd and 4th der.s   C
C                                                                    C
C____________________________________________________________________C
C Working Arrays   :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF1      : Array for Fourier transforming dim. (2*NPHPTS)     C
C     FTF2      : Array for Fourier transforming dim. (2*NPHPTS)     C
C     FTF3      : Array for Fourier transforming dim. (2*NPHPTS)     C
C     VFA       : Vector function. dimensions ( NPHPTS, NTHPTS, 3)   C
C     VFB       : Vector function. dimensions ( NPHPTS, NTHPTS, 3)   C
C     VFG       : Vector function. dimensions ( NPHPTS, NTHPTS, 3)   C
C     QST       : Obvious ... ( dimension  (  LH*( LH+2 ) , 3 )      C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPMVV0 ( N1, N2, LVEC, NR, NOH, NNH, KL, KU, KLE,
     1                   IMF, MOHT, MOHL, MOHM, MOHC, MNHT, MNHL,
     2                   MNHM, MNHC, LH, NPHPTS, NTHPTS, RI, RO, AMAT,
     3                   FACVV0, GAUX, GAUW, PA, DPA,
     4                   OLDVEC, ORD, FTF1, FTF2, FTF3, VFA, VFB,
     5                   VFG, QST )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, LVEC, NR, NOH, NNH, KL, KU, KLE, IMF,
     1        MOHT( NOH ), MOHM( NOH ), MOHC( NOH ), MOHL( NOH ),
     2        MNHT( NNH ), MNHM( NNH ), MNHC( NNH ), MNHL( NNH ),
     3        LH, NTHPTS, NPHPTS
      DOUBLE PRECISION RI, RO, AMAT( N1, N2 ), FACVV0,
     1                 GAUX( NTHPTS ), GAUW( NTHPTS ),
     2                 PA( (LH+1)*(LH+2)/2 , NTHPTS ),
     3                 DPA( (LH+1)*(LH+2)/2 , NTHPTS )
      DOUBLE PRECISION FTF1( 2*NPHPTS ),
     1                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ),
     2                 VFA( NPHPTS, NTHPTS, 3 ),
     3                 VFB( NPHPTS, NTHPTS, 3 ),
     4                 VFG( NPHPTS, NTHPTS, 3 ),
     5                 QST( LH*(LH+2), 3 )
      DOUBLE PRECISION OLDVEC( N2 )
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TOL, FAC, DPPARS( 1 ),
     1                 QABG, SABG, TABG
      EXTERNAL NLITS
      PARAMETER ( TOL = 1.0d-7 )
      INTEGER IALPHA, LA, MA, ICSA, IHMA, 
     1        IBETA, LB, MB, ICSB, IHMB, 
     2        IGAMMA, LG, MG, ICSG, IHMG, 
     3        INTPAR( 8 ), INDSHC, IFN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of N2 against NH and NR etc.
C
      IF ( N2.NE.NNH*NR ) THEN
        PRINT *,' Subroutine DPMVV0.  bad array size'
        PRINT *,' N2 (second array dimension) = ',N2
        PRINT *,' NNH = ',NNH,' NR = ',NR
        STOP
      ENDIF
C
C Check value of N2 against NH and NR etc.
C
      IF ( LVEC.NE.NOH*NR ) THEN
        PRINT *,' Subroutine DPMVV0.  bad vector length'
        PRINT *,' LVEC = ',LVEC
        PRINT *,' NOH = ',NOH,' NR = ',NR
        STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FACVV0 ).LT.TOL ) RETURN
C
C     --------------------------------------------------------
c     First let's loop around the toroidal harmonics IALPHA
      DO IALPHA = 1, NNH
        IF ( MNHT( IALPHA ).EQ.2 ) THEN
         LA   = MNHL ( IALPHA )
         MA   = MNHM ( IALPHA )
         ICSA = MNHC ( IALPHA )
         IHMA = INDSHC ( LA, MA, ICSA )
c        .
c        .           So work out the toroidal ialpha
c        .           vector in space in the function VFA
c        .
         IFN = 3
         CALL SHVECT ( IHMA, IFN, VFA, GAUX, PA, DPA, 
     1                 NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c        . Now we need to loop around all the vector spherical
c        . harmonics in the velocity. These will be denoted
c        . by the loop IBETA and we ofcourse need to treat the
c        . toroidal and poloidal cases separately.
c        .
c        . First do the poloidal harmonics
c        . ( and so we need to transform the equivalent
c        . toroidal vector spherical harmonic because of
c        . the curl of the 'beta' harmonic. )
c        .
         DO IBETA = 1, NOH
           IF ( MOHT( IBETA ).EQ.1 ) THEN
            LB   = MOHL ( IBETA )
            MB   = MOHM ( IBETA )
            ICSB = MOHC ( IBETA )
            IHMB = INDSHC( LB, MB, ICSB )
c           .
c           . So work out toroidal vector ibeta
c           . in space in function VFB
c           .
            IFN = 3
            CALL SHVECT ( IHMB, IFN, VFB, GAUX, PA, DPA, 
     1                NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c           .
c           . VFA contains the vector t_a
c           . VFB contains the vector t_b
c           .
            CALL VFCP ( VFA, VFB, VFG, NPHPTS, NTHPTS )
c           .
c           . VFG now contains the vector t_a x t_b
c           . Now let's transform into QST coefficients
c           .
            CALL VF2QSA ( QST, VFG, GAUX, GAUW, PA, DPA, 
     1                    FTF1, FTF2, FTF3, LH, LH,
     2                    NTHPTS, NTHPTS, NPHPTS, NPHPTS )
c           .
c           . QST( IHMG , 1 ) should now contain Q_abg^TT
c           . QST( IHMG , 2 ) should be zero.
c           . QST( IHMG , 3 ) should be zero.
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.1 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               QABG = QST( IHMG , 1 )
               IF ( ABS( QABG ).GT.TOL ) THEN
C
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new toroidal alpha, old poloidal beta
c        . output is toroidal gamma
c        .
            INTPAR( 1 ) = 13
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = QABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, ntorh loop ( Q_abg^TT terms )
           ENDIF
         ENDDO
c        .   end of ibeta = 1, npolh loop for toroidal ialpha
c        .
c        . Now do the toroidal harmonics
c        . ( So transform over the equivalent poloidal
c        . harmonics )
c        .
         DO IBETA = 1, NOH
           IF ( MOHT( IBETA ).EQ.2 ) THEN
            LB   = MOHL ( IBETA )
            MB   = MOHM ( IBETA )
            ICSB = MOHC ( IBETA )
            IHMB = INDSHC( LB, MB, ICSB )
c           .
c           . So work out scaloidal vector ibeta
c           . in space in function VFB
c           .
            IFN = 1
            CALL SHVECT ( IHMB, IFN, VFB, GAUX, PA, DPA, 
     1                NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c           .
c           . VFA contains the vector t_a
c           . VFB contains the vector q_b
c           .
            CALL VFCP ( VFA, VFB, VFG, NPHPTS, NTHPTS )
c           .
c           . VFG now contains the vector t_a x q_b
c           . Now let's transform into QST coefficients
c           .
            CALL VF2QSA ( QST, VFG, GAUX, GAUW, PA, DPA, 
     1                    FTF1, FTF2, FTF3, LH, LH,
     2                    NTHPTS, NTHPTS, NPHPTS, NPHPTS )
c           .
c           . QST( IHMG , 1 ) should be zero.
c           . QST( IHMG , 2 ) should contain S_abg^TQ
c           . QST( IHMG , 3 ) should contain T_abg^TQ
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.1 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               SABG = QST( IHMG , 2 )
               IF ( ABS( SABG ).GT.TOL ) THEN
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new toroidal alpha, old toroidal beta
c        . output is toroidal gamma
c        .
            INTPAR( 1 ) = 14
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = SABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, ntorh loop ( S_abg^TQ terms )
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.2 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               TABG = QST( IHMG , 3 )
               IF ( ABS( TABG ).GT.TOL ) THEN

c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new toroidal alpha, old toroidal beta
c        . output is poloidal gamma
c        .
            INTPAR( 1 ) = 15
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = TABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, npolh loop ( T_abg^TQ terms )

c           .
c           . So work out spheroidal vector ibeta
c           . in space in function VFB
c           .
            IFN = 2
            CALL SHVECT ( IHMB, IFN, VFB, GAUX, PA, DPA, 
     1                NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c           .
c           . VFA contains the vector t_a
c           . VFB contains the vector s_b
c           .
            CALL VFCP ( VFA, VFB, VFG, NPHPTS, NTHPTS )
c           .
c           . VFG now contains the vector t_a x s_b
c           . Now let's transform into QST coefficients
c           .
            CALL VF2QSA ( QST, VFG, GAUX, GAUW, PA, DPA, 
     1                    FTF1, FTF2, FTF3, LH, LH,
     2                    NTHPTS, NTHPTS, NPHPTS, NPHPTS )
c           .
c           . QST( IHMG , 1 ) should now contain Q_abg^TS
c           . QST( IHMG , 2 ) should be zero.
c           . QST( IHMG , 3 ) should be zero.
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.1 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               QABG = QST( IHMG , 1 )
               IF ( ABS( QABG ).GT.TOL ) THEN
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new toroidal alpha, old toroidal beta
c        . output is toroidal gamma
c        .
            INTPAR( 1 ) = 16
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = QABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, ntorh loop ( Q_abg^TS terms )
           ENDIF
         ENDDO
c        .   end of ibeta = 1, ntorh loop for toroidal ialpha
       ENDIF
      ENDDO
c     .   end of ialpha = 1, nh loop for toroidal harmonics
C     -----------------------------------------------------
c     Now let's loop around the poloidal harmonics IALPHA
      DO IALPHA = 1, NNH
        IF ( MNHT( IALPHA ).EQ.1 ) THEN
         LA   = MNHL ( IALPHA )
         MA   = MNHM ( IALPHA )
         ICSA = MNHC ( IALPHA )
         IHMA = INDSHC ( LA, MA, ICSA )
c        .
c        .           So work out the scaloidal ialpha
c        .           vector in space in the function VFA
c        .
         IFN = 1
         CALL SHVECT ( IHMA, IFN, VFA, GAUX, PA, DPA, 
     1                 NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c        . Now we need to loop around all the vector spherical
c        . harmonics in the velocity. These will be denoted
c        . by the loop IBETA and we ofcourse need to treat the
c        . toroidal and poloidal cases separately.
c        .
c        . Loop around poloidal harmonics beta and so perform
c        . transforms on the equivalent toroidal harmonics
c        .
         DO IBETA = 1, NOH
           IF ( MOHT( IBETA ).EQ.1 ) THEN
            LB   = MOHL ( IBETA )
            MB   = MOHM ( IBETA )
            ICSB = MOHC ( IBETA )
            IHMB = INDSHC( LB, MB, ICSB )
c           .
c           . So work out toroidal vector ibeta
c           . in space in function VFB
c           .
            IFN = 3
            CALL SHVECT ( IHMB, IFN, VFB, GAUX, PA, DPA, 
     1                NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c           .
c           . VFA contains the vector q_a
c           . VFB contains the vector t_b
c           .
            CALL VFCP ( VFA, VFB, VFG, NPHPTS, NTHPTS )
c           .
c           . VFG now contains the vector q_a x t_b
c           . Now let's transform into QST coefficients
c           .
            CALL VF2QSA ( QST, VFG, GAUX, GAUW, PA, DPA, 
     1                    FTF1, FTF2, FTF3, LH, LH,
     2                    NTHPTS, NTHPTS, NPHPTS, NPHPTS )
c           .
c           . QST( IHMG , 1 ) should be zero.
c           . QST( IHMG , 2 ) should now contain S_abg^QT
c           . QST( IHMG , 3 ) should now contain T_abg^QT
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.1 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               SABG = QST( IHMG , 2 )
               IF ( ABS( SABG ).GT.TOL ) THEN
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new poloidal alpha, old poloidal beta
c        . output is toroidal gamma
c        .
            INTPAR( 1 ) = 17
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = SABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, ntorh loop ( S_abg^QT terms )
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.2 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               TABG = QST( IHMG , 3 )
               IF ( ABS( TABG ).GT.TOL ) THEN
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new poloidal alpha, old poloidal beta
c        . output is poloidal gamma
c        .
            INTPAR( 1 ) = 18
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = TABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, npolh loop ( T_abg^QT terms )
           ENDIF
         ENDDO
c        .    End loop ibeta = 1, npolh for scaloidal ialpha

         DO IBETA = 1, NOH
           IF ( MOHT( IBETA ).EQ.2 ) THEN
            LB   = MOHL ( IBETA )
            MB   = MOHM ( IBETA )
            ICSB = MOHC ( IBETA )
            IHMB = INDSHC( LB, MB, ICSB )
c           .
c           . So work out spheroidal vector ibeta
c           . in space in function VFB
c           .
            IFN = 2
            CALL SHVECT ( IHMB, IFN, VFB, GAUX, PA, DPA, 
     1                NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c           .
c           . VFA contains the vector q_a
c           . VFB contains the vector s_b
c           .
            CALL VFCP ( VFA, VFB, VFG, NPHPTS, NTHPTS )
c           .
c           . VFG now contains the vector q_a x s_b
c           . Now let's transform into QST coefficients
c           .
            CALL VF2QSA ( QST, VFG, GAUX, GAUW, PA, DPA, 
     1                    FTF1, FTF2, FTF3, LH, LH,
     2                    NTHPTS, NTHPTS, NPHPTS, NPHPTS )
c           .
c           . QST( IHMG , 1 ) should be zero.
c           . QST( IHMG , 2 ) should contain S_abg^QS
c           . QST( IHMG , 3 ) should contain T_abg^QS
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.1 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               SABG = QST( IHMG , 2 )
               IF ( ABS( SABG ).GT.TOL ) THEN
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new poloidal alpha, old toroidal beta
c        . output is toroidal gamma
c        .
            INTPAR( 1 ) = 19
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = SABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, ntorh loop ( S_abg^QS terms )
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.2 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               TABG = QST( IHMG , 3 )
               IF ( ABS( TABG ).GT.TOL ) THEN
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new poloidal alpha, old toroidal beta
c        . output is poloidal gamma
c        .
            INTPAR( 1 ) = 20
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = TABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, npolh loop ( T_abg^QS terms )
          ENDIF
         ENDDO
c        .    End loop ibeta = 1, ntorh for scaloidal ialpha
c        .
c        .           So work out the spheroidal ialpha
c        .           vector in space in the function VFA
c        .
         IFN = 2
         CALL SHVECT ( IHMA, IFN, VFA, GAUX, PA, DPA, 
     1                 NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c        . Now we need to loop around all the vector spherical
c        . harmonics in the velocity. These will be denoted
c        . by the loop IBETA and we ofcourse need to treat the
c        . toroidal and poloidal cases separately.
c
c        . Loop around the poloidal harmonics ibeta and so
c        . transform 'curled' i.e. toroidal vectors
c        .
         DO IBETA = 1, NOH
           IF ( MOHT( IBETA ).EQ.1 ) THEN
            LB   = MOHL ( IBETA )
            MB   = MOHM ( IBETA )
            ICSB = MOHC ( IBETA )
            IHMB = INDSHC( LB, MB, ICSB )
c           .
c           . So work out toroidal vector ibeta
c           . in space in function VFB
c           .
            IFN = 3
            CALL SHVECT ( IHMB, IFN, VFB, GAUX, PA, DPA, 
     1                NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c           .
c           . VFA contains the vector s_a
c           . VFB contains the vector t_b
c           .
            CALL VFCP ( VFA, VFB, VFG, NPHPTS, NTHPTS )
c           .
c           . VFG now contains the vector s_a x t_b
c           .  Now let's transform into QST coefficients
c           .
            CALL VF2QSA ( QST, VFG, GAUX, GAUW, PA, DPA, 
     1                    FTF1, FTF2, FTF3, LH, LH,
     2                    NTHPTS, NTHPTS, NPHPTS, NPHPTS )
c           .
c           . QST( IHMG , 1 ) should contain Q_abg^ST
c           . QST( IHMG , 2 ) should be zero.
c           . QST( IHMG , 3 ) should be zero.
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.1 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               QABG = QST( IHMG , 1 )
               IF ( ABS( QABG ).GT.TOL ) THEN
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new poloidal alpha, old poloidal beta
c        . output is toroidal gamma
c        .
            INTPAR( 1 ) = 21
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = QABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, ntorh loop ( Q_abg^ST terms )
           ENDIF
         ENDDO
c        .    End loop ibeta = 1, npolh for scaloidal ialpha

         DO IBETA = 1, NOH
           IF ( MOHT( IBETA ).EQ.2 ) THEN
            LB   = MOHL ( IBETA )
            MB   = MOHM ( IBETA )
            ICSB = MOHC ( IBETA )
            IHMB = INDSHC( LB, MB, ICSB )
c           .
c           . So work out scaloidal vector ibeta
c           . in space in function VFB
c           .
            IFN = 1
            CALL SHVECT ( IHMB, IFN, VFB, GAUX, PA, DPA, 
     1                NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c           .
c           . VFA contains the vector s_a
c           . VFB contains the vector q_b
c           .
            CALL VFCP ( VFA, VFB, VFG, NPHPTS, NTHPTS )
c           .
c           . VFG now contains the vector s_a x q_b
c           . Now let's transform into QST coefficients
c           .
            CALL VF2QSA ( QST, VFG, GAUX, GAUW, PA, DPA, 
     1                    FTF1, FTF2, FTF3, LH, LH,
     2                    NTHPTS, NTHPTS, NPHPTS, NPHPTS )
c           .
c           . QST( IHMG , 1 ) should be zero.
c           . QST( IHMG , 2 ) should contain S_abg^SQ
c           . QST( IHMG , 3 ) should contain T_abg^SQ
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.1 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               SABG = QST( IHMG , 2 )
               IF ( ABS( SABG ).GT.TOL ) THEN
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new poloidal alpha, old toroidal beta
c        . output is toroidal gamma
c        .
            INTPAR( 1 ) = 22
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = SABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, ntorh loop ( S_abg^SQ terms )
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.2 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               TABG = QST( IHMG , 3 )
               IF ( ABS( TABG ).GT.TOL ) THEN
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new poloidal alpha, old toroidal beta
c        . output is poloidal gamma
c        .
            INTPAR( 1 ) = 23
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = TABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, npolh loop ( T_abg^SQ terms )

c           .
c           . So work out spheroidal vector ibeta
c           . in space in function VFB
c           .
            IFN = 2
            CALL SHVECT ( IHMB, IFN, VFB, GAUX, PA, DPA, 
     1                NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c
c           .
c           . VFA contains the vector s_a
c           . VFB contains the vector s_b
c           .
            CALL VFCP ( VFA, VFB, VFG, NPHPTS, NTHPTS )
c           .
c           . VFG now contains the vector s_a x s_b
c           .  Now let's transform into QST coefficients
c           .
            CALL VF2QSA ( QST, VFG, GAUX, GAUW, PA, DPA, 
     1                    FTF1, FTF2, FTF3, LH, LH,
     2                    NTHPTS, NTHPTS, NPHPTS, NPHPTS )
c           .
c           . QST( IHMG , 1 ) should now contain Q_abg^SS
c           . QST( IHMG , 2 ) should be zero.
c           . QST( IHMG , 3 ) should be zero.
c           .
            DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.1 ) THEN
               LG   = MNHL ( IGAMMA )
               MG   = MNHM ( IGAMMA )
               ICSG = MNHC ( IGAMMA )
               IHMG = INDSHC ( LG, MG, ICSG )
               QABG = QST( IHMG , 1 )
               IF ( ABS( QABG ).GT.TOL ) THEN
c        .
c        . Calculate curl ( v_new . Nabla ) v_old
c        . new poloidal alpha, old toroidal beta
c        . output is toroidal gamma
c        .
            INTPAR( 1 ) = 24
            INTPAR( 2 ) = LA
            INTPAR( 3 ) = LB
            INTPAR( 4 ) = LG
            INTPAR( 5 ) = NR
            INTPAR( 6 ) = NOH
            INTPAR( 7 ) = IBETA
            INTPAR( 8 ) = LVEC
C
            DPPARS( 1 ) = QABG
C
            FAC = FACVV0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, NLITS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
               ENDIF
              ENDIF
            ENDDO
c           .    end of igamma = 1, ntorh loop ( Q_abg^SS terms )
           ENDIF
         ENDDO
c        .    End loop ibeta = 1, ntorh for spheroidal ialpha
       ENDIF
      ENDDO
c     .   end of ialpha = 1, nh loop for poloidal harmonics
C     -----------------------------------------------------
C
      RETURN
      END
C*********************************************************************
