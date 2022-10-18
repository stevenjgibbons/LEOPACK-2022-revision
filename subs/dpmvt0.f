C*********************************************************************
C subroutine Double Precision Matrix Velocity . grad Theta_0 *********
C            -      -         -      -               -     - *********
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C                                                                    C
C This routine adds to the convection matrix the terms associated    C
C with the non-linear v . grad theta_0.                              C
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
C     FACVT0    : Multiplication factor for v . grad (theta_0)       C
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
C     FTF      : Array for Fourier transforming dim. (2*NPHPTS)      C
C     VELVF     : Vector function. dimensions ( NPHPTS, NTHPTS, 3)   C
C     THEVF     : Vector function. dimensions ( NPHPTS, NTHPTS, 3)   C
C     SF        : Vector function. dimensions ( NTHPTS, NPHPTS )     C
C     SHC       : Sph. Harm. Coeff.s - dim ( LH*(LH + 2) )           C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPMVT0 ( N1, N2, LVEC, NR, NOH, NNH, KL, KU, KLE,
     1                   IMF, MOHT, MOHL, MOHM, MOHC, MNHT, MNHL,
     1                   MNHM, MNHC, LH, NPHPTS, NTHPTS, RI, RO,
     2                   AMAT, FACVT0, GAUX, GAUW, PA, DPA,
     3                   OLDVEC, ORD, FTF, VELVF, THEVF, SF, SHC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, LVEC, NR, NOH, NNH, KL, KU, KLE, IMF,
     1        MOHT( NOH ), MOHM( NOH ), MOHC( NOH ), MOHL( NOH ),
     2        MNHT( NNH ), MNHM( NNH ), MNHC( NNH ), MNHL( NNH ),
     3        LH, NTHPTS, NPHPTS
      DOUBLE PRECISION RI, RO, AMAT( N1, N2 ), FACVT0,
     1                 GAUX( NTHPTS ), GAUW( NTHPTS ),
     2                 PA( (LH+1)*(LH+2)/2 , NTHPTS ),
     3                 DPA( (LH+1)*(LH+2)/2 , NTHPTS )
      DOUBLE PRECISION FTF( 2*NPHPTS ), 
     1                 VELVF( NPHPTS, NTHPTS, 3 ),
     2                 THEVF( NPHPTS, NTHPTS, 3 )
      DOUBLE PRECISION OLDVEC( N2 ), SHC( LH*(LH+2) ),
     1                 SF( NTHPTS, NPHPTS )
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TOL, FAC, ZCOEF, DPPARS( 1 ),
     1                 QAB, SAB, TAB
      EXTERNAL VDGTTS
      PARAMETER ( TOL = 1.0d-7 )
      INTEGER IALPHA, LA, MA, ICSA,
     1        IBETA, LB, MB, ICSB,
     2        IGAMMA, LG, MG, ICSG,
     3        INTPAR( 7 ), IND1, IND2, IND3, INDSHC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of N2 against NNH and NR etc.
C
      IF ( N2.NE.NNH*NR ) THEN
        PRINT *,' Subroutine DPMVT0, bad array size'
        PRINT *,' N2 (second array dimension) = ',N2
        PRINT *,' NNH = ',NNH,' NR = ',NR
        STOP
      ENDIF
C
C Check value of LVEC against NOH and NR etc.
C
      IF ( LVEC.NE.NOH*NR ) THEN
        PRINT *,' Subroutine DPMVT0, bad vector length'
        PRINT *,' LVEC = ',LVEC
        PRINT *,' NOH = ',NOH,' NR = ',NR
        STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FACVT0 ).LT.TOL ) RETURN
C
      DO IBETA = 1, NOH
         IF (     MOHT( IBETA ).EQ.3      .AND.
     1            MOHL( IBETA ).EQ.0      .AND.
     2            MOHM( IBETA ).EQ.0      .AND.
     3            MOHC( IBETA ).EQ.1            ) THEN
C
C ...... Before we do any of the more complex interactions, let us
C        first look at the l=0 temperature harmonic. For this,
C        the non-linear term is simply
C
C        La ( La + 1 ) Pa( r ) / r * d/dr ( Theta_0^0 ( r ) )
C
         DO IALPHA = 1, NNH
           IF ( MNHT( IALPHA ).EQ.1 ) THEN
            LA = MNHL( IALPHA )
c                 .
C                 ...... Fill in the v_new.grad(Theta_old) term ...
c                 .
             INTPAR( 1 ) = 1
             INTPAR( 2 ) = LA
             INTPAR( 3 ) = 0
             INTPAR( 4 ) = NR
             INTPAR( 5 ) = NOH
             INTPAR( 6 ) = IBETA
             INTPAR( 7 ) = LVEC
C
            FAC = FACVT0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IALPHA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
           ENDIF
         ENDDO
        ENDIF
      ENDDO
C     .     end IBETA = 1, NH (for ibeta = theta_0^0 ....

c        .
C        ............................................................
C        . First loop around the NPOLH Theta harmonics ( IBETA )    .
C        ............................................................

      DO IBETA = 1, NOH
        IF ( MOHT( IBETA ).EQ.3 .AND. MOHL( IBETA ).NE.0 ) THEN
         LB   = MOHL( IBETA )
         MB   = MOHM( IBETA )
         ICSB = MOHC( IBETA )
         IND1 = INDSHC( LB, MB, ICSB )

c        .
C        .... First evaluate the Q_{gamma}^{ab} terms
c        .
         CALL SHVECT ( IND1, 1, THEVF, GAUX, PA, DPA,
     1                 NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )

C .......... so THEVF now contains the scaloidal vector harmonic
C            IBETA evaluated over the sphere. Now let's loop around
C            the poloidal velocity vectors IALPHA in order to do
C            the scalar products required ..........................

C         .
          DO IALPHA = 1, NNH
            IF ( MNHT( IALPHA ).EQ.1 ) THEN
             LA = MNHL( IALPHA )
             MA = MNHM( IALPHA )
             ICSA = MNHC( IALPHA )
             IND2 = INDSHC( LA, MA, ICSA )
             CALL SHVECT ( IND2, 1, VELVF, GAUX, PA, DPA,
     1                 NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
C
c            ..... so going through the IALPHA, let's take the
C            .     scalar products of VELVF and THEVF and
C            .     expand the resulting scalars in spherical
C            .     harmonics.
C
             CALL VFDP( THEVF, VELVF, SF, NPHPTS, NTHPTS )
C
             CALL FORSSA ( SHC, SF, GAUW, PA,
     1                     FTF, LH, LH, NTHPTS, NTHPTS,
     2                     NPHPTS, NPHPTS, ZCOEF )
C
C            ..... loop around scalar harms IGAMMA ...............
             DO IGAMMA = 1, NNH
              IF ( MNHT( IGAMMA ).EQ.3 ) THEN
               LG = MNHL( IGAMMA )
               IF ( LG.GT.0 ) THEN
                 MG   = MNHM( IGAMMA )
                 ICSG = MNHC( IGAMMA )
                 IND3 = INDSHC( LG, MG, ICSG )
                 QAB = SHC( IND3 )
               ELSE
                 QAB = ZCOEF
               ENDIF
               IF ( ABS(QAB).LT.TOL ) GOTO 701
c              .
C              .... We have a non-zero Q_{gamma}^{alpha beta} term
c              .
C              ...... Fill in the v_new.grad(fg+Theta_old) term ...
c              .

          INTPAR( 1 ) = 3
          INTPAR( 2 ) = LA
          INTPAR( 3 ) = 0
          INTPAR( 4 ) = NR
          INTPAR( 5 ) = NOH
          INTPAR( 6 ) = IBETA
          INTPAR( 7 ) = LVEC
C
          DPPARS( 1 ) = QAB
C
          FAC = FACVT0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
C......................................................................
C
 701         CONTINUE
             ENDIF
             ENDDO
C               ..... ended looping around scalar harms IGAMMA ......

           ENDIF
          ENDDO
C            ...... ended looping around scaloidal harmonics IALPHA..
C            .

          ENDIF
         ENDDO
C        ..... Ended looping around Theta harmonics .. ( IBETA ) ....

C        ..... Finished with IBETA scaloidal so let's begin to loop
C        ..... around IBETA spheroidal.
c        .
         DO IBETA = 1, NOH
           IF ( MOHT( IBETA ).EQ.3 .AND. MOHL( IBETA ).NE.0 ) THEN
            LB   = MOHL( IBETA )
            MB   = MOHM( IBETA )
            ICSB = MOHC( IBETA )
            IND1 = INDSHC( LB, MB, ICSB )

c           .
C           .........................................................
C           .                                                       .
C           .  Now for each of the temperature harmonics IBETA,     .
C           .  let's evaluate the spheroidal vector over the sphere .
C           .                                                       .
C           .  So let's loop around IALPHA spheroidal and toroidal  .
C           .  harmonics to calculate the S_gamma^{alpha beta}      .
C           .  and T_gamma^{alpha beta} terms.                      .
C           .                                                       .
C           .........................................................
c           .
            CALL SHVECT( IND1, 2, THEVF, GAUX, PA, DPA,
     1                   NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )
c           .
C           ..... THEVF now contains the spheroidal vector IBETA
C           ..... evaluated over the sphere.
c           .
C           ........ Loop around Poloidal IALPHA harmonics to
C           .        evaluate SPHEROIDAL vectors and take scalar
C           .        products.
            DO IALPHA = 1, NNH
              IF ( MNHT( IALPHA ).EQ.1 ) THEN
               LA   = MNHL( IALPHA )
               MA   = MNHM( IALPHA )
               ICSA = MNHC( IALPHA )
               IND2 = INDSHC( LA, MA, ICSA )
               CALL SHVECT ( IND2, 2, VELVF, GAUX, PA, DPA,
     1                   NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )

c               ..... so going through the IALPHA, let's take the
C               .     scalar products of VELVF and THEVF and
C               .     expand the resulting scalars in spherical
C               .     harmonics.
C
                CALL VFDP( THEVF, VELVF, SF, NPHPTS, NTHPTS )
C
                CALL FORSSA ( SHC, SF, GAUW, PA,
     1                        FTF, LH, LH, NTHPTS, NTHPTS,
     2                        NPHPTS, NPHPTS, ZCOEF )
C               ..... loop around scalar harms IGAMMA ...............
                DO IGAMMA = 1, NNH
                 IF ( MNHT( IGAMMA ).EQ.3 ) THEN
                  LG = MNHL( IGAMMA )
                  IF ( LG.GT.0 ) THEN
                    MG   = MNHM( IGAMMA )
                    ICSG = MNHC( IGAMMA )
                    IND3 = INDSHC( LG, MG, ICSG )
                    SAB = SHC( IND3 )
                  ELSE
                    SAB = ZCOEF
                  ENDIF
                  IF ( ABS(SAB).LT.TOL ) GOTO 711
c                 .
C                 .... We have a non-zero S_{gamma}^{alpha beta} term
c                 .
C                 ...... Fill in the v_new.grad(fg+Theta_old) term ...
c                 .
             INTPAR( 1 ) = 5
             INTPAR( 2 ) = LA
             INTPAR( 3 ) = LB
             INTPAR( 4 ) = NR
             INTPAR( 5 ) = NOH
             INTPAR( 6 ) = IBETA
             INTPAR( 7 ) = LVEC
C
             DPPARS( 1 ) = SAB
C
             FAC = FACVT0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
C......................................................................
C
 711            CONTINUE
                 ENDIF
                ENDDO
C               ..... ended looping around scalar harms IGAMMA ......
             ENDIF
            ENDDO
C           .
C           ..... End loop around IALPHA spheroidal ....


c           .
C           ........ Loop around Toroidal IALPHA harmonics to
C           .        evaluate TOROIDAL vectors and take scalar
C           .        products.
            DO IALPHA = 1, NNH
              IF ( MNHT( IALPHA ).EQ.2 ) THEN
               LA   = MNHL( IALPHA )
               MA   = MNHM( IALPHA )
               ICSA = MNHC( IALPHA )
               IND2 = INDSHC( LA, MA, ICSA )
               CALL SHVECT ( IND2, 3, VELVF, GAUX, PA, DPA,
     1                   NTHPTS, NTHPTS, NPHPTS, NPHPTS, LH )

c               ..... so going through the IALPHA, let's take the
C               .     scalar products of VELVF and THEVF and
C               .     expand the resulting scalars in spherical
C               .     harmonics.
C
                CALL VFDP( THEVF, VELVF, SF, NPHPTS, NTHPTS )
C
                CALL FORSSA ( SHC, SF, GAUW, PA,
     1                        FTF, LH, LH, NTHPTS, NTHPTS,
     2                        NPHPTS, NPHPTS, ZCOEF )
C
C               ..... loop around scalar harms IGAMMA ...............
                DO IGAMMA = 1, NNH
                 IF ( MNHT( IGAMMA ).EQ.3 ) THEN
                  LG = MNHL( IGAMMA )
                  IF ( LG.GT.0 ) THEN
                    MG   = MNHM( IGAMMA )
                    ICSG = MNHC( IGAMMA )
                    IND3 = INDSHC( LG, MG, ICSG )
                    TAB = SHC( IND3 )
                  ELSE
                    TAB = ZCOEF
                  ENDIF
                  IF ( ABS(TAB).LT.TOL ) GOTO 721
c                 .
C                 .... We have a non-zero T_{gamma}^{alpha beta} term
c                 .
C                 ...... Fill in the v_new.grad(fg+Theta_old) term ...
c                 .
C .....................................................................
C
             INTPAR( 1 ) = 7
             INTPAR( 2 ) = LA
             INTPAR( 3 ) = LB
             INTPAR( 4 ) = NR
             INTPAR( 5 ) = NOH
             INTPAR( 6 ) = IBETA
             INTPAR( 7 ) = LVEC
C
             DPPARS( 1 ) = TAB
C
             FAC = FACVT0
C
            CALL DPMES ( N1, N2, NR, NNH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
C......................................................................
C
 721            CONTINUE
                 ENDIF
                ENDDO
C               ..... ended looping around scalar harms IGAMMA ......
             ENDIF
            ENDDO
C           .
C           ..... End loop around IALPHA toroidal ....
C
C
          ENDIF
         ENDDO
C        ..... Ended looping around Theta harmonics .. ( IBETA ) ....
C
      RETURN
      END
C*********************************************************************
