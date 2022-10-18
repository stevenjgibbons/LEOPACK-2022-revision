C*********************************************************************
C subroutine Double Precision Matrix Velocity . Grad Theta ***********
C            -      -         -      -          -    -     ***********
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C                                                                    C
C This routine adds to the convection matrix the terms associated    C
C with the non-linear v . grad theta. This is the version which      C
C assumes that the sets of harmonics for v_0 and theta_0 are         C
C identical. This saves time in evaluating the integrals. It is      C
C however of limited use for linear stability problems.              C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of harmonics (all types)                    C
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
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C     MHM       : Dimension ( NH )                                   C
C     MHC       : Dimension ( NH )                                   C
C                                                                    C
C  mht, mhl, mhm and mhc define the spherical harmonics present      C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
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
C     FACV0T    : Multiplication factor for v_0 . grad (theta)       C
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
      SUBROUTINE DPMVGT ( N1, N2, NR, NH, KL, KU, KLE, IMF, MHT, MHL,
     1                   MHM, MHC, LH, NPHPTS, NTHPTS, RI, RO, AMAT,
     2                   FACV0T, FACVT0, GAUX, GAUW, PA, DPA,
     3                   OLDVEC, ORD, FTF, VELVF, THEVF, SF, SHC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, NR, NH, KL, KU, KLE, IMF,
     1        MHT( NH ), MHM( NH ), MHC( NH ), MHL( NH ),
     2        LH, NTHPTS, NPHPTS
      DOUBLE PRECISION RI, RO, AMAT( N1, N2 ), FACV0T, FACVT0,
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
      LOGICAL OV0T, OVT0, OMNP
      INTEGER IALPHA, LA, MA, ICSA,
     1        IBETA, LB, MB, ICSB,
     2        IGAMMA, LG, MG, ICSG,
     3        INTPAR( 7 ), IND1, IND2, IND3, INDSHC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of N2 against NH and NR etc.
C
      IF ( N2.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPMVGT, bad array size'
        PRINT *,' N2 (second array dimension) = ',N2
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C Early exit for zero factor
C
      OV0T = .TRUE.
      OVT0 = .TRUE.
C
      IF ( ABS( FACV0T ).LT.TOL ) OV0T = .FALSE.
      IF ( ABS( FACVT0 ).LT.TOL ) OVT0 = .FALSE.
C
      IF ( (.NOT. OV0T) .AND. (.NOT. OVT0) ) RETURN
C
      OMNP = .FALSE.
      DO IBETA = 1, NH
         IF (     MHT( IBETA ).EQ.3      .AND.
     1            MHL( IBETA ).EQ.0      .AND.
     2            MHM( IBETA ).EQ.0      .AND.
     3            MHC( IBETA ).EQ.1            ) OMNP = .TRUE.
      ENDDO
C
      DO IBETA = 1, NH
         IF (     MHT( IBETA ).EQ.3      .AND.
     1            MHL( IBETA ).EQ.0      .AND.
     2            MHM( IBETA ).EQ.0      .AND.
     3            MHC( IBETA ).EQ.1            ) THEN
C
C ...... Before we do any of the more complex interactions, let us
C        first look at the l=0 temperature harmonic. For this,
C        the non-linear term is simply
C
C        La ( La + 1 ) Pa( r ) / r * d/dr ( Theta_0^0 ( r ) )
C
         DO IALPHA = 1, NH
           IF ( MHT( IALPHA ).EQ.1 ) THEN
            LA = MHL( IALPHA )

            IF ( .NOT. OVT0 ) GOTO 421
c                 .
C                 ...... Fill in the v_new.grad(Theta_old) term ...
c                 .
             INTPAR( 1 ) = 1
             INTPAR( 2 ) = LA
             INTPAR( 3 ) = 0
             INTPAR( 4 ) = NR
             INTPAR( 5 ) = NH
             INTPAR( 6 ) = IBETA
             INTPAR( 7 ) = N2
C
            FAC = FACVT0
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IALPHA,
     1                   IALPHA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
 421        CONTINUE

            IF ( .NOT. OV0T ) GOTO 422
c                 .
C                 ...... Fill in the v_old.grad(Theta_new) term ...
c                 .

             IF (.NOT. OMNP) GOTO 422
             INTPAR( 1 ) = 2
             INTPAR( 2 ) = LA
             INTPAR( 3 ) = 0
             INTPAR( 4 ) = NR
             INTPAR( 5 ) = NH
             INTPAR( 6 ) = IALPHA
             INTPAR( 7 ) = N2
C
            FAC = FACV0T
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IBETA,
     1                   IALPHA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
 422        CONTINUE

           ENDIF
         ENDDO
        ENDIF
      ENDDO
C     .     end IBETA = 1, NH (for ibeta = theta_0^0 ....

c        .
C        ............................................................
C        . First loop around the NPOLH Theta harmonics ( IBETA )    .
C        ............................................................

      DO IBETA = 1, NH
        IF ( MHT( IBETA ).EQ.3 .AND. MHL( IBETA ).NE.0 ) THEN
         LB   = MHL( IBETA )
         MB   = MHM( IBETA )
         ICSB = MHC( IBETA )
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
          DO IALPHA = 1, NH
            IF ( MHT( IALPHA ).EQ.1 ) THEN
             LA = MHL( IALPHA )
             MA = MHM( IALPHA )
             ICSA = MHC( IALPHA )
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
             DO IGAMMA = 1, NH
              IF ( MHT( IGAMMA ).EQ.3 ) THEN
               LG = MHL( IGAMMA )
               IF ( LG.GT.0 ) THEN
                 MG   = MHM( IGAMMA )
                 ICSG = MHC( IGAMMA )
                 IND3 = INDSHC( LG, MG, ICSG )
                 QAB = SHC( IND3 )
               ELSE
                 QAB = ZCOEF
               ENDIF
               IF ( ABS(QAB).LT.TOL ) GOTO 701
c              .
C              .... We have a non-zero Q_{gamma}^{alpha beta} term
c              .
               IF ( .NOT. OVT0 ) GOTO 702
c              .
C              ...... Fill in the v_new.grad(fg+Theta_old) term ...
c              .

          INTPAR( 1 ) = 3
          INTPAR( 2 ) = LA
          INTPAR( 3 ) = 0
          INTPAR( 4 ) = NR
          INTPAR( 5 ) = NH
          INTPAR( 6 ) = IBETA
          INTPAR( 7 ) = N2
C
          DPPARS( 1 ) = QAB
C
          FAC = FACVT0
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
C......................................................................
C
 702         CONTINUE
               IF ( .NOT. OV0T ) GOTO 701
c              .
C              ...... Fill in the v_old.grad(Theta_new) term .....
c              .

          INTPAR( 1 ) = 4
          INTPAR( 2 ) = LA
          INTPAR( 3 ) = 0
          INTPAR( 4 ) = NR
          INTPAR( 5 ) = NH
          INTPAR( 6 ) = IALPHA
          INTPAR( 7 ) = N2
C
          DPPARS( 1 ) = QAB
C
          FAC = FACV0T
C
          CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IBETA,
     1                 IGAMMA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                 DPPARS, OLDVEC, FAC )
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
         DO IBETA = 1, NH
           IF ( MHT( IBETA ).EQ.3 .AND. MHL( IBETA ).NE.0 ) THEN
            LB   = MHL( IBETA )
            MB   = MHM( IBETA )
            ICSB = MHC( IBETA )
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
            DO IALPHA = 1, NH
              IF ( MHT( IALPHA ).EQ.1 ) THEN
               LA   = MHL( IALPHA )
               MA   = MHM( IALPHA )
               ICSA = MHC( IALPHA )
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
                DO IGAMMA = 1, NH
                 IF ( MHT( IGAMMA ).EQ.3 ) THEN
                  LG = MHL( IGAMMA )
                  IF ( LG.GT.0 ) THEN
                    MG   = MHM( IGAMMA )
                    ICSG = MHC( IGAMMA )
                    IND3 = INDSHC( LG, MG, ICSG )
                    SAB = SHC( IND3 )
                  ELSE
                    SAB = ZCOEF
                  ENDIF
                  IF ( ABS(SAB).LT.TOL ) GOTO 711
c                 .
C                 .... We have a non-zero S_{gamma}^{alpha beta} term
c                 .
                  IF ( .NOT. OVT0 ) GOTO 712
c                 .
C                 ...... Fill in the v_new.grad(fg+Theta_old) term ...
c                 .
             INTPAR( 1 ) = 5
             INTPAR( 2 ) = LA
             INTPAR( 3 ) = LB
             INTPAR( 4 ) = NR
             INTPAR( 5 ) = NH
             INTPAR( 6 ) = IBETA
             INTPAR( 7 ) = N2
C
             DPPARS( 1 ) = SAB
C
             FAC = FACVT0
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
C......................................................................
C
 712            CONTINUE
                  IF ( .NOT. OV0T ) GOTO 711
c                 .
C                 ...... Fill in the v_old.grad(Theta_new) term .....
c                 .
             INTPAR( 1 ) = 6
             INTPAR( 2 ) = LA
             INTPAR( 3 ) = LB
             INTPAR( 4 ) = NR
             INTPAR( 5 ) = NH
             INTPAR( 6 ) = IALPHA
             INTPAR( 7 ) = N2
C
             DPPARS( 1 ) = SAB
C
             FAC = FACV0T
C
             CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IBETA,
     1                   IGAMMA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
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
            DO IALPHA = 1, NH
              IF ( MHT( IALPHA ).EQ.2 ) THEN
               LA   = MHL( IALPHA )
               MA   = MHM( IALPHA )
               ICSA = MHC( IALPHA )
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
                DO IGAMMA = 1, NH
                 IF ( MHT( IGAMMA ).EQ.3 ) THEN
                  LG = MHL( IGAMMA )
                  IF ( LG.GT.0 ) THEN
                    MG   = MHM( IGAMMA )
                    ICSG = MHC( IGAMMA )
                    IND3 = INDSHC( LG, MG, ICSG )
                    TAB = SHC( IND3 )
                  ELSE
                    TAB = ZCOEF
                  ENDIF
                  IF ( ABS(TAB).LT.TOL ) GOTO 721
c                 .
C                 .... We have a non-zero T_{gamma}^{alpha beta} term
c                 .
                  IF ( .NOT. OVT0 ) GOTO 722
c                 .
C                 ...... Fill in the v_new.grad(fg+Theta_old) term ...
c                 .
C .....................................................................
C
             INTPAR( 1 ) = 7
             INTPAR( 2 ) = LA
             INTPAR( 3 ) = LB
             INTPAR( 4 ) = NR
             INTPAR( 5 ) = NH
             INTPAR( 6 ) = IBETA
             INTPAR( 7 ) = N2
C
             DPPARS( 1 ) = TAB
C
             FAC = FACVT0
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IALPHA,
     1                   IGAMMA, INTPAR, VDGTTS, ORD, RI, RO, AMAT,
     2                   DPPARS, OLDVEC, FAC )
C
C......................................................................
C
 722            CONTINUE
                  IF ( .NOT. OV0T ) GOTO 721
c                 .
C                 ...... Fill in the v_old.grad(Theta_new) term .....
c                 .
C
             INTPAR( 1 ) = 8
             INTPAR( 2 ) = LA
             INTPAR( 3 ) = LB
             INTPAR( 4 ) = NR
             INTPAR( 5 ) = NH
             INTPAR( 6 ) = IALPHA
             INTPAR( 7 ) = N2
C
             DPPARS( 1 ) = TAB
C
             FAC = FACV0T
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IBETA,
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
