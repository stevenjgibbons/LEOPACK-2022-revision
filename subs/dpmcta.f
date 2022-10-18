C*********************************************************************
C subroutine Double Precision Matrix Coriolis Term Add ***************
C            -      -         -      -        -    -   ***************
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C Adds the term to the double precision matrix resulting from        C
C the curl ( k cross v ) term in the momentum equation.              C
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
C   MHT( I ) = 4 for a poloidal magnetic field vector                C
C   MHT( I ) = 5 for a poloidal magnetic field vector                C
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
C     FAC       : Multiplication factor for term                     C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2,NTHPTS)       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
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
C     VELVF     : Vector function. dimensions ( NPHPTS, NTHPTS, 3)   C
C     QST       : Obvious ... ( dimension  (  LH*( LH+2 ) , 3 )      C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPMCTA ( N1, N2, NR, NH, KL, KU, KLE, IMF, MHT,
     1                    MHL, MHM, MHC, RI, RO, AMAT, FAC,
     2                    LH, NTHPTS, NPHPTS, GAUX, GAUW, QST,
     3                    FTF1, FTF2, FTF3, VELVF, PA, DPA, ORD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, NR, NH, KL, KU, KLE, IMF,
     1        MHT( NH ), MHM( NH ), MHC( NH ), MHL( NH ),
     2        LH, NTHPTS, NPHPTS
      DOUBLE PRECISION RI, RO, AMAT( N1, N2 ), FAC,
     1                 GAUX( NTHPTS ), GAUW( NTHPTS ),
     2                 PA( (LH+1)*(LH+2)/2 , NTHPTS ),
     3                 DPA( (LH+1)*(LH+2)/2 , NTHPTS )
      DOUBLE PRECISION FTF1( 2*NPHPTS ), QST( LH*(LH+2), 3 ),
     1                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ),
     2                 VELVF( NPHPTS, NTHPTS, 3 )
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      EXTERNAL CORTS
      INTEGER IIH, IOH, LI, MI, ICSI, INTPAR( 3 ), IPT,
     1        INDSHC, IHMI, LO, MO, ICSO, IHMO
      DOUBLE PRECISION DPPARS( 1 ), VEC( 1 ), LOW, TOL,
     1        QAB, SAB, TAB
      PARAMETER ( LOW = 1.0d-6, TOL = 1.0d-8 )
C nb - vec is not referred to by LAPTS
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of N2 against NH and NR etc.
C
      IF ( N2.NE.NH*NR ) THEN
        PRINT *,' Subroutine DPMCTA, bad array size'
        PRINT *,' N2 (second array dimension) = ',N2
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FAC ).LT.TOL ) RETURN
C
C     ............................................................
C     .                                                          .
C     .  First go through the poloidal harmonics.                .
C     .                                                          .
C     ............................................................

c     .
C     ...... LET IIH be the seed poloidal harmonic
      DO IIH = 1, NH
         IPT  = MHT( IIH )
         IF ( IPT.NE.1 ) GOTO 500
         LI   = MHL( IIH )
         MI   = MHM( IIH )
         ICSI = MHC( IIH )
         IHMI = INDSHC ( LI, MI, ICSI )
c        .
C        .........................................................
C        . NOW CALL CORCOF with ITYPE = 1, and IHARM = IHMA      .
C        . This will give us Q^q_(alpha beta),                   .
C        .                   S^q_(alpha beta),                   .
C        .                   T^q_(alpha beta).                   .
C        .........................................................
c        .

         CALL CORCOF ( 1, IHMI, QST, LH, LH, NPHPTS,
     1                 NTHPTS, NTHPTS, FTF1, FTF2, FTF3, VELVF,
     2                 GAUX, GAUW, PA, DPA )

c Now the Q^q_(alpha beta) are ( SHOULD BE!! ) all zero so ignore ..

c Loop around the NTORH harmonics beta

C .... these terms are the T_{ab}^q ....
C .... contribution to poloidal beta from poloidal alpha
C .
C******************************************************************
C Now the contrib. to poloidal harm 'b' from poloidal harm 'a' is
C
C    LEVa * ( LEVa + 1 ) * T_ab^q
C    ----------------------------  *  Pa( r )
C      SQRLL1( LEVb ) * RAD
C
C******************************************************************
C .
         DO IOH = 1, NH
           IPT  = MHT( IOH )
           IF ( IPT.NE.2 ) GOTO 517
           LO   = MHL( IOH )
           MO   = MHM( IOH )
           ICSO = MHC( IOH )
           IHMO = INDSHC ( LO, MO, ICSO )
           TAB  = QST( IHMO, 3 )
           IF ( ABS(TAB).LT.LOW ) GOTO 517
c
            DPPARS( 1 ) = TAB
            INTPAR( 1 ) = 1
            INTPAR( 2 ) = LI
            INTPAR( 3 ) = LO
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                   IOH, INTPAR, CORTS, ORD, RI, RO, AMAT,
     2                   DPPARS, VEC, FAC )
C
 517     CONTINUE
         ENDDO

C ..... end loop around the NTORH harmonics beta

c        .
C ...... We're now doing the S_{ab}^q terms ....
C ...... contribution to toroidal beta from poloidal alpha
c        .
C******************************************************************
C Now the contribution to toroidal harmonic 'b' from poloidal harm.
C 'a' is
C
C    LEVa * ( LEVa + 1 ) * S_ab^q     d Pa( r )
C    ----------------------------  *  ---------
C      SQRLL1( LEVb ) * RAD             d r
C
C******************************************************************
c        .

c Now loop around the NPOLH harmonics beta

         DO IOH = 1, NH
           IPT  = MHT( IOH )
           IF ( IPT.NE.1 ) GOTO 518
           LO   = MHL( IOH )
           MO   = MHM( IOH )
           ICSO = MHC( IOH )
           IHMO = INDSHC ( LO, MO, ICSO )
           SAB  = QST( IHMO, 2 )
           IF ( ABS(SAB).LT.LOW ) GOTO 518
c           .
            DPPARS( 1 ) = SAB
            INTPAR( 1 ) = 2
            INTPAR( 2 ) = LI
            INTPAR( 3 ) = LO
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                   IOH, INTPAR, CORTS, ORD, RI, RO, AMAT,
     2                   DPPARS, VEC, FAC )
C
 518     CONTINUE
         ENDDO
C           .
C           ... end of loop around beta harmonics

c           .
C           .........................................................
C           . NOW CALL CORCOF with ITYPE = 2, and IHARM = IHMA      .
C           . This will give us Q^s_(alpha beta),                   .
C           .                   S^s_(alpha beta),                   .
C           .                   T^s_(alpha beta).                   .
C           .........................................................
c           .

         CALL CORCOF ( 2, IHMI, QST, LH, LH, NPHPTS,
     1                 NTHPTS, NTHPTS, FTF1, FTF2, FTF3, VELVF,
     2                 GAUX, GAUW, PA, DPA )

c first deal with the Q^s_(alpha beta)

c        .
C        ... These terms are the Q_{ab}^s terms ....
C        ... contribution to toroidal beta from poloidal alpha
C        .

C        .
C******************************************************************
C Now the contribution to toroidal harmonic 'b' from poloidal harm.
C 'a' is
C
C     SQRLL1( LEVa )            [ Pa( r )      d Pa ( r )  ]
C  -  -------------- * Q_ab^s * [ -------   +  ----------  ]
C      ( RAD )**1               [    r           d r       ]
C
C******************************************************************
C           .

c Loop around the NPOLH Q harmonics beta
         DO IOH = 1, NH
           IPT  = MHT( IOH )
           IF ( IPT.NE.1 ) GOTO 519
           LO   = MHL( IOH )
           MO   = MHM( IOH )
           ICSO = MHC( IOH )
           IHMO = INDSHC ( LO, MO, ICSO )
           QAB = QST( IHMO, 1 )
           IF ( ABS(QAB).LT.LOW ) GOTO 519
c           .
            DPPARS( 1 ) = QAB
            INTPAR( 1 ) = 3
            INTPAR( 2 ) = LI
            INTPAR( 3 ) = LO
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                   IOH, INTPAR, CORTS, ORD, RI, RO, AMAT,
     2                   DPPARS, VEC, FAC )
C
 519     CONTINUE
         ENDDO
c .... end loop around NPOLH Q harmonics beta ....

c    .
C    .... these terms are the S_{ab}^s terms ...
C    .... contribution to toroidal beta from poloidal alpha
C    .

C******************************************************************
C Now the contribution to toroidal harmonic 'b' from poloidal harm.
C 'a' is
C
C   SQRLL1( LA )            [               2             ]
C   ------------ * S_ab^s * [  P_a''(r)  + ----  P_a'(r)  ]
C   SQRLL1( LB )            [               r             ]
C
C******************************************************************

c Loop around the NPOLH S harmonics beta
         DO IOH = 1, NH
           IPT  = MHT( IOH )
           IF ( IPT.NE.1 ) GOTO 520
           LO   = MHL( IOH )
           MO   = MHM( IOH )
           ICSO = MHC( IOH )
           IHMO = INDSHC ( LO, MO, ICSO )
           SAB = QST( IHMO , 2 )
           IF ( ABS(SAB).LT.LOW ) GOTO 520
C
            DPPARS( 1 ) = SAB
            INTPAR( 1 ) = 4
            INTPAR( 2 ) = LI
            INTPAR( 3 ) = LO
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                   IOH, INTPAR, CORTS, ORD, RI, RO, AMAT,
     2                   DPPARS, VEC, FAC )
C
 520     CONTINUE
         ENDDO
c .... end loop around NPOLH S harmonics beta ....

C .... these terms are the T_{ab}^s ....
C .... contribution to poloidal beta from poloidal alpha
C .
C******************************************************************
C
C Now the contribution to poloidal harmonic 'b' from poloidal harm
C 'a' is
C            SQRLL1( la )   [  P_a(r)             ]
C  T_ab^s  -------------- * [  ------  +  P_a'(r) ]
C            SQRLL1( lb )   [    r                ]
C
C
C******************************************************************
C .
C Loop around the NTORH T harmonics beta
         DO IOH = 1, NH
           IPT  = MHT( IOH )
           IF ( IPT.NE.2 ) GOTO 521
           LO   = MHL( IOH )
           MO   = MHM( IOH )
           ICSO = MHC( IOH )
           IHMO = INDSHC ( LO, MO, ICSO )
           TAB = QST( IHMO , 3 )
           IF ( ABS(TAB).LT.LOW ) GOTO 521
c
            DPPARS( 1 ) = TAB
            INTPAR( 1 ) = 5
            INTPAR( 2 ) = LI
            INTPAR( 3 ) = LO
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                   IOH, INTPAR, CORTS, ORD, RI, RO, AMAT,
     2                   DPPARS, VEC, FAC )
C
 521     CONTINUE
         ENDDO
c .... end loop around NTORH T harmonics beta ....


 500     CONTINUE
C line 500 is at the end of the IIH poloidal loop
      ENDDO
C     end of loop around alpha harmonics

C        ............................................................
C        .                                                          .
C        .  Now go through the toroidal harmonics.                  .
C        .                                                          .
C        ............................................................

C        .... LET IALPHA be the seed toroidal harmonic
      DO IIH = 1, NH
         IPT  = MHT( IIH )
         IF ( IPT.NE.2 ) GOTO 501
         LI   = MHL( IIH )
         MI   = MHM( IIH )
         ICSI = MHC( IIH )
         IHMI = INDSHC ( LI, MI, ICSI )

c        .
C        .........................................................
C        . NOW CALL CORCOF with ITYPE = 3, and IHARM = IHMA      .
C        . This will give us Q^t_(alpha beta),                   .
C        .                   S^t_(alpha beta),                   .
C        .                   T^t_(alpha beta).                   .
C        .........................................................
c        .

         CALL CORCOF ( 3, IHMI, QST, LH, LH, NPHPTS,
     1                 NTHPTS, NTHPTS, FTF1, FTF2, FTF3, VELVF,
     2                 GAUX, GAUW, PA, DPA )

c   .
C   .... These terms are the Q_{ab}^t
C   .... contribution from toroidal alpha to toroidal beta
c   .
C******************************************************************
C
C  Now the contribution to the toroidal harm 'b' from
C  toroidal harm 'a' is
C
C  SQRLL1( La ) * Q_ab^t * tor_a(r) / rad
C
C******************************************************************

c
C ....      Loop around the NPOLH Q harmonics
         DO IOH = 1, NH
           IPT  = MHT( IOH )
           IF ( IPT.NE.1 ) GOTO 529
           LO   = MHL( IOH )
           MO   = MHM( IOH )
           ICSO = MHC( IOH )
           IHMO = INDSHC ( LO, MO, ICSO )
           QAB = QST( IHMO, 1 )
           IF ( ABS(QAB).LT.LOW ) GOTO 529
c           .
            DPPARS( 1 ) = QAB
            INTPAR( 1 ) = 6
            INTPAR( 2 ) = LI
            INTPAR( 3 ) = LO
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                   IOH, INTPAR, CORTS, ORD, RI, RO, AMAT,
     2                   DPPARS, VEC, FAC )
C
 529     CONTINUE
         ENDDO
c .... end loop around NPOLH Q harmonics beta ....

c   .
C   ...... These are the S_{ab}^t terms
C   ...... contr. to toroidal beta from toroidal alpha
c   .

c
C ....   Loop around the NPOLH S harmonics
         DO IOH = 1, NH
           IPT  = MHT( IOH )
           IF ( IPT.NE.1 ) GOTO 531
           LO   = MHL( IOH )
           MO   = MHM( IOH )
           ICSO = MHC( IOH )
           IHMO = INDSHC ( LO, MO, ICSO )
           SAB = QST( IHMO, 2 )
           IF ( ABS(SAB).LT.LOW ) GOTO 531
c
            DPPARS( 1 ) = SAB
            INTPAR( 1 ) = 7
            INTPAR( 2 ) = LI
            INTPAR( 3 ) = LO
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                   IOH, INTPAR, CORTS, ORD, RI, RO, AMAT,
     2                   DPPARS, VEC, FAC )
C
 531     CONTINUE
         ENDDO

c  .
C  .... These are the T_{ab}^t terms ......
C  .... contr. to poloidal beta from toroidal alpha
C  .

c
C ....   Loop around the NTORH T harmonics
         DO IOH = 1, NH
           IPT  = MHT( IOH )
           IF ( IPT.NE.2 ) GOTO 532
           LO   = MHL( IOH )
           MO   = MHM( IOH )
           ICSO = MHC( IOH )
           IHMO = INDSHC ( LO, MO, ICSO )
           TAB = QST( IHMO, 3 )
           IF ( ABS(TAB).LT.LOW ) GOTO 532

            DPPARS( 1 ) = TAB
            INTPAR( 1 ) = 8
            INTPAR( 2 ) = LI
            INTPAR( 3 ) = LO
C
            CALL DPMES ( N1, N2, NR, NH, KL, KU, KLE, IMF, IIH,
     1                   IOH, INTPAR, CORTS, ORD, RI, RO, AMAT,
     2                   DPPARS, VEC, FAC )
C
 532     CONTINUE
         ENDDO


 501     CONTINUE
C line 501 is at the end of the IIH toroidal loop
       ENDDO
C      end of loop around alpha harmonics
c      .

C
      RETURN
      END
C*********************************************************************
