C*********************************************************************
C subroutine Solid Body Rotation Coefficient Testing Subroutine ******
C            -     -    -        -           -       -          ******
C Steve Gibbons Tue Dec 14 11:12:08 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Takes an input of L, M, ICS and returns all the coefficients       C
C up to degree and order LH as given in equations (B.42) through to  C
C (B.50) in my thesis.                                               C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     L         : Spherical harmonic degree, l.                      C
C     M         : Spherical harmonic order, m.                       C
C     ICS       : 1 for cos m phi, 2 for sin m phi ...               C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NPHPTS    : The number of phi points.                          C
C     NTHPTS    : The number of theta points.                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Working Arrays   :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF1      : Array for Fourier transforming dim. (2*NPHPTS)     C
C     FTF2      : Array for Fourier transforming dim. (2*NPHPTS)     C
C     FTF3      : Array for Fourier transforming dim. (2*NPHPTS)     C
C     VF1       : Vector function. dimensions ( NPHPTS, NTHPTS, 3)   C
C     VF2       : Vector function. dimensions ( NPHPTS, NTHPTS, 3)   C
C     VF3       : Vector function. dimensions ( NPHPTS, NTHPTS, 3)   C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SBRCTS( L, M, ICS, LH, NTHPTS, NPHPTS, GAUX, GAUW,
     1       PA, DPA, FTF1, FTF2, FTF3, QST, VF1, VF2, VF3, KQSAG,
     2       KQTAG, KSQAG, KSSAG, KSTAG, KTQAG, KTSAG, KTTAG,
     3       QSSARG, QSTARG, SQSARG, SQTARG, SSQARG, TQSARG, TQTARG,
     4       TSQARG, TQSRAG, TQTRAG, TSQRAG, TTQRAG )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER L, M, ICS, LH, NPHPTS, NTHPTS
      DOUBLE PRECISION
     1                 KQSAG( LH*(LH+2) ), KQTAG( LH*(LH+2) ), 
     2                 KSQAG( LH*(LH+2) ), KSSAG( LH*(LH+2) ), 
     3                 KSTAG( LH*(LH+2) ), KTQAG( LH*(LH+2) ), 
     4                 KTSAG( LH*(LH+2) ), KTTAG( LH*(LH+2) )
      DOUBLE PRECISION
     1                 QSSARG( LH*(LH+2) ), QSTARG( LH*(LH+2) ),
     2                 SQSARG( LH*(LH+2) ), SQTARG( LH*(LH+2) ),
     3                 SSQARG( LH*(LH+2) ), TQSARG( LH*(LH+2) )
      DOUBLE PRECISION
     1                 TQTARG( LH*(LH+2) ), TSQARG( LH*(LH+2) ),
     2                 TQSRAG( LH*(LH+2) ), TQTRAG( LH*(LH+2) ),
     3                 TSQRAG( LH*(LH+2) ), TTQRAG( LH*(LH+2) )
      DOUBLE PRECISION
     1                 GAUX( NTHPTS ), GAUW( NTHPTS ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
      DOUBLE PRECISION
     1                 VF1( NPHPTS, NTHPTS, 3),
     2                 VF2( NPHPTS, NTHPTS, 3),
     3                 VF3( NPHPTS, NTHPTS, 3)
      DOUBLE PRECISION
     1                 FTF1( 2*NPHPTS ), FTF3( 2*NPHPTS ),
     2                 FTF2( 2*NPHPTS ), QST( LH*(LH+2), 3)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION ZERO, COEF, ZCOEF
      INTEGER IQST, IQSTS, LS, MS, ICSS, IOP, ILEN, IH
      PARAMETER ( LS = 1, MS = 0, ICSS = 1, IOP = 0, ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Zero all output arrays
C
      ILEN = LH*(LH+2)
      CALL VECOP ( KQSAG, ZERO, ILEN, IOP )
      CALL VECOP ( KQTAG, ZERO, ILEN, IOP )
      CALL VECOP ( KSQAG, ZERO, ILEN, IOP )
      CALL VECOP ( KSSAG, ZERO, ILEN, IOP )
      CALL VECOP ( KSTAG, ZERO, ILEN, IOP )
      CALL VECOP ( KTQAG, ZERO, ILEN, IOP )
      CALL VECOP ( KTSAG, ZERO, ILEN, IOP )
      CALL VECOP ( KTTAG, ZERO, ILEN, IOP )
C
      CALL VECOP ( QSSARG, ZERO, ILEN, IOP )
      CALL VECOP ( QSTARG, ZERO, ILEN, IOP )
      CALL VECOP ( SQSARG, ZERO, ILEN, IOP )
      CALL VECOP ( SQTARG, ZERO, ILEN, IOP )
      CALL VECOP ( SSQARG, ZERO, ILEN, IOP )
      CALL VECOP ( TQSARG, ZERO, ILEN, IOP )
      CALL VECOP ( TQTARG, ZERO, ILEN, IOP )
      CALL VECOP ( TSQARG, ZERO, ILEN, IOP )
      CALL VECOP ( TQSRAG, ZERO, ILEN, IOP )
      CALL VECOP ( TQTRAG, ZERO, ILEN, IOP )
      CALL VECOP ( TSQRAG, ZERO, ILEN, IOP )
      CALL VECOP ( TTQRAG, ZERO, ILEN, IOP )
C
C Firstly, the Coriolis coefficients are calculated very
C straightforwardly with a call to CORCOF.
C That is its exact purpose.
C
      DO IQST = 1, 3
        CALL CORCOF( IQST, L, M, ICS, QST, LH, NPHPTS, NTHPTS,
     1             FTF1, FTF2, FTF3, VF1, GAUX, GAUW, PA, DPA, LH )
        DO IQSTS = 1, 3
         DO IH = 1, ILEN
          COEF = QST( IH, IQSTS )
          IF ( IQST.EQ.1 .AND. IQSTS.EQ.2 ) KQSAG( IH ) = COEF
          IF ( IQST.EQ.1 .AND. IQSTS.EQ.3 ) KQTAG( IH ) = COEF
          IF ( IQST.EQ.2 .AND. IQSTS.EQ.1 ) KSQAG( IH ) = COEF
          IF ( IQST.EQ.2 .AND. IQSTS.EQ.2 ) KSSAG( IH ) = COEF
          IF ( IQST.EQ.2 .AND. IQSTS.EQ.3 ) KSTAG( IH ) = COEF
          IF ( IQST.EQ.3 .AND. IQSTS.EQ.1 ) KTQAG( IH ) = COEF
          IF ( IQST.EQ.3 .AND. IQSTS.EQ.2 ) KTSAG( IH ) = COEF
          IF ( IQST.EQ.3 .AND. IQSTS.EQ.3 ) KTTAG( IH ) = COEF
         ENDDO
        ENDDO
      ENDDO
C
C --**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1
C Evaluate q (user specified) in VF1
C --**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1
C
      IQST = 1
      CALL SHVECT ( L, M, ICS, IQST, VF1, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Evaluate s_1^0 in VF2
C
      IQSTS = 2
      CALL SHVECT ( LS, MS, ICSS, IQSTS, VF2, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Now calculate VF3 = VF1 x VF2
C
      CALL VFCP ( VF1, VF2, VF3, NPHPTS, NTHPTS )
C
C Transform into QST coefficients
C
      CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1             FTF3, ZCOEF, LH, NTHPTS, NPHPTS, LH )
C
C Loop around the harmonics filling QSSARG and QSTARG
C
      DO IH = 1, ILEN
        QSSARG( IH ) = QST( IH, 2 )
        QSTARG( IH ) = QST( IH, 3 )
      ENDDO
C
C --**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1
C Evaluate s (user specified) in VF1
C --**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1
C
      IQST = 2
      CALL SHVECT ( L, M, ICS, IQST, VF1, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Evaluate q_1^0 in VF2
C
      IQSTS = 1
      CALL SHVECT ( LS, MS, ICSS, IQSTS, VF2, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Now calculate VF3 = VF1 x VF2
C
      CALL VFCP ( VF1, VF2, VF3, NPHPTS, NTHPTS )
C
C Transform into QST coefficients
C
      CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1             FTF3, ZCOEF, LH, NTHPTS, NPHPTS, LH )
C
C Loop around the harmonics filling SQSARG and SQTARG
C
      DO IH = 1, ILEN
        SQSARG( IH ) = QST( IH, 2 )
        SQTARG( IH ) = QST( IH, 3 )
      ENDDO
C
C Evaluate s_1^0 in VF2
C
      IQSTS = 2
      CALL SHVECT ( LS, MS, ICSS, IQSTS, VF2, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Now calculate VF3 = VF1 x VF2
C
      CALL VFCP ( VF1, VF2, VF3, NPHPTS, NTHPTS )
C
C Transform into QST coefficients
C
      CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1             FTF3, ZCOEF, LH, NTHPTS, NPHPTS, LH )
C
C Loop around the harmonics filling SSQARG
C
      DO IH = 1, ILEN
        SSQARG( IH ) = QST( IH, 1 )
      ENDDO
C
C --**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1
C Evaluate t (user specified) in VF1
C --**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1
C
      IQST = 3
      CALL SHVECT ( L, M, ICS, IQST, VF1, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Evaluate q_1^0 in VF2
C
      IQSTS = 1
      CALL SHVECT ( LS, MS, ICSS, IQSTS, VF2, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Now calculate VF3 = VF1 x VF2
C
      CALL VFCP ( VF1, VF2, VF3, NPHPTS, NTHPTS )
C
C Transform into QST coefficients
C
      CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1             FTF3, ZCOEF, LH, NTHPTS, NPHPTS, LH )
C
C Loop around the harmonics filling TQSARG and TQTARG
C
      DO IH = 1, ILEN
        TQSARG( IH ) = QST( IH, 2 )
        TQTARG( IH ) = QST( IH, 3 )
      ENDDO
C
C Evaluate s_1^0 in VF2
C
      IQSTS = 2
      CALL SHVECT ( LS, MS, ICSS, IQSTS, VF2, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Now calculate VF3 = VF1 x VF2
C
      CALL VFCP ( VF1, VF2, VF3, NPHPTS, NTHPTS )
C
C Transform into QST coefficients
C
      CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1             FTF3, ZCOEF, LH, NTHPTS, NPHPTS, LH )
C
C Loop around the harmonics filling TSQARG
C
      DO IH = 1, ILEN
        TSQARG( IH ) = QST( IH, 1 )
      ENDDO
C
C --**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1
C Evaluate t_1^0 in VF1
C --**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1**--**VF1
C
      IQSTS = 3
      CALL SHVECT ( LS, MS, ICSS, IQSTS, VF1, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Evaluate q (user specified) in VF2
C
      IQST = 1
      CALL SHVECT ( L, M, ICS, IQST, VF2, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Now calculate VF3 = VF1 x VF2
C
      CALL VFCP ( VF1, VF2, VF3, NPHPTS, NTHPTS )
C
C Transform into QST coefficients
C
      CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1             FTF3, ZCOEF, LH, NTHPTS, NPHPTS, LH )
C
C Loop around the harmonics filling TQSRAG and TQTRAG
C
      DO IH = 1, ILEN
        TQSRAG( IH ) = QST( IH, 2 )
        TQTRAG( IH ) = QST( IH, 3 )
      ENDDO
C
C Evaluate s_1^0 in VF2
C
      IQST = 2
      CALL SHVECT ( L, M, ICS, IQST, VF2, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Now calculate VF3 = VF1 x VF2
C
      CALL VFCP ( VF1, VF2, VF3, NPHPTS, NTHPTS )
C
C Transform into QST coefficients
C
      CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1             FTF3, ZCOEF, LH, NTHPTS, NPHPTS, LH )
C
C Loop around the harmonics filling TSQRAG
C
      DO IH = 1, ILEN
        TSQRAG( IH ) = QST( IH, 1 )
      ENDDO
C
C Evaluate t_1^0 in VF2
C
      IQST = 3
      CALL SHVECT ( L, M, ICS, IQST, VF2, GAUX, PA, DPA,
     1              NTHPTS, NPHPTS, LH )
C
C Now calculate VF3 = VF1 x VF2
C
      CALL VFCP ( VF1, VF2, VF3, NPHPTS, NTHPTS )
C
C Transform into QST coefficients
C
      CALL VF2QST( QST, VF3, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1             FTF3, ZCOEF, LH, NTHPTS, NPHPTS, LH )
C
C Loop around the harmonics filling TTQRAG
C
      DO IH = 1, ILEN
        TTQRAG( IH ) = QST( IH, 1 )
      ENDDO
C
      RETURN
      END
C*********************************************************************

