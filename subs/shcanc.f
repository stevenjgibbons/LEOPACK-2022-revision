C*********************************************************************
C subroutine Spherical Harmonic Coef. Array Normalised Copy **********
C            -         -        -     -     -          -    **********
C Steve Gibbons Thu Mar 16 17:51:14 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Let the function g( \theta, \phi ) be expressed as the sum of      C
C Schmidt normalised spherical harmonics:                            C
C                                                                    C
C  g = \sum_{l,m} [  c_{l,mc} P_l^m( cos theta ) cos (m phi)         C
C               + c_{l,ms} P_l^m( cos theta ) sin (m phi)            C
C                                                                    C
C with the coefficients (ordered by the function INDSHC) given in    C
C the array SHC.                                                     C
C                                                                    C
C If all coefficients are zero, then SHCANC returns SHCN as a zero   C
C array. Otherwise, SHCN returns the coefficients scaled such that   C
C the spherical surface integral of g^2 is equal to VALN.            C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SHC       : Dim ( LH*(LH+2) ). Input array coefficients.       C
C     SHCN      : Dim ( LH*(LH+2) ). Output array coefficients.      C
C     VALN      : Value for normalisation.                           C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHCANC( LH, SHC, SHCN, VALN )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH
      DOUBLE PRECISION SHC( * ), SHCN( * ), VALN
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L, M, ICS, I, NH
      DOUBLE PRECISION DLOW, ZERO, PI, RNORM, COEF, FAC
      PARAMETER ( DLOW = 1.0d-8, ZERO = 0.0d0,
     1            PI = 3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NH = LH*(LH+2)
C
      RNORM = ZERO
      DO I = 1, NH
        CALL LMFIND( I, L, M, ICS )
        COEF = SHC( I )
        FAC  = 2.0d0*DBLE(L) + 1.0d0
        RNORM = RNORM + COEF*COEF*4.0d0*PI/FAC
      ENDDO
C
      IF ( RNORM.LT.DLOW ) THEN
        FAC = 0.0d0
      ELSE
        FAC = VALN/SQRT( RNORM )
      ENDIF
C
      DO I = 1, NH
        SHCN( I ) = SHC( I )*FAC
      ENDDO
C
      RETURN
      END
C*********************************************************************
