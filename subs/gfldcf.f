C*********************************************************************
C subroutine General First Left Difference Coefficients Fourth order C
C            -       -     -    -          -            -            C
C Steve Gibbons 28.08.97                                             C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     H         : Distance between radial grid nodes                 C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FIM1C     : Coeff.s for f( i - 1 ) Dimension ( 4 )             C
C     FIZOC     : Coeff.s for f(   i   ) Dimension ( 4 )             C
C     FIP1C     : Coeff.s for f( i + 1 ) Dimension ( 4 )             C
C     FIP2C     : Coeff.s for f( i + 2 ) Dimension ( 4 )             C
C     FIP3C     : Coeff.s for f( i + 3 ) Dimension ( 4 )             C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE GFLDCF ( H, FIM1C, FIZOC, FIP1C, FIP2C, FIP3C )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION H,  FIM1C( 4 ), FIZOC( 4 ), FIP1C( 4 ),
     1                 FIP2C( 4 ), FIP3C( 4 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION FAC1, FAC2, FAC3, FAC4
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check that H is Non-zero
      IF (H.LT.1d-10) THEN
         PRINT *,'Subroutine GFLDCF - Zero value for H received.'
         PRINT *,'Program aborted.'
         STOP
      ENDIF
C____________________________________________________________________C
      FAC1 = 1.0d0/(12.0d0*H)
      FAC2 = 1.0d0/(12.0d0*H*H)
      FAC3 = 1.0d0/(2.0d0*H*H*H)
      FAC4 = 1.0d0/(H*H*H*H)
C____________________________________________________________________C
C
      FIM1C( 1 ) = FAC1*(-3.0d0)
      FIM1C( 2 ) = FAC2*11.0d0
      FIM1C( 3 ) = FAC3*(-3.0d0)
      FIM1C( 4 ) = FAC4
C
      FIZOC( 1 ) = FAC1*(-10.0d0)
      FIZOC( 2 ) = FAC2*(-20.0d0)
      FIZOC( 3 ) = FAC3*10.0d0
      FIZOC( 4 ) = FAC4*(-4.0d0)
C
      FIP1C( 1 ) = FAC1*18.0d0
      FIP1C( 2 ) = FAC2*6.0d0
      FIP1C( 3 ) = FAC3*(-12.0d0)
      FIP1C( 4 ) = FAC4*6.0d0
C
      FIP2C( 1 ) = FAC1*(-6.0d0)
      FIP2C( 2 ) = FAC2*4.0d0
      FIP2C( 3 ) = FAC3*6.0d0
      FIP2C( 4 ) = FAC4*(-4.0d0)
C
      FIP3C( 1 ) = FAC1
      FIP3C( 2 ) = -FAC2
      FIP3C( 3 ) = -FAC3
      FIP3C( 4 ) = FAC4
C
      RETURN
      END
C*********************************************************************
