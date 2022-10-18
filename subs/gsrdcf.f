C*********************************************************************
C subroutine General Second Right Difference Coefficients Fourth order
C            -       -      -     -          -            -          C
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
C     FIM4C     : Coeff.s for f( i - 4 ) Dimension ( 4 )             C
C     FIM3C     : Coeff.s for f( i - 3 ) Dimension ( 4 )             C
C     FIM2C     : Coeff.s for f( i - 2 ) Dimension ( 4 )             C
C     FIM1C     : Coeff.s for f( i - 1 ) Dimension ( 4 )             C
C     FIZOC     : Coeff.s for f(   i   ) Dimension ( 4 )             C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE GSRDCF ( H, FIM4C, FIM3C, FIM2C, FIM1C, FIZOC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION H, FIM4C( 4 ) , FIM3C( 4 ) , FIM2C( 4 ) ,
     1                 FIM1C( 4 ) , FIZOC( 4 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION FAC1, FAC2, FAC3, FAC4
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check that H is Non-zero
      IF (H.LT.1d-10) THEN
         PRINT *,'Subroutine GSRDCF - Zero value for H received.'
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
      FIM4C( 1 ) = FAC1*3.0d0
      FIM4C( 2 ) = FAC2*11.0d0
      FIM4C( 3 ) = FAC3*3.0d0
      FIM4C( 4 ) = FAC4
C
      FIM3C( 1 ) = FAC1*(-16.0d0)
      FIM3C( 2 ) = FAC2*(-56.0d0)
      FIM3C( 3 ) = FAC3*(-14.0d0)
      FIM3C( 4 ) = FAC4*(-4.0d0)
C
      FIM2C( 1 ) = FAC1*36.0d0
      FIM2C( 2 ) = FAC2*114.0d0
      FIM2C( 3 ) = FAC3*24.0d0
      FIM2C( 4 ) = FAC4*6.0d0
C
      FIM1C( 1 ) = FAC1*(-48.0d0)
      FIM1C( 2 ) = FAC2*(-104.0d0)
      FIM1C( 3 ) = FAC3*(-18.0d0)
      FIM1C( 4 ) = FAC4*(-4.0d0)
C
      FIZOC( 1 ) = FAC1*25.0d0
      FIZOC( 2 ) = FAC2*35.0d0
      FIZOC( 3 ) = FAC3*5.0d0
      FIZOC( 4 ) = FAC4
C
      RETURN
      END
C*********************************************************************

