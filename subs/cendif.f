C*********************************************************************
C subroutine CENtral DIFference coefficients *************************
C            ---     ---                     *************************
C Steve Gibbons 18.11.98                                             C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     H         : Distance between radial grid nodes                 C
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
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     CM3       : Coeffs for f( i - 3 ) Dimension ( 4 )              C
C     CM2       : Coeffs for f( i - 2 ) Dimension ( 4 )              C
C     CM1       : Coeffs for f( i - 1 ) Dimension ( 4 )              C
C     C00       : Coeffs for f( i     ) Dimension ( 4 )              C
C     CP1       : Coeffs for f( i + 1 ) Dimension ( 4 )              C
C     CP2       : Coeffs for f( i + 2 ) Dimension ( 4 )              C
C     CP3       : Coeffs for f( i + 3 ) Dimension ( 4 )              C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CENDIF( H, ORD, CM3, CM2, CM1, C00, CP1, CP2, CP3 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION H,  CM3( 4 ), CM2( 4 ), CM1( 4 ),
     1                 C00( 4 ), CP1( 4 ), CP2( 4 ), CP3( 4 )
      CHARACTER *(2)   ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION FAC1, FAC2, FAC3, FAC4, ZERO
      INTEGER I
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check that H is Non-zero
      IF (H.LT.1d-10) THEN
         PRINT *,'Subroutine CENDIF - Zero value for H received.'
         PRINT *,'Program aborted.'
         STOP
      ENDIF
C____________________________________________________________________C
C
      IF ( ORD.NE.'SS' .AND. ORD.NE.'SF' .AND. ORD.NE.'O5' .AND.
     1     ORD.NE.'O7' ) THEN
         PRINT *,'Subroutine CENDIF - ORD = ',ORD
         PRINT *,'ORD must be either SS, SF, O5 or O7.'
         PRINT *,'Program aborted.'
         STOP
      ENDIF
C____________________________________________________________________C
C
      DO I = 1, 4
        CM3( I ) = ZERO
        CM2( I ) = ZERO
        CM1( I ) = ZERO
        C00( I ) = ZERO
        CP1( I ) = ZERO
        CP2( I ) = ZERO
        CP3( I ) = ZERO
      ENDDO
C____________________________________________________________________C
C
      IF ( ORD.EQ.'SS' ) THEN
c
         FAC1     = (1.0d0)/H
         FAC2     = FAC1/H
c
         CM1( 1 ) = (-0.5d0)*FAC1
         CP1( 1 ) = ( 0.5d0)*FAC1
c
         CM1( 2 ) = FAC2
         C00( 2 ) = (-2.0d0)*FAC2
         CP1( 2 ) = FAC2
c
      ENDIF
C____________________________________________________________________C
C
      IF ( ORD.EQ.'SS' .OR. ORD.EQ.'O5' ) THEN
c
         FAC3     = 1.0d0/(2.0d0*H*H*H)
         FAC4     = 1.0d0/(H*H*H*H)
c
         CM2( 3 ) = FAC3*(-1.0d0)
         CM1( 3 ) = FAC3*(2.0d0)
         CP1( 3 ) = FAC3*(-2.0d0)
         CP2( 3 ) = FAC3
c
         CM2( 4 ) = FAC4
         CM1( 4 ) = FAC4*(-4.0d0)
         C00( 4 ) = FAC4*( 6.0d0)
         CP1( 4 ) = FAC4*(-4.0d0)
         CP2( 4 ) = FAC4
c
      ENDIF
C____________________________________________________________________C
C
      IF ( ORD.EQ.'SF' .OR. ORD.EQ.'O5' ) THEN
c
         FAC1     = 1.0d0/(12.0d0*H)
         FAC2     = 1.0d0/(12.0d0*H*H)
c
         CM2( 1 ) = FAC1
         CM1( 1 ) = FAC1*(-8.0d0)
         CP1( 1 ) = FAC1*( 8.0d0)
         CP2( 1 ) = FAC1*(-1.0d0)
c
         CM2( 2 ) = FAC2*(-1.0d0)
         CM1( 2 ) = FAC2*(16.0d0)
         C00( 2 ) = FAC2*(-30.0d0)
         CP1( 2 ) = FAC2*(16.0d0)
         CP2( 2 ) = FAC2*(-1.0d0)
c
      ENDIF
C____________________________________________________________________C
C
      IF ( ORD.EQ.'SF' .OR. ORD.EQ.'O7' ) THEN
c
         FAC3     = 1.0d0/(8.0d0*H*H*H)
         FAC4     = 1.0d0/(6.0d0*H*H*H*H)
c
         CM3( 3 ) = FAC3
         CM2( 3 ) = FAC3*( -8.0d0)
         CM1( 3 ) = FAC3*( 13.0d0)
         CP1( 3 ) = FAC3*(-13.0d0)
         CP2( 3 ) = FAC3*(  8.0d0)
         CP3( 3 ) = FAC3*( -1.0d0)
c
         CM3( 4 ) = FAC4*( -1.0d0)
         CM2( 4 ) = FAC4*( 12.0d0)
         CM1( 4 ) = FAC4*(-39.0d0)
         C00( 4 ) = FAC4*( 56.0d0)
         CP1( 4 ) = FAC4*(-39.0d0)
         CP2( 4 ) = FAC4*( 12.0d0)
         CP3( 4 ) = FAC4*( -1.0d0)
c
      ENDIF
C____________________________________________________________________C
C
      IF ( ORD.EQ.'O7' ) THEN
c
         FAC1     = 1.0d0/(60.0d0*H)
         FAC2     = 1.0d0/(180.0d0*H*H)
c
         CM3( 1 ) = FAC1*( -1.0d0)
         CM2( 1 ) = FAC1*(  9.0d0)
         CM1( 1 ) = FAC1*(-45.0d0)
         CP1( 1 ) = FAC1*( 45.0d0)
         CP2( 1 ) = FAC1*( -9.0d0)
         CP3( 1 ) = FAC1
c
         CM3( 2 ) = FAC2*(  2.0d0)
         CM2( 2 ) = FAC2*(-27.0d0)
         CM1( 2 ) = FAC2*(270.0d0)
         C00( 2 ) = FAC2*(-490.0d0)
         CP1( 2 ) = FAC2*(270.0d0)
         CP2( 2 ) = FAC2*(-27.0d0)
         CP3( 2 ) = FAC2*(  2.0d0)
c
      ENDIF
C
      RETURN
      END
C*********************************************************************
