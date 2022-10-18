C*********************************************************************
C double precision function SECDER (second derivative) ***************
C Steve Gibbons 17.4.97                                              C
C____________________________________________________________________C
C Calculates d^2 P_l^m(theta) / d(theta)^2                           C
C by Legendre's equation ...                                         C
C                                                                    C
C  = [  { m*m-l*(l+1)*SIN^2(theta) }*P_l^m(theta) -                  C
C      COS(theta)*SIN(theta)* dP_l^m(theta)/d(theta). ]/sin^2(theta) C
C____________________________________________________________________C
C                                                                    C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     L         : Spherical harmonic degree, l.                      C
C     M         : Spherical harmonic order, m.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     X         : Cos(theta)                                         C
C     S         : Sin(theta)                                         C
C     POL       : P_l^m(theta)                                       C
C     DER1      : d P_l^m(theta)/ d theta                            C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION SECDER ( L, M, C, S, POL, DER1 )
      IMPLICIT NONE
      DOUBLE PRECISION SECDER
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER L, M
      DOUBLE PRECISION C, S, POL, DER1
C____________________________________________________________________C
C Variable declarations - Working Variables .........................C
      DOUBLE PRECISION RM, RL, DLOW
      PARAMETER ( DLOW = 1.0d-8 )
C____________________________________________________________________C
      RM = DBLE( M )
      RL = DBLE( L )
      IF ( DABS( S ).LT.DLOW ) THEN
        PRINT *,' Function SECDER'
        PRINT *,' S = ',S,' Div. by zero error imminent.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      SECDER = RM*RM - RL*( RL+1.0d0 )*S*S
      SECDER = SECDER*POL  - C*S*DER1
      SECDER = SECDER/(S*S)
      RETURN
      END
C*********************************************************************
