C*********************************************************************
C subroutine Changing Epsilon 1 Periodic Kumar Roberts Flow **********
C            -        -       - -        -     -       -    **********
C Steve Gibbons Wed Mar 28 10:21:55 MET DST 2001                     C
C____________________________________________________________________C
C                                                                    C
C We wish to parametrise a Kumar Roberts velocity to have an         C
C alternative value for EPS1 which varies as a cosine                C
C such that when CONST*DTIME = n for an integer n, EPS1 = EPS1A, and C
C when CONST*DTIME = ( n + 0.5 ), for an integer n, EPS1 = EPS1B.    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     CONST     : The "speed" of the change in MPAR.                 C
C     DTIME     : The "real" time of integration.                    C
C     EPS1A     : Value of EPS1 for CONST*DTIME = n (for an integer) C
C     EPS1B     : Value of EPS1 for CONST*DTIME = ( n + 0.5 )        C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     E1PAR     : eps1 -                                             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CE1KRF( CONST, DTIME, E1PAR, EPS1A, EPS1B )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION CONST, DTIME, E1PAR, EPS1A, EPS1B
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      DOUBLE PRECISION DPI, TWOPI, FAC1, FAC2, FAC3
      PARAMETER ( DPI = 3.14159265358979312D0, TWOPI = 2.0d0*DPI )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Calculate E1PAR
C
      FAC1 = TWOPI*CONST*DTIME
      FAC2 = EPS1A - EPS1B
      FAC3 = EPS1A + EPS1B
      E1PAR = 0.5d0*FAC2*DCOS( FAC1 ) + 0.5d0*FAC3
C
      RETURN
      END
C*********************************************************************
