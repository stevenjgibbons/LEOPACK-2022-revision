C*********************************************************************
C subroutine Constant D Periodic Kumar Roberts Flow ******************
C            -        - -        -     -       -    ******************
C Steve Gibbons Wed Feb 28 16:39:19 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C We wish to parametrise a Kumar Roberts velocity to have a          C
C constant DPAR value and an MPAR value which varies as a cosine     C
C such that when CONST*DTIME = n for an integer n, MPAR = DM1, and   C
C when CONST*DTIME = ( n + 0.5 ), for an integer n, MPAR = DM2.      C
C                                                                    C
C DM1 and DM2 must be distinct. For a constant velocity, set CONST   C
C to zero - this will keep the flow perpetually at MPAR = DM1.       C
C                                                                    C
C RMPAR - the magnetic Reynolds number is RM1 at DM1 and is RM2 at   C
C DM2. RM1 and RM2 may be identical.                                 C
C                                                                    C
C CDPKRF outputs E0PAR, E1PAR, E2PAR, E3PAR and PPAR.                C
C                                                                    C
C PPAR is ALWAYS set to 3.0d0!                                       C
C                                                                    C
C E0PAR, E1PAR, E2PAR and E3PAR are scaled by RMPAR.                 C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DPAR      : The constant value of DPAR.                        C
C     CONST     : The "speed" of the change in MPAR.                 C
C     DTIME     : The "real" time of integration.                    C
C     DM1       : Value of MPAR for CONST*DTIME = n (for an integer) C
C     DM2       : Value of MPAR for CONST*DTIME = ( n + 0.5 )        C
C     RM1       : Value of RMPAR corresponding to MPAR = M1          C
C     RM2       : Value of RMPAR corresponding to MPAR = M2          C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RMPAR     : Magnetic reynolds number (should NOT be used to    C
C                  scale velocity parameters)                        C
C                                                                    C
C     PPAR      : Always set to 3.0d0                                C
C                                                                    C
C     E0PAR     : eps0 - scaled by magnetic reynolds number          C
C     E1PAR     : eps1 - scaled by magnetic reynolds number          C
C     E2PAR     : eps2 - scaled by magnetic reynolds number          C
C     E3PAR     : eps3 - scaled by magnetic reynolds number          C
C                                                                    C
C     MPAR      : Value of "M" :                                     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CDPKRF( DPAR, CONST, DTIME, DM1, DM2, RM1, RM2,
     1                RMPAR, PPAR, E0PAR, E1PAR, E2PAR, E3PAR, MPAR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION DPAR, CONST, DTIME, DM1, DM2, RM1, RM2,
     1                RMPAR, PPAR, E0PAR, E1PAR, E2PAR, E3PAR, MPAR
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      DOUBLE PRECISION DM2MM1, DTOL, DPI, TWOPI, FAC1,
     1                 DM1MM2, DM1PM2
      PARAMETER ( DTOL = 1.0d-9, DPI = 3.14159265358979312D0,
     1            TWOPI = 2.0d0*DPI )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First check that m1 and m2 are distinct
C
      IF ( DABS( DM1 - DM2 ).LT.DTOL ) THEN
        PRINT *,' Subroutine CDPKRF.'
        PRINT *,' DM1 = ', DM1
        PRINT *,' DM2 = ', DM2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Calculate MPAR
C
      DM2MM1 = DM2 - DM1
      DM1MM2 = DM1 - DM2
      DM1PM2 = DM1 + DM2
      FAC1 = TWOPI*CONST*DTIME
      MPAR = 0.5d0*DM1MM2*DCOS( FAC1 ) + 0.5d0*DM1PM2
C
C Calculate RMPAR
C
      RMPAR = (RM2-RM1)*MPAR/DM2MM1 + (RM1*DM2-RM2*DM1)/DM2MM1
C
C Calculate E0PAR, E1PAR, E2PAR, E3PAR
C
      CALL DM2E( DPAR, MPAR, E0PAR, E1PAR, E2PAR, E3PAR )
C
C Scale coefficients by Magnetic Reynolds number
C
      E0PAR = E0PAR * RMPAR
      E1PAR = E1PAR * RMPAR
      E2PAR = E2PAR * RMPAR
      E3PAR = E3PAR * RMPAR
C
      PPAR = 3.0d0
C
      RETURN
      END
C*********************************************************************
