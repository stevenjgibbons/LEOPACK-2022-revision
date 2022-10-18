C*********************************************************************
C subroutine Constant D Cosine Bell Kumar Roberts flow ***************
C            -        - -      -    -     -            ***************
C Steve Gibbons Wed Nov  7 09:46:50 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C We wish to parametrise a Kumar Roberts velocity to have a          C
C constant DPAR value and an MPAR value which varies as follows:     C
C MPAR will be constant and equal to DM1 until the time              C
C                                                                    C
C    DTIME = DPEAKT - 0.5d0/CONST                                    C
C                                                                    C
C when MPAR will begin to vary in a cosine bell, reaching the        C
C value DM2 at DTIME = DPEAKT.                                       C
C By the time DTIME = DPEAKT + 0.5d0/CONST, MPAR will have reached   C
C the value DM1 again where it will stay.                            C
C                                                                    C
C DM1 and DM2 must be distinct. For a constant velocity, set CONST   C
C to zero - this will keep the flow perpetually at MPAR = DM1.       C
C                                                                    C
C RMPAR - the magnetic Reynolds number is RM1 at DM1 and is RM2 at   C
C DM2. RM1 and RM2 may be identical.                                 C
C                                                                    C
C CDCBKR outputs E0PAR, E1PAR, E2PAR, E3PAR and PPAR.                C
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
C     DPEAKT    : Time at which value M = DM2 must be reached.       C
C     DPAR      : The constant value of DPAR.                        C
C     CONST     : The "speed" of the change in MPAR.                 C
C     DTIME     : The "real" time of integration.                    C
C     DM1       : Value of MPAR for CONST*DTIME = n (for an integer) C
C                   before time DTIME = DPEAKT - 0.5d0/CONST         C
C                   and after time DTIME = DPEAKT + 0.5d0/CONST.     C
C                                                                    C
C     DM2       : Value of MPAR at time DTIME = DPEAKT               C
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
      SUBROUTINE CDCBKR( DPEAKT, DPAR, CONST, DTIME, DM1, DM2, RM1,
     1         RM2, RMPAR, PPAR, E0PAR, E1PAR, E2PAR, E3PAR, MPAR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION DPEAKT, DPAR, CONST, DTIME, DM1, DM2, RM1,
     1           RM2, RMPAR, PPAR, E0PAR, E1PAR, E2PAR, E3PAR, MPAR
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      DOUBLE PRECISION DM2MM1, DTOL, DPI, TWOPI, FAC1,
     1                 DM1MM2, DM1PM2, DTIMED, DELTAT
      PARAMETER ( DTOL = 1.0d-9, DPI = 3.14159265358979312D0,
     1            TWOPI = 2.0d0*DPI )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First check that m1 and m2 are distinct
C and that const is non-zero ...
C
      IF ( DABS( DM1 - DM2 ).LT.DTOL    .OR.
     1     CONST.LT.DTOL      ) THEN
        PRINT *,' Subroutine CDCBKR.'
        PRINT *,' DM1    = ', DM1
        PRINT *,' DM2    = ', DM2
        PRINT *,' CONST  = ', CONST
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Calculate DTIMED, DELTAT
C
      DTIMED = DTIME - DPEAKT
      DELTAT = 0.5d0/CONST
C
C Calculate MPAR
C
      DM2MM1 = DM2 - DM1
      DM1MM2 = DM1 - DM2
      DM1PM2 = DM1 + DM2
      IF ( DABS( DTIMED ).LT.DELTAT ) THEN
        FAC1 = TWOPI*CONST*(DTIMED - DELTAT)
      ELSE
        FAC1 = 0.0d0
      ENDIF
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
