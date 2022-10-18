C*********************************************************************
C subroutine Periodic Function 2 Sin Only Function *******************
C            -        -        - -   -    -        *******************
C Steve Gibbons Mon Mar 20 11:43:18 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Let F^{M} be a periodic function of the form                       C
C                                                                    C
C  F^{M}( phi ) = FMC cos ( m phi )  +  FMS sin ( m phi )            C
C                                                                    C
C PF2SOF returns the values GMS and TAU such that                    C
C                                                                    C
C  F^{M}( phi ) = GMS sin ( m phi + m tau )                          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     M        : Wavenumber, m.                                      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     FMC      : See above equation. (Input)                         C
C     FMS      : See above equation. (Input)                         C
C     GMS      : See above equation. (Output)                        C
C     TAU      : See above equation. (Output)                        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PF2SOF( M, FMC, FMS, GMS, TAU )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION FMC, FMS, GMS, TAU
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION PI, DM, DLOW, DMTAU, DFRAC
      PARAMETER ( PI=3.14159265358979312D0, DLOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      IF ( M.EQ.0 ) THEN
        PRINT *,' Subroutine PF2SOF. M = 0 !!!'
        PRINT *,' You will not find a sin(m phi)'
        PRINT *,' function. Program aborted.'
        STOP
      ENDIF
C     .
C     . Case where both coefficients are zero
C     .
      IF (        DABS( FMC ).LT.DLOW        .AND.
     1            DABS( FMS ).LT.DLOW        )    THEN
        GMS = 0.0d0
        TAU = 0.0d0
        RETURN
      ENDIF
C     .
C     . Case where we have a purely cosine function
C     .
      IF (        DABS( FMS ).LT.DLOW        )    THEN
        GMS = FMC
        TAU = 0.5d0*PI
        RETURN
      ENDIF
C     .
C     . We can now do the general case without fear
C     . of dividing by zero ...
C     .
      DM    = DBLE( M )
      DFRAC = FMC/FMS
      DMTAU = DATAN( DFRAC )
      TAU   = DMTAU/DM
      GMS   = FMS/DCOS( DMTAU )
C     .
      RETURN
      END
C*********************************************************************
