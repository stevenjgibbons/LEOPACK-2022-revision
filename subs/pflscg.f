C*********************************************************************
C subroutine Periodic Function Longitude Shift Coefficient Give ******
C            -        -        -         -     -           -    ******
C Steve Gibbons Mon Mar 20 13:29:49 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Let F^{M} be a periodic function of the form                       C
C                                                                    C
C  F^{M}( phi ) = FMC cos ( m phi )  +  FMS sin ( m phi )            C
C                                                                    C
C For a given longitude, tau, PFLSCG returns GMC and GMS such that   C
C                                                                    C
C  F^{M}( phi ) =     GMS sin ( m phi + m tau )                      C
C                  +  GMC cos ( m phi + m tau )                      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     M        : Wavenumber, m.                                      C
C     TAU      : See above equation. (Input)                         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     FMC      : See above equation. (Input)                         C
C     FMS      : See above equation. (Input)                         C
C     GMC      : See above equation. (Output)                        C
C     GMS      : See above equation. (Output)                        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PFLSCG( M, TAU, FMC, FMS, GMC, GMS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION FMC, FMS, GMC, GMS, TAU
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IPIV( 2 ), ITWO, INFO, IONE
      DOUBLE PRECISION DM, DMTAU, DSMT, DCMT,
     1                 RHS( 2, 1 ), DMAT( 2, 2 )
      CHARACTER *(1) TRANS
      PARAMETER ( IONE = 1, ITWO = 2, TRANS = 'N' )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      DM    = DBLE( M )
      DMTAU = DM*TAU
      DSMT  = DSIN( DMTAU )
      DCMT  = DCOS( DMTAU )
C     .
C     . We need to solve a 2x2 matrix equation.
C     .
      RHS( 1, 1 ) = FMC
      RHS( 2, 1 ) = FMS
C     .
      DMAT( 1, 1 ) = DCMT
      DMAT( 1, 2 ) = DSMT
      DMAT( 2, 1 ) = (-1.0d0)*DSMT
      DMAT( 2, 2 ) = DCMT
C     .
C     . Perform LU decomposition on DMAT
C     .
      CALL DGETRF( ITWO, ITWO, DMAT, ITWO, IPIV, INFO )
C     .
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine PFLSCG.'
        PRINT *,' LAPACK subroutine DGETRF has returned'
        PRINT *,' the value INFO = ', INFO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Now solve the linear system
C     .
      CALL DGETRS( TRANS, ITWO, IONE, DMAT, ITWO, IPIV, RHS,
     1             ITWO, INFO )
C     .
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine PFLSCG.'
        PRINT *,' LAPACK subroutine DGETRS has returned'
        PRINT *,' the value INFO = ', INFO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      GMC = RHS( 1, 1 )
      GMS = RHS( 2, 1 )
C     .
      RETURN
      END
C*********************************************************************
