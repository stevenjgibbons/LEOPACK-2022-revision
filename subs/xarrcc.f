C*********************************************************************
C subroutine XARR Compatibility Check ********************************
C            ---- -             -     ********************************
C Steve Gibbons Mon Apr  3 09:30:19 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C For programs with conducting inner cores, we need two arrays of    C
C r values for the two solution vectors. This routine simply checks  C
C that the last NR nodes of the XARRM arrays correspond to the NR    C
C nodes of the XARR array.                                           C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes (outer core only).            C
C     NRMF      : Number of grid nodes (magnetic field).             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     XARRM     : Array of dimension (  NRMF  ).                     C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XARRCC( NR, XARR, NRMF, XARRM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NRMF
      DOUBLE PRECISION XARR( NR ), XARRM( NRMF )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR, NRIC, IRMF
      DOUBLE PRECISION TOL, X1, X2
      PARAMETER ( TOL = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NRIC = NRMF - NR
      IF ( NRIC.LT.0 ) THEN
        PRINT *,' Subroutine XARRCC.'
        PRINT *,' NR = ', NR,' NRMF = ', NRMF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO IR = 1, NR
        IRMF = IR + NRIC
        X1   = XARR( IR )
        X2   = XARRM( IRMF )
        IF ( DABS( X1 - X2 ).GT.TOL ) THEN
          PRINT *,' Subroutine XARRCC.'
          PRINT *,' There are ',NRIC,' inner core nodes.'
          PRINT *,' xarr( ',ir,') = ',X1
          PRINT *,' xarrm(',irmf,') = ',X2
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDDO
C
      RETURN
      END
C*********************************************************************
