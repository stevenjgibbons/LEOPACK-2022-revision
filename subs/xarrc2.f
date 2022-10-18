C*********************************************************************
C subroutine XARR Compatibility check 2 ******************************
C            ---- -                   - ******************************
C Steve Gibbons Thu Dec  6 10:06:39 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C For programs with conducting inner cores, we need two arrays of    C
C r values for the two solution vectors. This routine checks that    C
C NR nodes of the XARRM array corresponds to the NR nodes of the     C
C XARRV array. The arrays can ofcourse be equal, which means that    C
C there is an insulating core and mantle. If NRM is greater than NRV C
C then there is a conducting region outside of the outer core.       C
C This may be a conducting inner core, a conducting layer at the     C
C outer boundary, or both.                                           C
C                                                                    C
C XARRC2 returns NRIC which is the number of grid nodes in the       C
C inner core (not including the node which is the IC/OC boundary.    C
C                                                                    C
C XARRM( NRIC + 1   ) = XARRV( 1 )                                   C
C XARRM( NRIC + NRV ) = XARRV( NRV )                                 C
C                                                                    C
C XARRM( i ) with i .gt. (NRIC + NRV) is in the mantle.              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NRV       : Number of grid nodes (outer core only).            C
C     NRM       : Number of grid nodes (magnetic field).             C
C     NRIC      : (Output) Number of grid nodes in inner core.       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARRV     : Array of dimension (  NRV  ).                      C
C                  XARRV( i ) contains the value of x or r at the    C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     XARRM     : Array of dimension (  NRM  ).                      C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XARRC2( NRV, XARRV, NRM, XARRM, NRIC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NRV, NRM, NRIC
      DOUBLE PRECISION XARRV( NRV ), XARRM( NRM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRV, IRM, NRCS
      DOUBLE PRECISION DTOL, X1, X2, RICB
      PARAMETER ( DTOL = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C NRCS is the number of nodes in a conducting solid
C
      NRCS = NRM - NRV
      IF ( NRCS.LT.0 ) THEN
        PRINT *,' Subroutine XARRC2.'
        PRINT *,' NRV = ', NRV,' NRM = ', NRM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RICB = XARRV( 1 )
C
      DO IRM = 1, NRM
        X1   = XARRM( IRM )
        IF ( DABS( X1 - RICB ).LT.DTOL ) THEN
          NRIC = IRM - 1
          GOTO 50
        ENDIF
      ENDDO
 50   CONTINUE
C
      DO IRV = 1, NRV
        IRM  = IRV + NRIC
        X1   = XARRV( IRV )
        X2   = XARRM( IRM )
        IF ( DABS( X1 - X2 ).GT.DTOL ) THEN
          PRINT *,' Subroutine XARRC2.'
          PRINT *,' There are ',NRIC,' inner core nodes.'
          PRINT *,' xarrv(',irv,') = ',X1
          PRINT *,' xarrm(',irm,') = ',X2
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDDO
C
      RETURN
      END
C*********************************************************************
