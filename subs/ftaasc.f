C*********************************************************************
C subroutine Fourier Transform Array Azimuthal Symmetry Complete *****
C            -       -         -     -         -        -        *****
C Steve Gibbons Sat Jul 22 11:16:34 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C The array FTF ( 2*NPHPR ) is zero except for the elements          C
C ind = ( 2*i - 1 ) for i = 1, nphpv.                                C
C                                                                    C
C Now, NPHPR must be a multiple of NPHPV.                            C
C                                                                    C
C FTAASC completes the array FTF in the same pattern but             C
C repeated for every NPHPV grid points.                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NPHPR     : Total number of phi points needed to perform       C
C                 Fourier transform.                                 C
C     NPHPV     : Number of phi points required to fully recreate    C
C                 the behaviour of the solution in a subspace of the C
C                 total physical space.                              C
C                                                                    C
C                 NPHPR must be a multiple, M0 * NPHPV.              C
C                                                                    C
C                 This is checked for.                               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     FTF       : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPR )                        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FTAASC( NPHPR, NPHPV, FTF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NPHPR, NPHPV
      DOUBLE PRECISION  FTF( 2*NPHPR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER           M0, I, IND, J, IND2
      DOUBLE PRECISION  VAL
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NPHPV.LE.0 .OR. NPHPR.LE.0 ) THEN
        PRINT *,' Subroutine FTAASC.'
        PRINT *,' NPHPV = ', NPHPV
        PRINT *,' NPHPR = ', NPHPR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      M0 = NPHPR/NPHPV
C
C Check that NPHPV divides NPHPR ...
C
      IF ( NPHPV*M0.NE.NPHPR ) THEN
        PRINT *,' Subroutine FTAASC.'
        PRINT *,' M0    = ', M0
        PRINT *,' NPHPV = ', NPHPV
        PRINT *,' NPHPR = ', NPHPR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Exit if solution already has lowest possible symmetry.
C
      IF ( M0.EQ.1 ) RETURN
C
C So M0 is atleast 2.
C
      DO I = 1, NPHPV
        IND  = 2*I - 1
        VAL  = FTF( IND )
        IND2 = IND
        DO J = 2, M0
          IND2        = IND2 + NPHPV + NPHPV
          FTF( IND2 ) = VAL
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
