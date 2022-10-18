C*********************************************************************
C subroutine MMAX FinD ***********************************************
C            ---- -  - ***********************************************
C Steve Gibbons Tue Feb 15 09:48:45 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we have NH harmonics with harmonic ih having a wavenumber, m,   C
C given as IABS( MHM( ih ) ) then MMAXFD returns MMAX (the highest   C
C wavenumber used) and M0 which is the lowest non zero integer which C
C divides every wavenumber, m.                                       C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NH        : Number of spherical harmonics.                     C
C     MHM       : Dim ( * ). MHM( ih ) = m for cos m phi dependence. C
C                            MHM( ih ) = -m for sin m phi dependence C
C                                                                    C
C     MMAX      : Maximum value of m.                                C
C     M0        : Minimum value of non-zero m such that              C
C                  MOD( m, M0 ) = 0 for all m.                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MMAXFD( NH, MHM, MMAX, M0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, MHM( * ), MMAX, M0
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH, M
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      MMAX = 0
      M0   = 10000
C
      DO IH = 1, NH
        M = IABS( MHM( IH ) )
        IF ( M.GT.MMAX ) MMAX = M
        IF ( M.NE.0 .AND. M.LT.M0 ) M0 = M
      ENDDO
C     .
      DO IH = 1, NH
        M = IABS( MHM( IH ) )
        IF ( MOD( M, M0 ).NE.0 ) THEN
          M0 = 1
          RETURN
        ENDIF
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
