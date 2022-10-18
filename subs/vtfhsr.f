C*********************************************************************
C subroutine Velocity and Temp. Floquet Harmonic Selection Routine ***
C            -            -     -       -        -         -       ***
C Steve Gibbons Thu Jun  8 13:45:08 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C This routine assumes that we have a steady velocity with modes     C
C 0 up to MMAX through increments of MINC.                           C
C VTFHSR will choose all harmonics of the correct symmetry such      C
C that m is only selected if (m-mfloq) is divisible by minc          C
C or (m+mfloq) is divisible by minc.                                 C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ISYM      : ISYM.eq.1 --> only equatorially symmetric harm.s   C
C                 ISYM.eq.2 --> only equatorially symmetric harm.s   C
C                 ISYM.eq.3 --> all equatorially symmetries          C
C                                                                    C
C     NH        : Output only integer giving the total number of     C
C                  harmonics selected.                               C
C                                                                    C
C     MLOW      : Lowest wavenumber, m.                              C
C     MINC      : Increment in wavenumber, m.                        C
C     MMAX      : Highest wavenumber, m.                             C
C                                                                    C
C     NHMAX     : Maximum number of harmonics permitted. If the      C
C                  specified parameters demand more harmonics than   C
C                   permitted by this bound then HMINDA calculates   C
C                    how many harmonics are necessary and then       C
C                     aborts with an appropriate message.            C
C                                                                    C
C     MHT       : MHT( ih ) contains itype for harmonic 'ih'         C
C                                                                    C
C     MHL       : MHL( ih ) contains degree, l, for harmonic 'ih'    C
C                                                                    C
C     MHM       : MHM( ih ) contains order, m, for harmonic 'ih' if  C
C                  ih has cos m phi dependency and (-m) if ih has    C
C                   sin m phi dependency.                            C
C                                                                    C
C     LH        : Maximum requested degree, l, of a spherical harm.  C
C                                                                    C
C     LHMAX     : Global maximum permitted degree, l, of a spherical C
C                  harmonic.                                         C
C                                                                    C
C     MFLOQ     : Floquet integer.                                   C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VTFHSR( ISYM, NH, MLOW, MINC, MMAX, NHMAX, MHT, MHL,
     1                   MHM, LH, LHMAX, MFLOQ )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER ISYM, NH, MLOW, MINC, MMAX, NHMAX, MHT( * ), MHL( * ),
     1        MHM( * ), LH, LHMAX, MFLOQ
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L, M, LMIN, ICS, ITYPE, IS, M2, MFAC
      LOGICAL OK, OES, OEA
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      IF ( LH.GT.LHMAX ) THEN
        PRINT *,' Subroutine VTFHSR.'
        PRINT *,' LH    = ', LH
        PRINT *,' LHMAX = ', LHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
c
      IF ( ISYM.NE.1 .AND. ISYM.NE.2 .AND. ISYM.NE.3 ) THEN
        PRINT *,' Subroutine VTFHSR.'
        PRINT *,' ISYM = ', ISYM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
c
      IF ( MINC.EQ.0 ) THEN
        PRINT *,' Subroutine VTFHSR.'
        PRINT *,' MINC = ', MINC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
c
      NH = 0
      OK = .TRUE.
c
      OES = .FALSE.
      OEA = .FALSE.
      IF ( ISYM.EQ.1 .OR. ISYM.EQ.3 ) OES = .TRUE.
      IF ( ISYM.EQ.2 .OR. ISYM.EQ.3 ) OEA = .TRUE.
c
      DO M = MLOW, MMAX
        MFAC = MOD( M, MINC )
        IF ( ( .NOT. MFLOQ.EQ.MFAC )             .AND.
     1       ( .NOT. (MFAC+MFLOQ-MINC).EQ.0 )   ) GOTO 50
        DO ITYPE = 1, 3
          IF ( ITYPE.EQ.1 ) M2 = 0
          IF ( ITYPE.EQ.2 ) M2 = 1
          IF ( ITYPE.EQ.3 ) M2 = 0
          DO IS = 1, 2
            IF ( .NOT. OES  .AND. IS.EQ.1 ) GOTO 49
            IF ( .NOT. OEA  .AND. IS.EQ.2 ) GOTO 49
            LMIN = M
            IF ( M.EQ.0 .AND. ( ITYPE.NE.3 .OR. IS.NE.1 ) ) LMIN = 1
            DO L = LMIN, LH
              IF ( IS.EQ.1 .AND. MOD( (L-M), 2 ).NE.M2 ) GOTO 48
              IF ( IS.EQ.2 .AND. MOD( (L-M), 2 ).EQ.M2 ) GOTO 48
              DO ICS = 1, 2
                IF ( M.EQ.0 .AND. ICS.EQ.2 ) GOTO 48
c               . ok this harmonic DOES go in
c               . (provided we have enough room)
c               .
                NH = NH + 1
                IF ( NH.GT.NHMAX ) OK = .FALSE.
                IF ( OK ) MHT( NH ) = ITYPE
                IF ( OK ) MHL( NH ) = L
                IF ( OK .AND. ICS.EQ.1 ) MHM( NH ) = M
                IF ( OK .AND. ICS.EQ.2 ) MHM( NH ) = -M
c               .
              ENDDO
 48         CONTINUE
            ENDDO
 49       CONTINUE
          ENDDO
        ENDDO
 50   CONTINUE
      ENDDO
c
      IF ( OK ) RETURN
      PRINT *,' Subroutine VTFHSR. Your specifications'
      PRINT *,' require ',NH,' harmonics. Maximum was set'
      PRINT *,' at ',NHMAX,'. Program aborted.'
c
      STOP
      END
C*********************************************************************
