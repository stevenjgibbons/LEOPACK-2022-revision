C*********************************************************************
C subroutine HarMonic INDex Allocation subroutine ********************
C            -  -     ---   -                     ********************
C Steve Gibbons Sat Sep 25 11:52:41 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Fills the arrays MHT, MHL and MHM according to user-specified      C
C conditions on symmetry (equatorial and rotational) dimension and   C
C choice of variables.                                               C
C                                                                    C
C HMINDA will allocate indices to these arrays in ascending order    C
C of modes with symmetries separate, regardless of type of function. C
C                                                                    C
C The maximum degree of spherical harmonic, LH, can be set           C
C differently for each kind of function and they can have different  C
C specified equatorial symmetries. HMINDA does not check on whether  C
C the supplied combination is a physically sensible or meaningful    C
C choice.                                                            C
C                                                                    C
C  For example, set LH( 1 ) = 8      (poloidal velocity)             C
C                   LH( 2 ) = 8      (toroidal velocity)             C
C                   LH( 3 ) = 8      (temperature)                   C
C                   LH( 4 ) = 0      (poloidal magnetic field)       C
C                   LH( 5 ) = 0      (toroidal magnetic field)       C
C                                                                    C
C  Set NMODES = 2 and MMODES( 1 ) = 0                                C
C                     MMODES( 2 ) = 4                                C
C                                                                    C
C  Set ISYM( 1 ) = 3  This will select both equatorially symmetric   C
C      ISYM( 2 ) = 3  and equatorially anti-symmetric                C
C      ISYM( 3 ) = 3    harmonics for first 3 components             C
C      ISYM( 4 ) = 0    There are not to be any magnetic field       C
C      ISYM( 5 ) = 0      harmonics.                                 C
C                                                                    C
C  HMINDA will loop around the wavenumbers (imode)                   C
C  We then loop around the equatorial symmetries first ES then EA    C
C  If either is not wanted, we move straight onto the next.          C
C  It will then loop around ITYPE ( = 1, 2, 3, 4 , 5)                C
C  If LH( itype ).lt.MMODES( imode ) then we move onto the           C
C  next type.                                                        C
C  Finally, if MMODES( imode ).ne.0, we loop around ICS = cos        C
C                                                     to ICS = sin   C
C                                                                    C
C  In the following, ITYPE will always refer to                      C
C                                                                    C
C         ITYPE = 1 for a poloidal velocity harmonic.                C
C         ITYPE = 2 for a toroidal velocity harmonic.                C
C         ITYPE = 3 for a temperature harmonic.                      C
C         ITYPE = 4 for a poloidal magnetic field harmonic.          C
C         ITYPE = 5 for a toroidal magnetic field harmonic.          C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LHARR     : Dim ( 5 ). LHARR( itype ) contains maximum degree  C
C                 of spherical harmonic for type ITYPE.              C
C                  LHARR( itype ) not referred to if                 C
C                   ISYMA( itype ).eq.0.                             C
C                                                                    C
C     ISYMA     : Dim ( 5 ). Specifies equatorial symmetry property  C
C                 of type ITYPE.                                     C
C                                                                    C
C                ISYMA( itype ).eq.0 --> we are simply not including C
C                                        any harmonics of this type. C
C                ISYMA( itype ).eq.1 --> only equatorially symmetric C
C                ISYMA( itype ).eq.2 --> only equatorially anti-sym. C
C                ISYMA( itype ).eq.3 --> all eq. symmetries          C
C                                                                    C
C     NMODES    : Number of wavenumbers ( m ) to be included in the  C
C                   solution vector.                                 C
C                                                                    C
C     MMODES    : Integer array Dim.( NMODES )                       C
C                  NMODES( im ) contains m for this mode             C
C                                                                    C
C     NH        : Output only integer giving the total number of     C
C                  harmonics selected.                               C
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
C     LHMAX     : Global maximum permitted degree, l, of a spherical C
C                  harmonic.                                         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE HMINDA( LHARR, ISYMA, NMODES, MMODES, NH, NHMAX,
     1                   MHT, MHL, MHM, LHMAX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LHARR( 5 ), ISYMA( 5 ), NMODES, MMODES( NMODES ),
     1        NH, NHMAX, MHT( * ), MHL( * ), MHM( * ), LHMAX
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L, M, IMODE, LMIN, LMAX, ISYM, ICS, M2, ITYPE
      LOGICAL OES( 5 ), OEA( 5 ), OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      DO ITYPE = 1, 5
C
        OES( ITYPE ) = .FALSE.
        OEA( ITYPE ) = .FALSE.
        IF (           ISYMA( ITYPE ).EQ.1 .OR.
     1                 ISYMA( ITYPE ).EQ.3      )
     2                         OES( ITYPE ) = .TRUE.
        IF (           ISYMA( ITYPE ).EQ.2 .OR.
     1                 ISYMA( ITYPE ).EQ.3      )
     2                         OEA( ITYPE ) = .TRUE.
C
        IF ( LHARR( ITYPE ).GT.LHMAX ) THEN
          PRINT *,' Subroutine HMINDA.'
          PRINT *,' LHARR(',ITYPE,') = ',LHARR( ITYPE )
          PRINT *,' LHMAX = ', LHMAX
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
      ENDDO
c
      NH = 0
      OK = .TRUE.
c
      DO IMODE = 1, NMODES
        M = MMODES( IMODE )
C
C isym = 1 will consider the equatorially symmetric
C harmonics, isym = 2 the equatorially antisymmetric
C
        DO ISYM = 1, 2
          DO ITYPE = 1, 5
C
            IF ( ITYPE.EQ.1 ) M2 = 0
            IF ( ITYPE.EQ.2 ) M2 = 1
            IF ( ITYPE.EQ.3 ) M2 = 0
            IF ( ITYPE.EQ.4 ) M2 = 0
            IF ( ITYPE.EQ.5 ) M2 = 1
C
C If an ES function has (l-m) even then M2 is 0
C If an ES function has (l-m) odd then M2 is 1
C
            IF ( ISYMA( ITYPE ).EQ.0 ) GOTO 50
            IF ( .NOT. OES( ITYPE ) .AND. ISYM.EQ.1 ) GOTO 50
            IF ( .NOT. OEA( ITYPE ) .AND. ISYM.EQ.2 ) GOTO 50
C
C set LMIN and LMAX - first checking for the monopole term
C
            LMIN = MMODES( IMODE )
            IF (     LMIN.EQ.0 .AND. 
     1             ( ISYM.NE.1 .OR. ITYPE.NE.3 ) ) LMIN = 1
            LMAX = LHARR( ITYPE )
C
C Loop around L from LMIN, LMAX
C
            DO L = LMIN, LMAX
c             .
c             . ignore if our symmetry is not satisfied
c             .
              IF (     ISYM.EQ.1 .AND.
     1              MOD( (L-M), 2 ).NE.M2 
     2                                       ) GOTO 49
              IF (     ISYM.EQ.2 .AND.
     1              MOD( (L-M), 2 ).EQ.M2 
     2                                       ) GOTO 49
c             .
              DO ICS = 1, 2
                IF ( M.EQ.0 .AND. ICS.EQ.2 ) GOTO 49
c               .
c               . ok this harmonic DOES go in
c               . (provided we have enough room)
c               .
                NH = NH + 1
                IF ( NH.GT.NHMAX ) OK = .FALSE.
                IF ( OK ) MHT( NH ) = ITYPE
                IF ( OK ) MHL( NH ) = L
                IF ( OK .AND. ICS.EQ.1 ) MHM( NH ) = M
                IF ( OK .AND. ICS.EQ.2 ) MHM( NH ) = -M
              ENDDO
c             .
 49         CONTINUE
            ENDDO
C
 50       CONTINUE
          ENDDO
        ENDDO
      ENDDO
c
      IF ( OK ) RETURN
      PRINT *,' Subroutine HMINDA. Your specifications'
      PRINT *,' require ',NH,' harmonics. Maximum was set'
      PRINT *,' at ',NHMAX,'. Program aborted.'
c
      STOP
      END
C*********************************************************************
