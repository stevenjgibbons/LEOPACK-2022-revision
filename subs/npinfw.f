C*********************************************************************
C subroutine NPlot INteger File Write ********************************
C            --    --      -    -     ********************************
C Steve Gibbons Wed Dec  8 08:17:15 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Writes out the integers file to be read by the nplot program.      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LU        : Logical unit number of file.                       C
C     IWR       : Write flag. = 2 for caution, =3 for overwrite.     C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     NH        : Number of harmonics.                               C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : Name of output file.                               C
C                                                                    C
C     CHVMFF    : Velocity/magnetic field select.                    C
C                 Set to either 'VEL' or 'MAG'                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NPINFW( LU, FNAME, IWR, MHT, MHL, MHM, NH, CHVMFF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LU, IWR, MHT( * ), MHL( * ), MHM( * ), NH
      CHARACTER *(3) CHVMFF
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IL, IM, IICS, NT, IPOL, ITOR, IH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IWR.NE.2 .AND. IWR.NE.3 ) THEN
        PRINT *,' Subroutine NPINFW.'
        PRINT *,' IWR = ', IWR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Let's check validity of CHVMFF
C     .
      IF ( CHVMFF.NE.'VEL' .AND. CHVMFF.NE.'vel' .AND.
     1     CHVMFF.NE.'Vel' .AND. CHVMFF.NE.'MAG' .AND.
     2     CHVMFF.NE.'mag' .AND. CHVMFF.NE.'Mag'  ) THEN
        PRINT *,' Subroutine NPINFW.'
        PRINT *,' CHVMFF = ', CHVMFF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( CHVMFF.EQ.'VEL' .OR. CHVMFF.EQ.'vel' .OR.
     1     CHVMFF.EQ.'Vel' ) THEN
        IPOL = 1
        ITOR = 2
      ENDIF
C     .
      IF ( CHVMFF.EQ.'MAG' .OR. CHVMFF.EQ.'mag' .OR.
     1     CHVMFF.EQ.'Mag' ) THEN
        IPOL = 4
        ITOR = 5
      ENDIF
C
C Open file
C
      CALL FOPEN ( LU, FNAME, IWR )
C
      DO IH = 1, NH
        NT = MHT( IH )
        IF ( NT.NE.IPOL .AND. NT.NE.ITOR ) GOTO 50
        IL = MHL( IH )
C
C (note that the definition of ICS in nplot is the
C opposite to that in my original codes)
C
        IF ( MHM( IH ).LT.0 ) THEN
          IM = -MHM( IH )
          IICS = 1
        ELSE
          IM = MHM( IH )
          IICS = 2
        ENDIF
        IF ( NT.EQ.IPOL ) NT = 1
        IF ( NT.EQ.ITOR ) NT = 2
         WRITE ( LU, 791 ) NT, IL, IM, IICS
 50   CONTINUE
      ENDDO
C
C Close file.
C
      CALL FCLOSE ( LU, FNAME, 'Error in closing file.')
C
 791  FORMAT (4i2)
      RETURN
      END
C*********************************************************************
