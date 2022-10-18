C*********************************************************************
C subroutine HarMonic Type Array Check *******************************
C            -  -     -    -     -     *******************************
C Steve Gibbons Mon Apr  3 12:34:38 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C The element ih of the MHT array indicates to which scalar function C
C the ih^{th} harmonic belongs.                                      C
C                                                                    C
C MHT( ih ) = 1 --> harmonic is poloidal velocity                    C
C MHT( ih ) = 2 --> harmonic is toroidal velocity                    C 
C MHT( ih ) = 3 --> harmonic is temperature                          C
C MHT( ih ) = 4 --> harmonic is poloidal magnetic field              C
C MHT( ih ) = 5 --> harmonic is toroidal magnetic field              C
C                                                                    C
C We may want to ensure that a vector contains only harmonics of     C
C a certain type. HMTAC returns an error message if MHT contains     C
C an element which is not allowed by the specification CHOPT.        C
C                                                                    C
C  CHOPT = 'VEL' only allows 1 and 2.                                C
C  CHOPT = 'VTF' only allows 1, 2 and 3.                             C
C  CHOPT = 'MAG' only allows 4 and 5.                                C
C                                                                    C
C  Other options may be allowed at a later date.                     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Number of vector spherical harmonics.              C
C     MHT       : Array length ( * ) - atleast length NHMAX          C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHOPT     : *(*) Option. See above.                            C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE HMTAC( NH, MHT, CHOPT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, MHT( * )
      CHARACTER *(*) CHOPT
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH, ITYPE
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( CHOPT(1:3).EQ.'VEL' ) GOTO 50
      IF ( CHOPT(1:3).EQ.'VTF' ) GOTO 50
      IF ( CHOPT(1:3).EQ.'MAG' ) GOTO 50
      PRINT *,' Subroutine HMTAC.'
      PRINT *,' CHOPT = ', CHOPT
      PRINT *,' Option not yet allowed.'
      PRINT *,' PROGRAM ABORTED.'
      STOP
C
 50   CONTINUE
C
      DO IH = 1, NH
C       .
        ITYPE = MHT( IH )
C       .
        IF ( CHOPT(1:3).EQ.'VEL' .AND. ITYPE.NE.1 .AND.
     1       ITYPE.NE.2                                 ) THEN
          PRINT *,' Subroutine HMTAC.'
          PRINT *,' Harmonic ',IH,' has type ',ITYPE
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
        IF ( CHOPT(1:3).EQ.'VTF' .AND. ITYPE.NE.1 .AND.
     1       ITYPE.NE.2   .AND. ITYPE.NE.3              ) THEN
          PRINT *,' Subroutine HMTAC.'
          PRINT *,' Harmonic ',IH,' has type ',ITYPE
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
        IF ( CHOPT(1:3).EQ.'MAG' .AND. ITYPE.NE.4 .AND.
     1       ITYPE.NE.5                                 ) THEN
          PRINT *,' Subroutine HMTAC.'
          PRINT *,' Harmonic ',IH,' has type ',ITYPE
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
      ENDDO
C
      RETURN
      END
C*********************************************************************
