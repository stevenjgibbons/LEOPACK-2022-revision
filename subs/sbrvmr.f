C*********************************************************************
C subroutine Solid Body Rotation Vector Modification Routine *********
C            -     -    -        -      -            -       *********
C Steve Gibbons Mon Oct 15 09:43:10 WEST 2001                        C
C____________________________________________________________________C
C                                                                    C
C VEC0 is a double precision vector array which has spherical        C
C harmonic radial functions defined by the integer arrays IN0, MT0,  C
C ML0, MM0. The radial spacings are defined by the array XARR.       C
C SBRVMR sets a second array VEC0M to REY*VEC0 for all elements of   C
C VEC0 except for elements corresponding to the solid body rotation  C
C harmonic MT0( ih ) = 2, ML0( ih ) = 1, MM0( ih ) = 0, which is     C
C modified to VEC0M = REY*VEC0 + REYSBR*RAD                          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IN0       : See INDFUN.                                        C
C                                                                    C
C     MT0       : Type of harmonic formed. Only MT0( i ) = 2 is      C
C                 likely to result from this routine - toroidal vel. C
C                                                                    C
C     ML0       : Spherical harmonic degree l.                       C
C     MM0       : Spherical harmonic order, m. (Set to zero).        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC0      : Dim ( * ). Atleast dimension NR*IN0( 3 )           C
C                 (Normalised vector array)                          C
C                                                                    C
C     VEC0M     : Dim ( * ). Atleast dimension NR*IN0( 3 )           C
C                 (Modified vector array)                            C
C                                                                    C
C     XARR      : Dim ( NR ). Values of radius.                      C
C                                                                    C
C     REY       : Multiplication factor for total vector             C
C     REYSBR    : Multiplication factor for solid body rotation      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SBRVMR( IN0, MT0, ML0, MM0, VEC0, VEC0M, XARR,
     1                   REY, REYSBR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          IN0( * ), MT0( * ), ML0( * ), MM0( * )
      DOUBLE PRECISION VEC0( * ), VEC0M( * ), XARR( * ), REY, REYSBR
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IH, IR, NH, NR, IND, INDFUN
      EXTERNAL         INDFUN
      DOUBLE PRECISION RAD
      LOGICAL          OSBR
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      OSBR = .FALSE.
C
      NR   = IN0( 2 )
      NH   = IN0( 3 )
C
      DO IH = 1, NH
C       .
C       . First case of solid body rotation harmonic
C       .
        IF (        MT0( IH ).EQ.2     .AND.
     1              ML0( IH ).EQ.1     .AND.
     2              MM0( IH ).EQ.0                  )   THEN
C         .
          IF ( OSBR ) THEN
            PRINT *,' Subroutine SBRVMR.'
            PRINT *,' Solid body rotation harmonic already found.'
            PRINT *,' Program aborted.'
            STOP
          ENDIF
          OSBR = .TRUE.
          DO IR = 1, NR
            IND = INDFUN( IR, IH, IN0 )
            RAD = XARR( IR )
            VEC0M( IND ) = REY*VEC0( IND ) + REYSBR*RAD
          ENDDO
C         .
        ELSE
C         .
C         . This harmonic does not involve a solid body rotation
C         .
          DO IR = 1, NR
            IND = INDFUN( IR, IH, IN0 )
            VEC0M( IND ) = REY*VEC0( IND )
          ENDDO
C         .
        ENDIF
      ENDDO
C
      IF ( .NOT. OSBR ) THEN
        PRINT *,' Subroutine SBRVMR.'
        PRINT *,' Solid body rotation harmonic not found.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
