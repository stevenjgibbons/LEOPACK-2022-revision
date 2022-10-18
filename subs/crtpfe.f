C*********************************************************************
C double precision function Chebyshev R, Theta, Phi Function Eval. ***
C                           -         -  -      -   -        -     ***
C Steve Gibbons Sat Feb 19 11:40:52 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C This routine is identical to RTPFCE except that it acts upon       C
C a set of Chebyshev polynomial coefficients instead of finite       C
C difference function values.                                        C
C                                                                    C
C Given a solution vector, VEC, with format given by the INARR array C
C (see INDFUN) harmonics type given by MHT (a harmonic, ih, has      C
C                                                                    C
C      MHT( ih ) = 1 --> poloidal velocity                           C
C      MHT( ih ) = 2 --> toroidal velocity                           C
C      MHT( ih ) = 3 --> temperature                                 C
C      MHT( ih ) = 4 --> poloidal magnetic field                     C
C      MHT( ih ) = 5 --> toroidal magnetic field )                   C
C                                                                    C
C MHL( ih ) is the spherical harmonic degree, l.                     C
C MHM( ih ) is the spherical harmonic oder, m, for cos m phi         C
C           dependence and, -m, for sin m phi dependence.            C
C                                                                    C
C then CRTPFE will evaluate the value of a given component selected  C
C by the integer ICOMP at the co-ordinates RAD (distance from the    C
C centre of the sphere), THETA (standard spherical co-latitude in    C
C radians) and PHI (standard spherical longitude in radians).        C
C                                                                    C
C ICOMP can take one of the following values:-                       C
C                                                                    C
C   1 :  v_r - radial component of velocity                          C
C   2 :  v_{theta} theta component of velocity                       C
C   3 :  v_{phi} phicomponent of velocity                            C
C   4 :  B_r - radial component of magnetic field.                   C
C   5 :  B_{theta} theta component of magnetic field.                C
C   6 :  B_{phi} phicomponent of magnetic field.                     C
C   7 :  Theta - this is the homogeneous part of the temperature.    C
C        If an inhomogeneous part exists, it must be calculated      C
C        separately.                                                 C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ICOMP     : See key above.                                     C
C     INARR     : Integer array dimension ( * ).                     C
C                 See CHINDF.                                        C
C                                                                    C
C     MHT       : Dim (*). See above.                                C
C     MHL       : Dim (*). See above.                                C
C     MHM       : Dim (*). See above.                                C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Radial position.                                   C
C     THE       : Colatitude (radians).                              C
C     PHI       : Longitude (radians).                               C
C                                                                    C
C     CCV       : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION CRTPFE( ICOMP, INARR, MHT, MHL, MHM, RAD,
     1                 THE, PHI, CCV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER ICOMP, INARR( * ), MHT( * ), MHL( * ), MHM( * )
      DOUBLE PRECISION RAD, THE, PHI, CCV( * ), CRTPFE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH, IPOL, ITOR, ITYPE, L, M, NH
      DOUBLE PRECISION COSTH, P, DP, PHTRM, DPHTRM, DLOW, TERM,
     1                 SINTH, SHMPLG, SHDPLG, Y, DY, D2Y
      LOGICAL OVEL, OMAG, OTEM, OVEC
      PARAMETER ( DLOW = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NH    = INARR( 3 )
C
      COSTH = DCOS( THE )
      SINTH = DSQRT( 1.0d0 - COSTH*COSTH )
C
      CRTPFE = 0.0d0
C
      IPOL   = 1
      ITOR   = 2
C
      OVEL   = .FALSE.
      OVEC   = .FALSE.
      OMAG   = .FALSE.
      OTEM   = .FALSE.
C
      IF ( ICOMP.EQ.1 .OR. ICOMP.EQ.2 .OR. ICOMP.EQ.3 ) THEN
        OVEL = .TRUE.
      ENDIF
C
      IF ( ICOMP.EQ.4 .OR. ICOMP.EQ.5 .OR. ICOMP.EQ.6 ) THEN
        OMAG = .TRUE.
      ENDIF
C
      IF ( ICOMP.EQ.7 ) THEN
        OTEM = .TRUE.
      ENDIF
C
      IF ( (.NOT. OMAG) .AND. (.NOT. OVEL) .AND. (.NOT. OTEM) ) THEN
        PRINT *,' Function CRTPFE.'
        PRINT *,' ICOMP = ', ICOMP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( OMAG .OR. OVEL ) OVEC = .TRUE.
C
      IF ( OMAG ) THEN
        IPOL = 4
        ITOR = 5
      ENDIF
C
      DO IH = 1, NH
C       .
        ITYPE = MHT( IH )
C       .
C       . Go straight onto next harmonic if type is wrong
C       .
        IF ( OVEC .AND. ITYPE.NE.IPOL .AND. ITYPE.NE.ITOR ) GOTO 50
        IF ( ICOMP.EQ.7 .AND. ITYPE.NE.3 ) GOTO 50
        IF ( ICOMP.EQ.1 .AND. ITYPE.NE.IPOL ) GOTO 50
        IF ( ICOMP.EQ.4 .AND. ITYPE.NE.IPOL ) GOTO 50
C       .
C       . OK - we definitely need to consider this harmonic
C       .
c       CALL SVRINT( RAD, VEC, XARR, INARR, IH, NNDS, WORK1,
c    1               IWORK, WORK2, COEFM )
        CALL CCVSPE( INARR, IH, RAD, CCV, Y, DY, D2Y )
C       .
C       . Y now contains f(r)
C       . DY now contains df(r)/dr etc.
C       .
        L        = MHL( IH )
        IF (     MHM( IH ).LT.0     ) THEN
          M      = -MHM( IH )
          PHTRM  = DSIN( M*PHI )
          DPHTRM = DBLE( M )*DCOS( M*PHI )
        ELSE
          M      = MHM( IH )
          PHTRM  = DCOS( M*PHI )
          DPHTRM = DBLE( M )*(-1.0d0)*DSIN( M*PHI )
        ENDIF
C       .
        IF ( DABS( PHTRM ).LT.DLOW ) GOTO 50
C       .
C       . First do case of temperature harmonic
C       .
        IF ( ICOMP.EQ.7 ) THEN
          P      = SHMPLG ( L, M, COSTH )
          CRTPFE = CRTPFE + Y*P*PHTRM
          GOTO 50
        ENDIF
C       .
C       . Case of radial vector component
C       .
        IF ( ICOMP.EQ.1 .OR. ICOMP.EQ.4 ) THEN
          P      = SHMPLG ( L, M, COSTH )
          TERM   = DBLE( L*L + L )*Y/RAD
          CRTPFE = CRTPFE + TERM*P*PHTRM
          GOTO 50
        ENDIF
C       .
C       . Case of theta vector component
C       .
        IF ( ICOMP.EQ.2 .OR. ICOMP.EQ.5 ) THEN
          IF ( ITYPE.EQ.IPOL ) THEN
            TERM   = Y/RAD + DY
            DP     = SHDPLG ( L, M, COSTH )
            CRTPFE = CRTPFE + TERM*DP*PHTRM
          ELSE
            IF ( DABS( SINTH ).LT.DLOW ) THEN
              PRINT *,' Function CRTPFE.'
              PRINT *,' SINTH = ',SINTH
              PRINT *,' Division be zero imminent.'
              PRINT *,' Program aborted.'
              STOP
            ENDIF
            TERM   = Y/SINTH
            P      = SHMPLG ( L, M, COSTH )
            CRTPFE = CRTPFE + TERM*P*DPHTRM
          ENDIF
          GOTO 50
        ENDIF
C       .
C       . Case of phi vector component
C       .
        IF ( ICOMP.EQ.3 .OR. ICOMP.EQ.6 ) THEN
          IF ( ITYPE.EQ.IPOL ) THEN
            IF ( DABS( SINTH ).LT.DLOW ) THEN
              PRINT *,' Function CRTPFE.'
              PRINT *,' SINTH = ',SINTH
              PRINT *,' Division be zero imminent.'
              PRINT *,' Program aborted.'
              STOP
            ENDIF
            TERM   = Y/RAD + DY
            P      = SHMPLG ( L, M, COSTH )
            CRTPFE = CRTPFE + TERM*P*DPHTRM/SINTH
          ELSE
            TERM   = (-1.0d0)*Y
            DP     = SHDPLG ( L, M, COSTH )
            CRTPFE = CRTPFE + TERM*DP*PHTRM
          ENDIF
          GOTO 50
        ENDIF
C       .
 50   CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************
