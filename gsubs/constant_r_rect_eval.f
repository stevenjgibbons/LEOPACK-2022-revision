C*********************************************************************
C subroutine CONSTANT R RECTangular EVALuate *************************
C            -------- - ----        ----     *************************
C Steve Gibbons Tue Jan 23 13:23:38 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C                                                                    C
C F is a real array of dimensions ( NPHI, NTHE ).                    C
C This routine ADDS to F the values of a function at constant radius C
C RAD. Theta ranges from THETA1 to THETA2, and PHI ranges from       C
C PHI1 to PHI2                                                       C
C                                                                    C
C All of the solution is contained within the vector VEC with        C
C format given by INARR (see INDFUN)                                 C
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
C The component is given by the integer flag ICOMP where             C
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
C   8 :  Heat flux: -(dT/dr) where T is as in 7.                     C
C   9 :  Upwelling: -(dv_r/dr)                                       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NPHI      : Number of phi points.                              C
C     NTHE      : Number of theta points.                            C
C                                                                    C
C     ICOMP     : See key above.                                     C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NR. Number of radial grid nodes.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C     NNDS      : Number of nodes to be used in interpolation.       C
C                 Must be atleast 3.                                 C
C     IWORK     : Dimension ( NNDS ). Work array.                    C
C                                                                    C
C     MHT       : Dim (*). See above.                                C
C     MHL       : Dim (*). See above.                                C
C     MHM       : Dim (*). See above.                                C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     IAS       : Axisymmetric switch. = 1 --> only m=0 harmonics.   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Radius at which functions are to be evaluated.     C
C                                                                    C
C     VEC       : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C     XARR      : Dim ( * ) but length atleast NR. Location of       C
C                  radial grid nodes.                                C
C                                                                    C
C     WORK1     : Dimension ( NNDS ). Work array.                    C
C     WORK2     : Dimension ( NNDS ). Work array.                    C
C     COEFM     : Dimension ( NNDS, NNDS ). Work array.              C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     THETA1    : Lowest value of theta.                             C
C     THETA2    : Highest value of theta.                            C
C                                                                    C
C     PHI1      : Lowest value of phi.                               C
C     PHI2      : Highest value of phi.                              C
C                                                                    C
C     CRPAV     : Dim ( NPHI, NTHE) output with evaluated functions. C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CONSTANT_R_RECT_EVAL( ICOMP, INARR, NNDS, IWORK,
     1          MHT, MHL, MHM, LH, IAS, RAD, VEC, XARR, WORK1, WORK2,
     2          COEFM, CRPAV, NPHI, NTHE, THETA1, THETA2,
     3          PHI1, PHI2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER ICOMP, INARR( * ), NNDS, IWORK( NNDS ), MHT( * ),
     1        MHL( * ), MHM( * ), LH, IAS, NPHI, NTHE
      DOUBLE PRECISION RAD, VEC( * ), XARR( * ), WORK1( NNDS ),
     1                 WORK2( NNDS ), COEFM( NNDS, NNDS )
      REAL             CRPAV( NPHI, NTHE ), THETA1, THETA2,
     1                 PHI1, PHI2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IPHI, ITHE, NH, ITYPE, IPOL, ITOR, L, M, IH
      DOUBLE PRECISION PHI, PHTRM, DPHTRM, P, DP, THE,
     1                 DLOW, SINTH, COSTH, TERM, SHMPLG, SHDPLG
      LOGICAL OVEL, OMAG, OVEC
      REAL             P1, T1, DELP1, DELT1
      PARAMETER ( DLOW = 2.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check that ICOMP is valid ...
C
      IF ( ICOMP.LT.1 .OR. ICOMP.GT.9 ) THEN
        PRINT *,' Subroutine CONSTANT_R_RECT_EVAL.'
        PRINT *,' ICOMP = ', ICOMP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NH = INARR( 3 )
C
      IPOL   = 1
      ITOR   = 2
C
      OVEL   = .FALSE.
      OVEC   = .FALSE.
      OMAG   = .FALSE.
C
      IF ( ICOMP.EQ.1 .OR. ICOMP.EQ.2 .OR. ICOMP.EQ.3 ) THEN
        OVEL = .TRUE.
      ENDIF
C
      IF ( ICOMP.EQ.4 .OR. ICOMP.EQ.5 .OR. ICOMP.EQ.6 ) THEN
        OMAG = .TRUE.
      ENDIF
C
      IF ( OMAG .OR. OVEL ) OVEC = .TRUE.
C
      IF ( OMAG ) THEN
        IPOL = 4
        ITOR = 5
      ENDIF
C
      DELP1 = (PHI2-PHI1)/REAL(NPHI-1)
      DELT1 = (THETA2-THETA1)/REAL(NTHE-1)
      T1    = THETA1 - DELT1
      DO ITHE = 1, NTHE
c       T1  = THETA1 + (THETA2-THETA1)*REAL(ITHE-1)/REAL(NTHE-1)
        T1  = T1 + DELT1
        THE   = DBLE( T1 )
        COSTH = DCOS( THE )
        SINTH = DSIN( THE )
C
        IF ( DABS( SINTH ).LT.DLOW ) THEN
c         PRINT *,' Subroutine CONSTANT_R_RECT_EVAL.'
c         PRINT *,' SINTH = ',SINTH
c         PRINT *,' Division be zero imminent.'
c         PRINT *,' Program aborted.'
c         STOP
          SINTH = DLOW
        ENDIF
C       .
C       . Loop around spherical harmonics
C       .
        DO IH = 1, NH
          ITYPE = MHT( IH )
C         .
C         . Go straight onto next harmonic if type is wrong
C         .
          IF ( OVEC .AND. ITYPE.NE.IPOL .AND.
     1                         ITYPE.NE.ITOR ) GOTO 50
          IF ( ICOMP.EQ.7 .AND. ITYPE.NE.3 ) GOTO 50
          IF ( ICOMP.EQ.8 .AND. ITYPE.NE.3 ) GOTO 50
          IF ( ICOMP.EQ.9 .AND. ITYPE.NE.1 ) GOTO 50
          IF ( ICOMP.EQ.1 .AND. ITYPE.NE.IPOL ) GOTO 50
          IF ( ICOMP.EQ.4 .AND. ITYPE.NE.IPOL ) GOTO 50
C         .
C         . OK - so we definitely need to consider this harmonic
C         .
          CALL SVRINT( RAD, VEC, XARR, INARR, IH, NNDS, WORK1,
     1             IWORK, WORK2, COEFM )
C         .
C         . WORK1( 1 ) now contains f(r)
C         . WORK1( 2 ) now contains df(r)/dr etc.
C         .
          L        = MHL( IH )
          IF ( L.GT.LH ) THEN
            PRINT *,' Subroutine CONSTANT_R_RECT_EVAL.'
            PRINT *,' MHL(',IH,') = ', L,' LH = ', LH
            PRINT *,' Program aborted.'
            STOP
          ENDIF
C         .
          IF (     MHM( IH ).LT.0     ) THEN
            M      = -MHM( IH )
          ELSE
            M      = MHM( IH )
          ENDIF
          IF ( IAS.EQ.1 .AND. M.NE.0 ) GOTO 50
C         .
C         . First do case of temperature harmonic
C         .
          IF ( ICOMP.EQ.7 ) THEN
            P      = SHMPLG ( L, M, COSTH )
            P1     = PHI1 - DELP1
            DO IPHI = 1, NPHI
c             P1 = PHI1 + (PHI2-PHI1)*REAL(IPHI-1)/REAL(NPHI-1)
              P1 = P1 + DELP1
              PHI = DBLE( P1 )
              IF (     MHM( IH ).LT.0     ) THEN
                PHTRM  = DSIN( M*PHI )
              ELSE
                PHTRM  = DCOS( M*PHI )
              ENDIF
              CRPAV( IPHI, ITHE ) = CRPAV( IPHI, ITHE ) +
     1                       REAL(     WORK1( 1 )*P*PHTRM  )
            ENDDO
            GOTO 50
          ENDIF
C         .
C         . Now do case of temperature harmonic (heat flux)
C         .
          IF ( ICOMP.EQ.8 ) THEN
            P      = SHMPLG ( L, M, COSTH )
            P1     = PHI1 - DELP1
            DO IPHI = 1, NPHI
c             P1 = PHI1 + (PHI2-PHI1)*REAL(IPHI-1)/REAL(NPHI-1)
              P1 = P1 + DELP1
              PHI = DBLE( P1 )
              IF (     MHM( IH ).LT.0     ) THEN
                PHTRM  = DSIN( M*PHI )
              ELSE
                PHTRM  = DCOS( M*PHI )
              ENDIF
              CRPAV( IPHI, ITHE ) = CRPAV( IPHI, ITHE ) -
     1                       REAL(     WORK1( 2 )*P*PHTRM  )
            ENDDO
            GOTO 50
          ENDIF
C         .
C         . Case of radial vector component
C         .
          IF ( ICOMP.EQ.1 .OR. ICOMP.EQ.4 ) THEN
            P      = SHMPLG ( L, M, COSTH )
            TERM   = DBLE( L*L + L )*WORK1( 1 )/RAD
            P1     = PHI1 - DELP1
            DO IPHI = 1, NPHI
c             P1 = PHI1 + (PHI2-PHI1)*REAL(IPHI-1)/REAL(NPHI-1)
              P1 = P1 + DELP1
              PHI = DBLE( P1 )
              IF (     MHM( IH ).LT.0     ) THEN
                PHTRM  = DSIN( M*PHI )
              ELSE
                PHTRM  = DCOS( M*PHI )
              ENDIF
              CRPAV( IPHI, ITHE ) = CRPAV( IPHI, ITHE ) +
     1                          REAL(   TERM*P*PHTRM  )
            ENDDO
            GOTO 50
          ENDIF
C         .
C         . Case of upwelling
C         .
          IF ( ICOMP.EQ.9 ) THEN
            P      = SHMPLG ( L, M, COSTH )
            TERM   = DBLE( L*L + L )*(WORK1( 2 )-WORK1( 1 )/RAD)/RAD
            P1     = PHI1 - DELP1
            DO IPHI = 1, NPHI
c             P1 = PHI1 + (PHI2-PHI1)*REAL(IPHI-1)/REAL(NPHI-1)
              P1 = P1 + DELP1
              PHI = DBLE( P1 )
              IF (     MHM( IH ).LT.0     ) THEN
                PHTRM  = DSIN( M*PHI )
              ELSE
                PHTRM  = DCOS( M*PHI )
              ENDIF
              CRPAV( IPHI, ITHE ) = CRPAV( IPHI, ITHE ) -
     1                          REAL(   TERM*P*PHTRM  )
            ENDDO
            GOTO 50
          ENDIF
C         .
C         . Case of theta vector component
C         .
          IF ( ICOMP.EQ.2 .OR. ICOMP.EQ.5 ) THEN
            IF ( ITYPE.EQ.IPOL ) THEN
              TERM   = WORK1( 1 )/RAD + WORK1( 2 )
              DP     = SHDPLG ( L, M, COSTH )
              P1     = PHI1 - DELP1
              DO IPHI = 1, NPHI
c             P1 = PHI1 + (PHI2-PHI1)*REAL(IPHI-1)/REAL(NPHI-1)
              P1 = P1 + DELP1
                PHI = DBLE( P1 )
                IF (     MHM( IH ).LT.0     ) THEN
                  PHTRM  = DSIN( M*PHI )
                ELSE
                  PHTRM  = DCOS( M*PHI )
                ENDIF
                CRPAV( IPHI, ITHE ) = CRPAV( IPHI, ITHE ) +
     1                          REAL(   TERM*DP*PHTRM  )
              ENDDO
            ELSE
              TERM   = WORK1( 1 )/SINTH
              P      = SHMPLG ( L, M, COSTH )
              P1     = PHI1 - DELP1
              DO IPHI = 1, NPHI
c             P1 = PHI1 + (PHI2-PHI1)*REAL(IPHI-1)/REAL(NPHI-1)
              P1 = P1 + DELP1
                PHI = DBLE( P1 )
                IF (     MHM( IH ).LT.0     ) THEN
                  DPHTRM = DBLE( M )*DCOS( M*PHI )
                ELSE
                  DPHTRM = DBLE( M )*(-1.0d0)*DSIN( M*PHI )
                ENDIF
                CRPAV( IPHI, ITHE ) = CRPAV( IPHI, ITHE ) +
     1                          REAL(   TERM*P*DPHTRM  )
              ENDDO
            ENDIF
            GOTO 50
          ENDIF
C         .
C         . Case of phi vector component
C         .
          IF ( ICOMP.EQ.3 .OR. ICOMP.EQ.6 ) THEN
            IF ( ITYPE.EQ.IPOL ) THEN
              TERM   = WORK1( 1 )/RAD + WORK1( 2 )
              P      = SHMPLG ( L, M, COSTH )
              P1     = PHI1 - DELP1
              DO IPHI = 1, NPHI
c             P1 = PHI1 + (PHI2-PHI1)*REAL(IPHI-1)/REAL(NPHI-1)
              P1 = P1 + DELP1
                PHI = DBLE( P1 )
                IF (     MHM( IH ).LT.0     ) THEN
                  DPHTRM = DBLE( M )*DCOS( M*PHI )
                ELSE
                  DPHTRM = DBLE( M )*(-1.0d0)*DSIN( M*PHI )
                ENDIF
                CRPAV( IPHI, ITHE ) = CRPAV( IPHI, ITHE ) +
     1                       REAL(   TERM*P*DPHTRM/SINTH  )
              ENDDO
            ELSE
              TERM   = (-1.0d0)*WORK1( 1 )
              DP     = SHDPLG ( L, M, COSTH )
              P1     = PHI1 - DELP1
              DO IPHI = 1, NPHI
c             P1 = PHI1 + (PHI2-PHI1)*REAL(IPHI-1)/REAL(NPHI-1)
              P1 = P1 + DELP1
                PHI = DBLE( P1 )
                IF (     MHM( IH ).LT.0     ) THEN
                  PHTRM  = DSIN( M*PHI )
                ELSE
                  PHTRM  = DCOS( M*PHI )
                ENDIF
                CRPAV( IPHI, ITHE ) = CRPAV( IPHI, ITHE ) +
     1                       REAL(   TERM*DP*PHTRM  )
              ENDDO
            ENDIF
            GOTO 50
          ENDIF
C         .
 50     CONTINUE
        ENDDO
C       .
C       . Finish loop around harmonics
C       .
      ENDDO
C     .
C     . End loop around ITHE
C     .
      RETURN
      END
C*********************************************************************
