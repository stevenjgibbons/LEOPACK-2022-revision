C*********************************************************************
C subroutine EQuatorial SECtion POLAR array EVALuation ***************
C            --         ---     -----       ----       ***************
C Steve Gibbons Tue Jan 23 08:39:42 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C F is a real array of dimensions ( NRAD, NTHE ).                    C
C EQ_SEC_POLAR_EVAL adds to it the values of a function evaluated    C
C in the equatorial section with the = THE.                          C
C The coordinates are defined by a COMMON block:                     C
C                                                                    C
C      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST   C
C                                                                    C
C It consists of two INTEGERs                                        C
C NRAD (number of equally spaced radial grid nodes) and              C
C NTHE (number of equally spaced phi grid nodes), and four REAL      C
C variables RFIRST, RLAST, TFIRST, TLAST                             C
C                                                                    C
C RFIRST is the inner boundary radius                                C
C RLAST is the outer boundary radius                                 C
C                                                                    C
C TFIRST is the first phi point, in radians.                         C
C TLAST is the last phi point, in radians.                           C
C                                                                    C
C The value of theta is kept constant at THE (specified in radians)  C
C Although the most natural value to use is THE = pi/2, any value    C
C of theta strictly greater than 0 and strictly less than pi will    C
C work.                                                              C
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
C  11 :  F_theta streamfunction (magnetic field)                     C
C  12 :  F_theta streamfunction (velocity field)                     C
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
C     NRAD      : Number of radial values for evaluation.            C
C     NTHE      : Number of phi values for evaluation.               C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     IAS       : Axisymmetric switch. = 1 --> only m=0 harmonics.   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     THE       : Co-latitude (radians).                             C
C                 ( pi/2 for equatorial section is most likely ).    C
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
C     VMISS     : Value to be used if out of range.                  C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     ESPAV     : Dim ( NRAD, NTHE) output with evaluated functions. C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE EQ_SEC_POLAR_EVAL( ICOMP, INARR, NNDS, IWORK,
     1           MHT, MHL, MHM, LH, IAS, THE, VEC, XARR, WORK1,
     2           WORK2, COEFM, VMISS, ESPAV )
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER ICOMP, INARR( * ), NNDS, IWORK( NNDS ), MHT( * ),
     1        MHL( * ), MHM( * ), LH, IAS
      DOUBLE PRECISION THE, VEC( * ), XARR( * ), WORK1( NNDS ),
     1                 WORK2( NNDS ), COEFM( NNDS, NNDS ), VMISS
      REAL             ESPAV( NRAD, NTHE )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRAD, ITHE, NR, NH, ITYPE, IPOL, ITOR, L, M, IH
      DOUBLE PRECISION PHI, RAD, RI, RO, PHTRM, DPHTRM, P, DP,
     1                 DLOW, SINTH, COSTH, TERM, SHMPLG, SHDPLG
      LOGICAL OVEL, OMAG, OVEC
      PARAMETER ( DLOW = 2.0d-3 )
      REAL    R1, T1, DELR1, DELT1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check that ICOMP is valid ...
C
      IF ( ICOMP.LT.1 .OR. ICOMP.GT.12 ) THEN
        PRINT *,' Subroutine EQ_SEC_POLAR_EVAL.'
        PRINT *,' ICOMP = ', ICOMP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NR = INARR( 2 )
      NH = INARR( 3 )
      RI = XARR(  1 )
      RO = XARR( NR )
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
      DELR1 = (RLAST-RFIRST)/REAL(NRAD-1)
      DELT1 = (TLAST-TFIRST)/REAL(NTHE-1)
C
      COSTH = DCOS( THE )
      SINTH = DSIN( THE )
      IF ( DABS( SINTH ).LT.DLOW ) THEN
c       PRINT *,' Subroutine EQ_SEC_POLAR_EVAL.'
c       PRINT *,' SINTH = ',SINTH
c       PRINT *,' Division be zero imminent.'
c       PRINT *,' Program aborted.'
c       STOP
        SINTH = DLOW
      ENDIF
      R1 = RFIRST - DELR1
      DO IRAD = 1, NRAD
c       R1  = RFIRST + (RLAST-RFIRST)*REAL(IRAD-1)/REAL(NRAD-1)
        R1  = R1 + DELR1
        RAD = DBLE( R1 )
        IF ( (RI-RAD).GT.DLOW .OR. (RAD-RO).GT.DLOW
     1                        .OR. RAD.LT.DLOW ) THEN
          DO ITHE = 1, NTHE
            ESPAV( IRAD, ITHE ) = REAL( VMISS )
          ENDDO
        ELSE
          IF ( RAD.LT.RI ) RAD = RI
          IF ( RAD.GT.RO ) RAD = RO
C         .
C         . Loop around spherical harmonics
C         .
          DO IH = 1, NH
            ITYPE = MHT( IH )
C           .
C           . Go straight onto next harmonic if type is wrong
C           .
            IF ( OVEC .AND. ITYPE.NE.IPOL .AND.
     1                         ITYPE.NE.ITOR ) GOTO 50
            IF ( ICOMP.EQ.7 .AND. ITYPE.NE.3 ) GOTO 50
            IF ( ICOMP.EQ.1 .AND. ITYPE.NE.IPOL ) GOTO 50
            IF ( ICOMP.EQ.4 .AND. ITYPE.NE.IPOL ) GOTO 50
            IF ( ICOMP.EQ.11 .AND. ITYPE.NE.4 ) GOTO 50
            IF ( ICOMP.EQ.12 .AND. ITYPE.NE.1 ) GOTO 50
C           .
C           . OK - so we definitely need to consider this harmonic
C           .
            CALL SVRINT( RAD, VEC, XARR, INARR, IH, NNDS, WORK1,
     1             IWORK, WORK2, COEFM )
C           .
C           . WORK1( 1 ) now contains f(r)
C           . WORK1( 2 ) now contains df(r)/dr etc.
C           .
            L        = MHL( IH )
            IF ( L.GT.LH ) THEN
              PRINT *,' Subroutine EQ_SEC_POLAR_EVAL.'
              PRINT *,' MHL(',IH,') = ', L,' LH = ', LH
              PRINT *,' Program aborted.'
              STOP
            ENDIF
C           .
            IF (     MHM( IH ).LT.0     ) THEN
              M      = -MHM( IH )
            ELSE
              M      = MHM( IH )
            ENDIF
            IF ( IAS.EQ.1 .AND. M.NE.0 ) GOTO 50
C           .
C           . First do case of temperature harmonic
C           .
            IF ( ICOMP.EQ.7 ) THEN
              P      = SHMPLG ( L, M, COSTH )
              T1     = TFIRST - DELT1
              DO ITHE = 1, NTHE
c               T1 = TFIRST + (TLAST-TFIRST)*REAL(ITHE-1)/REAL(NTHE-1)
                T1 = T1 + DELT1
                PHI = DBLE( T1 )
                IF (     MHM( IH ).LT.0     ) THEN
                  PHTRM  = DSIN( M*PHI )
                ELSE
                  PHTRM  = DCOS( M*PHI )
                ENDIF
                ESPAV( IRAD, ITHE ) = ESPAV( IRAD, ITHE ) +
     1                          REAL(  WORK1( 1 )*P*PHTRM  )
              ENDDO
              GOTO 50
            ENDIF
C           .
C           . Case of radial vector component
C           .
            IF ( ICOMP.EQ.1 .OR. ICOMP.EQ.4 ) THEN
              P      = SHMPLG ( L, M, COSTH )
              TERM   = DBLE( L*L + L )*WORK1( 1 )/RAD
              T1     = TFIRST - DELT1
              DO ITHE = 1, NTHE
c               T1 = TFIRST + (TLAST-TFIRST)*REAL(ITHE-1)/REAL(NTHE-1)
                T1 = T1 + DELT1
                PHI = DBLE( T1 )
                IF (     MHM( IH ).LT.0     ) THEN
                  PHTRM  = DSIN( M*PHI )
                ELSE
                  PHTRM  = DCOS( M*PHI )
                ENDIF
                ESPAV( IRAD, ITHE ) = ESPAV( IRAD, ITHE ) +
     1                           REAL(  TERM*P*PHTRM  )
              ENDDO
              GOTO 50
            ENDIF
C           .
C           . Case of F_theta streamfunction
C           .
            IF ( ICOMP.EQ.11 .OR. ICOMP.EQ.12 ) THEN
              P      = SHMPLG ( L, M, COSTH )
              TERM   = RAD*WORK1( 1 )/SINTH
              T1     = TFIRST - DELT1
              DO ITHE = 1, NTHE
                T1 = T1 + DELT1
                PHI = DBLE( T1 )
                IF (     MHM( IH ).LT.0     ) THEN
                  DPHTRM = DBLE( M )*DCOS( M*PHI )
                ELSE
                  DPHTRM = DBLE( M )*(-1.0d0)*DSIN( M*PHI )
                ENDIF
                ESPAV( IRAD, ITHE ) = ESPAV( IRAD, ITHE ) +
     1                           REAL(  TERM*P*DPHTRM  )
              ENDDO
              GOTO 50
            ENDIF
C           .
C           . Case of theta vector component
C           .
            IF ( ICOMP.EQ.2 .OR. ICOMP.EQ.5 ) THEN
              IF ( ITYPE.EQ.IPOL ) THEN
                TERM   = WORK1( 1 )/RAD + WORK1( 2 )
                DP     = SHDPLG ( L, M, COSTH )
                T1     = TFIRST - DELT1
                DO ITHE = 1, NTHE
c               T1 = TFIRST + (TLAST-TFIRST)*REAL(ITHE-1)/REAL(NTHE-1)
                T1 = T1 + DELT1
                  PHI = DBLE( T1 )
                  IF (     MHM( IH ).LT.0     ) THEN
                    PHTRM  = DSIN( M*PHI )
                  ELSE
                    PHTRM  = DCOS( M*PHI )
                  ENDIF
                  ESPAV( IRAD, ITHE ) = ESPAV( IRAD, ITHE ) +
     1                           REAL(  TERM*DP*PHTRM  )
                ENDDO
              ELSE
                TERM   = WORK1( 1 )/SINTH
                P      = SHMPLG ( L, M, COSTH )
                T1     = TFIRST - DELT1
                DO ITHE = 1, NTHE
c               T1 = TFIRST + (TLAST-TFIRST)*REAL(ITHE-1)/REAL(NTHE-1)
                T1 = T1 + DELT1
                  PHI = DBLE( T1 )
                  IF (     MHM( IH ).LT.0     ) THEN
                    DPHTRM = DBLE( M )*DCOS( M*PHI )
                  ELSE
                    DPHTRM = DBLE( M )*(-1.0d0)*DSIN( M*PHI )
                  ENDIF
                  ESPAV( IRAD, ITHE ) = ESPAV( IRAD, ITHE ) +
     1                            REAL( TERM*P*DPHTRM  )
                ENDDO
              ENDIF
              GOTO 50
            ENDIF
C           .
C           . Case of phi vector component
C           .
            IF ( ICOMP.EQ.3 .OR. ICOMP.EQ.6 ) THEN
              IF ( ITYPE.EQ.IPOL ) THEN
                TERM   = WORK1( 1 )/RAD + WORK1( 2 )
                P      = SHMPLG ( L, M, COSTH )
                T1     = TFIRST - DELT1
                DO ITHE = 1, NTHE
c               T1 = TFIRST + (TLAST-TFIRST)*REAL(ITHE-1)/REAL(NTHE-1)
                T1 = T1 + DELT1
                  PHI = DBLE( T1 )
                  IF (     MHM( IH ).LT.0     ) THEN
                    DPHTRM = DBLE( M )*DCOS( M*PHI )
                  ELSE
                    DPHTRM = DBLE( M )*(-1.0d0)*DSIN( M*PHI )
                  ENDIF
                  ESPAV( IRAD, ITHE ) = ESPAV( IRAD, ITHE ) +
     1                         REAL( TERM*P*DPHTRM/SINTH  )
                ENDDO
              ELSE
                TERM   = (-1.0d0)*WORK1( 1 )
                DP     = SHDPLG ( L, M, COSTH )
                T1     = TFIRST - DELT1
                DO ITHE = 1, NTHE
c               T1 = TFIRST + (TLAST-TFIRST)*REAL(ITHE-1)/REAL(NTHE-1)
                T1 = T1 + DELT1
                  PHI = DBLE( T1 )
                  IF (     MHM( IH ).LT.0     ) THEN
                    PHTRM  = DSIN( M*PHI )
                  ELSE
                    PHTRM  = DCOS( M*PHI )
                  ENDIF
                  ESPAV( IRAD, ITHE ) = ESPAV( IRAD, ITHE ) +
     1                        REAL(  TERM*DP*PHTRM  )
                ENDDO
              ENDIF
              GOTO 50
            ENDIF
C           .
 50       CONTINUE
          ENDDO
C         .
C         . Finish loop around harmonics
C         .
        ENDIF
      ENDDO
C     .
C     . End loop around IRAD
C     .
      RETURN
      END
C*********************************************************************
