C*********************************************************************
C subroutine Constant Z Polar Array Value Find ***********************
C            -        - -     -     -     -    ***********************
C Steve Gibbons Tue Aug 22 16:27:48 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given NESS values of S in the array SVALS and NPHI values of       C
C phi in the array PVALS (given in RADIANS!), CZPAVF will fill       C
C the array CZPAV( dim[ NESS, NPHI ]) with values corresponding to   C
C a given component of the solution. Now it is customary for         C
C                                                                    C
C The value of z (height above equatorial plane) is kept constant at C
C ZED.                                                               C
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
C     NESS      : Number of s values for evaluation.                 C
C                 's' is the distance from rotation axis.            C
C     NPHI      : Number of phi values for evaluation.               C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     IAS       : Axisymmetric switch. = 1 --> only m=0 harmonics.   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     ZED       : Height above equatorial plane.                     C
C                                                                    C
C     VEC       : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C     XARR      : Dim ( * ) but length atleast NR. Location of       C
C                  radial grid nodes.                                C
C                                                                    C
C     SVALS     : Radius values for evaluation. Dim (NESS)           C
C     PVALS     : Theta values for evaluation in radians. Dim (NPHI) C
C                                                                    C
C     WORK1     : Dimension ( NNDS ). Work array.                    C
C     WORK2     : Dimension ( NNDS ). Work array.                    C
C     COEFM     : Dimension ( NNDS, NNDS ). Work array.              C
C                                                                    C
C     VMISS     : Value to be used if out of range.                  C
C                                                                    C
C     CZPAV     : Dim ( NESS, NPHI) output with evaluated functions. C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CZPAVF( ICOMP, INARR, NNDS, IWORK, MHT, MHL, MHM, 
     1                   NESS, NPHI, LH, ZED, VEC, XARR, SVALS,
     2                   PVALS, WORK1, WORK2, COEFM, VMISS, CZPAV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER ICOMP, INARR( * ), NNDS, IWORK( NNDS ), MHT( * ),
     1        MHL( * ), MHM( * ), NESS, NPHI, LH
      DOUBLE PRECISION ZED, VEC( * ), XARR( * ), SVALS( NESS ),
     1                 PVALS( NPHI ), WORK1( NNDS ),
     2                 WORK2( NNDS ), COEFM( NNDS, NNDS ),
     3                 CZPAV( NESS, NPHI ), VMISS
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IESS, IPHI, NR, NH, ITYPE, IPOL, ITOR, L, M, IH
      DOUBLE PRECISION PHI, RAD, RI, RO, PHTRM, DPHTRM, P, DP, THE,
     1                 DLOW, SINTH, COSTH, TERM, SHMPLG, SHDPLG, S
      LOGICAL OVEL, OMAG, OVEC
      PARAMETER ( DLOW = 2.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check that ICOMP is valid ...
C
      IF ( ICOMP.LT.1 .OR. ICOMP.GT.7 ) THEN
        PRINT *,' Subroutine CZPAVF.'
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
      IF ( DABS( ZED ).GT.RO ) THEN
        PRINT *,' Subroutine CZPAVF.'
        PRINT *,' RO = ', RO,' ZED = ', ZED
        PRINT *,' Program aborted.'
        STOP
      ENDIF
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
      DO IESS = 1, NESS
        S     = SVALS( IESS )
        THE   = DATAN2( S, ZED )
        COSTH = DCOS( THE )
        SINTH = DSIN( THE )
        RAD   = DSQRT( S*S + ZED*ZED )
        IF ( RAD.LT.RI .OR. RAD.GT.RO .OR. RAD.LT.DLOW ) THEN
          DO IPHI = 1, NPHI
            CZPAV( IESS, IPHI ) = VMISS
          ENDDO
        ELSE
          DO IPHI = 1, NPHI
            CZPAV( IESS, IPHI ) = 0.0d0
          ENDDO
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
              PRINT *,' Subroutine CZPAVF.'
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
C           .
C           . First do case of temperature harmonic
C           .
            IF ( ICOMP.EQ.7 ) THEN
              P      = SHMPLG ( L, M, COSTH )
              DO IPHI = 1, NPHI
                PHI = PVALS( IPHI )
                IF (     MHM( IH ).LT.0     ) THEN
                  PHTRM  = DSIN( M*PHI )
                ELSE
                  PHTRM  = DCOS( M*PHI )
                ENDIF
                CZPAV( IESS, IPHI ) = CZPAV( IESS, IPHI ) +
     1                                 WORK1( 1 )*P*PHTRM
              ENDDO
              GOTO 50
            ENDIF
C           .
C           . Case of radial vector component
C           .
            IF ( ICOMP.EQ.1 .OR. ICOMP.EQ.4 ) THEN
              P      = SHMPLG ( L, M, COSTH )
              TERM   = DBLE( L*L + L )*WORK1( 1 )/RAD
              DO IPHI = 1, NPHI
                PHI = PVALS( IPHI )
                IF (     MHM( IH ).LT.0     ) THEN
                  PHTRM  = DSIN( M*PHI )
                ELSE
                  PHTRM  = DCOS( M*PHI )
                ENDIF
                CZPAV( IESS, IPHI ) = CZPAV( IESS, IPHI ) +
     1                                  TERM*P*PHTRM
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
                DO IPHI = 1, NPHI
                  PHI = PVALS( IPHI )
                  IF (     MHM( IH ).LT.0     ) THEN
                    PHTRM  = DSIN( M*PHI )
                  ELSE
                    PHTRM  = DCOS( M*PHI )
                  ENDIF
                  CZPAV( IESS, IPHI ) = CZPAV( IESS, IPHI ) +
     1                                  TERM*DP*PHTRM
                ENDDO
              ELSE
                IF ( DABS( SINTH ).LT.DLOW ) THEN
                  PRINT *,' Subroutine MSPAVF.'
                  PRINT *,' SINTH = ',SINTH
                  PRINT *,' Division be zero imminent.'
                  PRINT *,' Program aborted.'
                  STOP
                ENDIF
                TERM   = WORK1( 1 )/SINTH
                P      = SHMPLG ( L, M, COSTH )
                DO IPHI = 1, NPHI
                  PHI = PVALS( IPHI )
                  IF (     MHM( IH ).LT.0     ) THEN
                    DPHTRM = DBLE( M )*DCOS( M*PHI )
                  ELSE
                    DPHTRM = DBLE( M )*(-1.0d0)*DSIN( M*PHI )
                  ENDIF
                  CZPAV( IESS, IPHI ) = CZPAV( IESS, IPHI ) +
     1                                  TERM*P*DPHTRM
                ENDDO
              ENDIF
              GOTO 50
            ENDIF
C           .
C           . Case of phi vector component
C           .
            IF ( ICOMP.EQ.3 .OR. ICOMP.EQ.6 ) THEN
              IF ( ITYPE.EQ.IPOL ) THEN
                IF ( DABS( SINTH ).LT.DLOW ) THEN
                  PRINT *,' Subroutine MSPAVF.'
                  PRINT *,' SINTH = ',SINTH
                  PRINT *,' Division be zero imminent.'
                  PRINT *,' Program aborted.'
                  STOP
                ENDIF
                TERM   = WORK1( 1 )/RAD + WORK1( 2 )
                P      = SHMPLG ( L, M, COSTH )
                DO IPHI = 1, NPHI
                  PHI = PVALS( IPHI )
                  IF (     MHM( IH ).LT.0     ) THEN
                    DPHTRM = DBLE( M )*DCOS( M*PHI )
                  ELSE
                    DPHTRM = DBLE( M )*(-1.0d0)*DSIN( M*PHI )
                  ENDIF
                  CZPAV( IESS, IPHI ) = CZPAV( IESS, IPHI ) +
     1                               TERM*P*DPHTRM/SINTH
                ENDDO
              ELSE
                TERM   = (-1.0d0)*WORK1( 1 )
                DP     = SHDPLG ( L, M, COSTH )
                DO IPHI = 1, NPHI
                  PHI = PVALS( IPHI )
                  IF (     MHM( IH ).LT.0     ) THEN
                    PHTRM  = DSIN( M*PHI )
                  ELSE
                    PHTRM  = DCOS( M*PHI )
                  ENDIF
                  CZPAV( IESS, IPHI ) = CZPAV( IESS, IPHI ) +
     1                               TERM*DP*PHTRM
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
C     . End loop around IESS
C     .
      RETURN
      END
C*********************************************************************
