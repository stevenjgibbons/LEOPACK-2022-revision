C*********************************************************************
C subroutine Meridian Section Polar Array Value Find *****************
C            -        -       -     -     -     -    *****************
C Steve Gibbons Mon Jul 10 17:29:41 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given NRAD values of R in the array RVALS and NTHE values of       C
C theta in the array TVALS (given in RADIANS!), MSPAVF will fill     C
C the array MSPAV( dim[ NRAD, NTHE ]) with values corresponding to   C
C a given component of the solution. Now it is customary for         C
C RVALS( 1 ) to be zero; and all j with RVALS( j ).lt.XARR( 1 )      C
C will have MSPAV filled with VMISS                                  C
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
C IAS is the axisymmetric switch. If IAS = 1, then only harmonics    C
C with m=0 are counted. Otherwise all are counted.                   C
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
C     NTHE      : Number of theta values for evaluation.             C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     IAS       : Axisymmetric switch. = 1 --> only m=0 harmonics.   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     PHI       : Longitude (radians).                               C
C                                                                    C
C     VEC       : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C     XARR      : Dim ( * ) but length atleast NR. Location of       C
C                  radial grid nodes.                                C
C                                                                    C
C     RVALS     : Radius values for evaluation. Dim (NRAD)           C
C     TVALS     : Theta values for evaluation in radians. Dim (NTHE) C
C                                                                    C
C     WORK1     : Dimension ( NNDS ). Work array.                    C
C     WORK2     : Dimension ( NNDS ). Work array.                    C
C     COEFM     : Dimension ( NNDS, NNDS ). Work array.              C
C                                                                    C
C     WORK3     : Dimension ( NTHE ). Work array.                    C
C     PA        : Work array dim ( (LH+1)*(LH+2)/2 , NTHE )          C
C     DPA       : Work array dim ( (LH+1)*(LH+2)/2 , NTHE )          C
C     VMISS     : Value to be used if out of range.                  C
C                                                                    C
C     MSPAV     : Dim (NRAD, NTHE ) output with evaluated functions. C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MSPAVF( ICOMP, INARR, NNDS, IWORK, MHT, MHL, MHM, 
     1                   NRAD, NTHE, LH, IAS, PHI, VEC, XARR, RVALS,
     2                   TVALS, WORK1, WORK2, COEFM, WORK3, PA, DPA,
     3                   VMISS, MSPAV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER ICOMP, INARR( * ), NNDS, IWORK( NNDS ), MHT( * ),
     1        MHL( * ), MHM( * ), NRAD, NTHE, LH, IAS
      DOUBLE PRECISION PHI, VEC( * ), XARR( * ), RVALS( NRAD ),
     1                 TVALS( NTHE ), WORK1( NNDS ),
     2                 WORK2( NNDS ), COEFM( NNDS, NNDS ), VMISS,
     3                 WORK3( NTHE ), PA( (LH+1)*(LH+2)/2, NTHE),
     4                 DPA( (LH+1)*(LH+2)/2, NTHE)
      DOUBLE PRECISION MSPAV( NRAD, NTHE)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRAD, ITHE, NR, NH, ITYPE, IPOL, ITOR, L, M,
     1        INDKLM, IH
      DOUBLE PRECISION THETA, RAD, RI, RO, PHTRM, DPHTRM, P, DP,
     1                 DLOW, SINTH, VALUE, TERM
      LOGICAL OVEL, OMAG, OVEC
      PARAMETER ( DLOW = 2.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check that ICOMP is valid ...
C
      IF ( ICOMP.LT.1 .OR. ICOMP.GT.7 ) THEN
        PRINT *,' Subroutine MSPAVF.'
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
C Now calculate all the Legendre polynomials etc.
C First fill WORK3 with the cosines of the theta values
C
      DO ITHE = 1, NTHE
        THETA = TVALS( ITHE )
        WORK3( ITHE ) = DCOS( THETA )
      ENDDO
C
C Now call SCHNLA to calculate P_l^m and derivatives
C
      CALL SCHNLA( PA, DPA, WORK3, LH, NTHE )
C
      DO IRAD = 1, NRAD
        RAD = RVALS( IRAD )
        IF ( RAD.LT.RI .OR. RAD.GT.RO .OR. RAD.LT.DLOW ) THEN
          DO ITHE = 1, NTHE
            MSPAV( IRAD, ITHE ) = VMISS
          ENDDO
        ELSE
          DO ITHE = 1, NTHE
            VALUE = 0.0d0
            THETA = TVALS( ITHE )
            SINTH = DSIN( THETA )
C           .
C           . Loop around spherical harmonics
C           .
            DO IH = 1, NH
              ITYPE = MHT( IH )
C             .
C             . Go straight onto next harmonic if type is wrong
C             .
              IF ( OVEC .AND. ITYPE.NE.IPOL .AND.
     1                           ITYPE.NE.ITOR ) GOTO 50
              IF ( ICOMP.EQ.7 .AND. ITYPE.NE.3 ) GOTO 50
              IF ( ICOMP.EQ.1 .AND. ITYPE.NE.IPOL ) GOTO 50
              IF ( ICOMP.EQ.4 .AND. ITYPE.NE.IPOL ) GOTO 50
C             .
C             . OK - so we definitely need to consider this harmonic
C             .
              CALL SVRINT( RAD, VEC, XARR, INARR, IH, NNDS, WORK1,
     1               IWORK, WORK2, COEFM )
C             .
C             . WORK1( 1 ) now contains f(r)
C             . WORK1( 2 ) now contains df(r)/dr etc.
C             .
              L        = MHL( IH )
              IF ( L.GT.LH ) THEN
                PRINT *,' Subroutine MSPAVF.'
                PRINT *,' MHL(',IH,') = ', L,' LH = ', LH
                PRINT *,' Program aborted.'
                STOP
              ENDIF
C             .
              IF (     MHM( IH ).LT.0     ) THEN
                M      = -MHM( IH )
                PHTRM  = DSIN( M*PHI )
                DPHTRM = DBLE( M )*DCOS( M*PHI )
              ELSE
                M      = MHM( IH )
                PHTRM  = DCOS( M*PHI )
                DPHTRM = DBLE( M )*(-1.0d0)*DSIN( M*PHI )
              ENDIF
              INDKLM   = L*(L+1)/2+M+1
              IF ( IAS.EQ.1 .AND. M.NE.0 ) GOTO 50
C             .
C             . First do case of temperature harmonic
C             .
              IF ( ICOMP.EQ.7 ) THEN
                P      = PA( INDKLM, ITHE )
                VALUE  = VALUE + WORK1( 1 )*P*PHTRM
                GOTO 50
              ENDIF
C             .
C             . Case of radial vector component
C             .
              IF ( ICOMP.EQ.1 .OR. ICOMP.EQ.4 ) THEN
                P      = PA( INDKLM, ITHE )
                TERM   = DBLE( L*L + L )*WORK1( 1 )/RAD
                VALUE  = VALUE + TERM*P*PHTRM
                GOTO 50
              ENDIF
C             .
C             . Case of theta vector component
C             .
              IF ( ICOMP.EQ.2 .OR. ICOMP.EQ.5 ) THEN
                IF ( ITYPE.EQ.IPOL ) THEN
                  TERM   = WORK1( 1 )/RAD + WORK1( 2 )
                  DP     = DPA( INDKLM, ITHE )
                  VALUE  = VALUE + TERM*DP*PHTRM
                ELSE
                  IF ( DABS( SINTH ).LT.DLOW ) THEN
                    PRINT *,' Subroutine MSPAVF.'
                    PRINT *,' SINTH = ',SINTH
                    PRINT *,' Division be zero imminent.'
                    PRINT *,' Program aborted.'
                    STOP
                  ENDIF
                  TERM   = WORK1( 1 )/SINTH
                  P      = PA( INDKLM, ITHE )
                  VALUE  = VALUE + TERM*P*DPHTRM
                ENDIF
                GOTO 50
              ENDIF
C             .
C             . Case of phi vector component
C             .
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
                  P      = PA( INDKLM, ITHE )
                  VALUE  = VALUE + TERM*P*DPHTRM/SINTH
                ELSE
                  TERM   = (-1.0d0)*WORK1( 1 )
                  DP     = DPA( INDKLM, ITHE )
                  VALUE  = VALUE + TERM*DP*PHTRM
                ENDIF
                GOTO 50
              ENDIF
C             .
 50         CONTINUE
            ENDDO
C           .
C           . Finish loop around harmonics
C           .
          MSPAV( IRAD, ITHE ) = VALUE
          ENDDO
C         .
C         . End loop around ITHE
C         .
        ENDIF
      ENDDO
C     .
C     . End loop around IRAD
C     .
      RETURN
      END
C*********************************************************************
