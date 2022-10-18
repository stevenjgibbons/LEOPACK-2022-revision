C*********************************************************************
C subroutine Solution Vector Longitudinal Shift Subroutine ***********
C            -        -      -            -     -          ***********
C Steve Gibbons Mon Mar 20 14:55:39 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Let SV be a solution vector containing radial functions of the     C
C standard spherical harmonics with respect to longitude phi         C
C                                                                    C
C SVLSS returns SSV as the same vector but as a function of          C
C (phi + tau) where tau is an angle of longitude.                    C
C                                                                    C
C The only information necessary about the vector are the arrays     C
C MHM [ mhm( ih ) = m, the wavenumber for cos m phi and mhm(ih) = -m C
C for sin m phi dependence] and INARR (see INDFUN) which specifies   C
C the format.                                                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Int. parameter array corresponding to SV.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NR     See INDFUN for details         C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     MHT       : Harmonic type (poloidal velocity etc.)             C
C     MHL       : Spherical harmonic degree, l.                      C
C                                                                    C
C     MHM       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     TAU       : Longitudinal shift angle.                          C
C                                                                    C
C     SV        : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NR*NH                       C
C                                                                    C
C     SSV       : Shifted solution vector. Dim ( * ) ( input )       C
C                 Length must be atleast NR*NH                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVLSS( INARR, MHT, MHL, MHM, TAU, SV, SSV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * )
      DOUBLE PRECISION TAU, SV( * ), SSV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IPIV( 2 ), ITWO, INFO, IONE, NR, NH, IR, IHC, IHS,
     1        INDC, INDS, M, INDFUN, IH
      DOUBLE PRECISION DM, DMTAU, DSMT, DCMT,
     1                 RHS( 2, 1 ), DMAT( 2, 2 )
      CHARACTER *(1) TRANS
      PARAMETER ( IONE = 1, ITWO = 2, TRANS = 'N' )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      NR = INARR( 2 )
      NH = INARR( 3 )
C     .
      DO IHC = 1, NH
        M = MHM( IHC )
C       .
C       . First case of axisymmetric harmonics
C       .
        IF ( M.EQ.0 ) THEN
          DO IR = 1, NR
            INDC = INDFUN( IR, IHC, INARR )
            SSV( INDC ) = SV( INDC )
          ENDDO
        ENDIF
C       .
C       . Case of cos ( m phi ) harmonic with m non-zero
C       .
        IF ( M.GT.0 ) THEN
C          .
C          . First find corresponding sin ( m phi ) harmonic
C          .
           IHS = 0
           DO IH = 1, NH
             IF (  MHT( IH ).EQ.MHT( IHC )      .AND.
     1             MHL( IH ).EQ.MHL( IHC )      .AND.
     2           ( MHM( IH ) + MHM( IHC ) ).EQ.0      ) THEN
               IHS = IH
             ENDIF
           ENDDO
           IF ( IHS.EQ.0 ) THEN
             PRINT *,' Subroutine SVLSS.'
             PRINT *,' Harmonic with type = ', MHT( IHC )
             PRINT *,' l= ',MHL( IHC ),' and m = ',MHM( IHC )
             PRINT *,' does not appear to have a sin harmonic.'
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . OK we have both IHC and IHS.
C          . Now calculate transformation matrix for
C          . coefficients.
C          .
           DM = DBLE( M )
           DMTAU = DM*TAU
           DSMT  = DSIN( DMTAU )
           DCMT  = DCOS( DMTAU )
C          .
           DMAT( 1, 1 ) = DCMT
           DMAT( 1, 2 ) = DSMT
           DMAT( 2, 1 ) = (-1.0d0)*DSMT
           DMAT( 2, 2 ) = DCMT
C          .
C          . Perform LU decomposition on DMAT
C          .
           CALL DGETRF( ITWO, ITWO, DMAT, ITWO, IPIV, INFO )
C          .
           IF ( INFO.NE.0 ) THEN
             PRINT *,' Subroutine SVLSS.'
             PRINT *,' LAPACK subroutine DGETRF has returned'
             PRINT *,' the value INFO = ', INFO
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . OK - now loop around the radial grid nodes.
C          .
           DO IR = 1, NR
             INDC        = INDFUN( IR, IHC, INARR )
             INDS        = INDFUN( IR, IHS, INARR )
             RHS( 1, 1 ) = SV( INDC )
             RHS( 2, 1 ) = SV( INDS )
C            .
C            . Solve linear system for transformed co-ordinates
C            .
             CALL DGETRS( TRANS, ITWO, IONE, DMAT, ITWO, IPIV, RHS,
     1             ITWO, INFO )
C            .
             IF ( INFO.NE.0 ) THEN
               PRINT *,' Subroutine SVLSS.'
               PRINT *,' LAPACK subroutine DGETRS has returned'
               PRINT *,' the value INFO = ', INFO
               PRINT *,' Program aborted.'
               STOP
             ENDIF
C            .
             SSV( INDC ) = RHS( 1, 1 )
             SSV( INDS ) = RHS( 2, 1 )
C            .
           ENDDO
C          .
        ENDIF
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
