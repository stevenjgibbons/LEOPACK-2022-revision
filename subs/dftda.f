C*********************************************************************
C subroutine Drifting Frame Time Derivative Add **********************
C            -        -     -    -          -   **********************
C Steve Gibbons Wed Feb 23 09:13:59 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C If Y_a^c and Y_a^s are defined by                                  C
C                                                                    C
C  Y_a^c = P_{la}^{mc}(cos theta) cos {ma} phi     and               C
C  Y_a^s = P_{la}^{mc}(cos theta) sin {ma} phi         then          C
C                                                                    C
C in a rotating frame with phi' = phi - ct,                          C
C                                                                    C
C  d Y_a^c/dt = c {ma} Y_a^s    and                                  C
C  d Y_a^s/dt = (-c) {ma} Y_a^c                                      C
C                                                                    C
C If the time derivatives in the heat, vorticity and induction       C
C equations are respectively  CA, CE and CK, then DFTDA adds a       C
C multiple of FAC*this time derivative to the forcing term.          C
C                                                                    C
C Zeroth, first and second derivatives have been calculated in       C
C advance by CASVDR and are stored in VI0, VI1 and VI2.              C
C                                                                    C
C The drift velocity, c, is assumed to be 1.0. Departures            C
C from this are ofcourse taken into account in FAC.                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Int. parameter array corresponding to VI.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C 
C                 INARR( 2 ) = NR      See INDFUN for details        C
C                 INARR( 3 ) = NH                                    C
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
C     MHTR      : Array length ( * ) - atleast length NH             C
C                 Equivalent to MHT but for RHS (see CINDSW).        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VI0       : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NR*NH                       C
C     VI1       : First derivs of VI0. Dim ( * ) ( input )           C
C     VI2       : Second derivs of VI0. Dim ( * ) ( input )          C
C     RHS       : Output vector. Dim ( * ).                          C
C                                                                    C
C     FAC       : Linear scaling of term.                            C
C     CA        : Coeff. of time derivative in heat equation.        C
C     CE        : Coeff. of time derivative in vortcity equation.    C
C     CK        : Coeff. of time derivative in induction equation.   C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DFTDA( INARR, MHT, MHL, MHM, MHTR, VI0, VI1, VI2,
     1                  RHS, FAC, CA, CE, CK, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * ), MHTR( * )
      DOUBLE PRECISION VI0( * ), VI1( * ), VI2( * ), RHS( * ),
     1                 FAC, CA, CE, CK, XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NH, NR, IHI, IHO, IR, INDO, INDI, INDFUN, L, ITYPEI,
     1        ITYPEO, ILNR, IRNR
      DOUBLE PRECISION LOW, RAD, DL, D0F, D1F, D2F, COEF, TERM, DM
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     . early exit?
C
      IF ( ABS(FAC).LT.LOW ) RETURN
C     .
      NR = INARR( 2 )
      NH = INARR( 3 )
C     .
C     . Loop around the input harmonics
C     .
      DO IHI = 1, NH
        IF ( MHT( IHI ).LT.1 .OR. MHT( IHI ).GT.5 ) GOTO 60
        IF ( MHM( IHI ).EQ.0 ) GOTO 60
        L  = MHL( IHI )
        DM = DBLE(  MHM( IHI )  )
C       .
C       . Ok - this is a non-axisymmetric harmonic
C       .
        ITYPEI = MHT( IHI )
        ILNR   = 2
        IRNR   = NR - 1
        IF ( ITYPEI.EQ.1 ) THEN
          ILNR   = 3
          IRNR   = NR - 2
        ENDIF
C       .
        COEF   = CA
        IF ( ITYPEI.EQ.1 .OR. ITYPEI.EQ.2 ) COEF = CE
        IF ( ITYPEI.EQ.4 .OR. ITYPEI.EQ.5 ) COEF = CK
C       .
        ITYPEO = ITYPEI
        IF ( ITYPEI.EQ.1 ) ITYPEO = 2
        IF ( ITYPEI.EQ.2 ) ITYPEO = 1
C       .
C       . Loop around output harmonics ...
C       .
        DO IHO = 1, NH
          IF (      MHTR( IHO ).EQ.ITYPEO        .AND.
     1              MHL( IHO ).EQ.L              .AND.
     2              (MHM( IHO )+MHM( IHI )).EQ.0     ) THEN
C           .
C           . We have found the appropriate harmonic
C           .
            DO IR = ILNR, IRNR
              INDI = INDFUN( IR, IHI, INARR )
              INDO = INDFUN( IR, IHO, INARR )
              D0F  = VI0( INDI )
              IF ( ITYPEI.EQ.1 ) THEN
                RAD  = XARR( IR )
                D1F  = VI1( INDI )
                D2F  = VI2( INDI )
                TERM = (-1.0d0)*DL( L, RAD, D0F, D1F, D2F )
              ELSE
                TERM = D0F
              ENDIF
              RHS( INDO ) = RHS( INDO ) + FAC*DM*COEF*TERM
            ENDDO
C           .
          ENDIF
        ENDDO
C       .
 60   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
