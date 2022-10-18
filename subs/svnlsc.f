C*********************************************************************
C subroutine Solution Vector Necessary Longitudinal Shift Calculate **
C            -        -      -         -            -     -         **
C Steve Gibbons Mon Mar 20 16:02:48 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Let SV be a solution vector containing radial functions of the     C
C standard spherical harmonics with respect to longitude phi         C
C                                                                    C
C The user supplies a type (1,2,3,4 or 5 - usual system),            C
C a degree, l, a wavenumber, m, and a grid node number IGN  -        C
C SVNLSC then calculates how big a shift in latitude must be applied C
C (TAU) such that the cos ( m phi ) harmonic vanishes at this        C
C grid node. It also returns as ITCDC, the number of the cos (m phi) C
C harmonic.                                                          C
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
C     IT        : Type: 1, 2, 3, 4 or 5.                             C
C     IL        : Degree, l.                                         C
C     IM        : Order, m.                                          C
C     IGN       : Grid node.                                         C
C                                                                    C
C     ITCDC     : Selected cos (M phi) harmonic number.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     TAU       : Longitudinal shift angle.  (Output)                C
C                                                                    C
C     SV        : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NR*NH                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVNLSC( INARR, MHT, MHL, MHM, IT, IL, IM, IGN,
     1                   ITCDC, TAU, SV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * ), IT, IL, IM,
     1        IGN, ITCDC
      DOUBLE PRECISION TAU, SV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NR, NH, IHC, IHS,
     1        INDC, INDS, INDFUN, IH
      DOUBLE PRECISION FMS, FMC, GMS
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      OK = .FALSE.
C     .
      NR = INARR( 2 )
      NH = INARR( 3 )
C     .
      IF ( IGN.LT.1 .OR. IGN.GT.NR ) THEN
        PRINT *,' Subroutine SVNLSC.'
        PRINT *,' Suggested grid node = ', IGN
        PRINT *,' Number of nodes = ', NR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      DO IHC = 1, NH
        IF ( MHT( IHC ).EQ.IT .AND. MHL( IHC ).EQ.IL .AND.
     1       MHM( IHC ).EQ.IABS( IM )  ) THEN
C          .
           OK    = .TRUE.
           ITCDC = IHC
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
             PRINT *,' Subroutine SVNLSC.'
             PRINT *,' Harmonic with type = ', MHT( IHC )
             PRINT *,' l= ',MHL( IHC ),' and m = ',MHM( IHC )
             PRINT *,' does not appear to have a sin harmonic.'
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . OK we have both IHC and IHS.
C          .
           INDC    = INDFUN( IGN, IHC, INARR )
           INDS    = INDFUN( IGN, IHS, INARR )
           FMC     = SV( INDC )
           FMS     = SV( INDS )
C          .
           CALL PF2SOF( IM, FMC, FMS, GMS, TAU )
C          .
        ENDIF
C       .
      ENDDO
C     .
      IF ( .NOT. OK ) THEN
        PRINT *,' Subroutine SVNLSC.'
        PRINT *,' Failiure to find requested harmonic.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      RETURN
      END
C*********************************************************************
