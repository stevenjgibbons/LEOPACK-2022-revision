C*********************************************************************
C subroutine Single Harmonic VECtor Fill *****************************
C            -      -        ---    -    *****************************
C Steve Gibbons Thu Nov 18 10:41:52 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Takes a harmonic number IH and fills that harmonic with a          C
C (very arbitrarily generated!) function which satisfies the         C
C appropriate boundary conditions.                                   C
C                                                                    C
C  The arrays MHIBC and MHOBC instruct SHVECF how to manipulate      C
C  the finite difference coefficients at the boundaries.             C
C                                                                    C
C  MHIBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( ih ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                        where L = MHL( ih )                         C
C                                                                    C
C  Similarly, at the outer boundary:-                                C
C                                                                    C
C  MHOBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( ih ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                        where L = MHL( ih )                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IH        : Number of harmonic.                                C
C                                                                    C
C     INARR     : Int. parameter array corresponding to V.           C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR     See INDFUN for details        C
C                 INARR( 3 ) = NH      nrr must equal nr             C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C     MHT( ih ) = 1 if 'ih' is poloidal velocity harmonic.           C
C     MHT( ih ) = 2 if 'ih' is toroidal velocity harmonic.           C
C     MHT( ih ) = 3 if 'ih' is temperature harmonic.                 C
C     MHT( ih ) = 4 if 'ih' is poloidal magnetic field harmonic.     C
C     MHT( ih ) = 5 if 'ih' is toroidal magnetic field harmonic.     C
C                                                                    C
C     MHL       : Spherical harmonic degree, l, of harmonic 'ih'.    C
C                                                                    C
C     MHM       : Spherical harmonic degree, m, of harmonic 'ih'     C
C                if harmonic has (cos m phi dependence) - otherwise  C
C                MHM( ih ) = -m                                      C
C     MHP       : Pointer array to finite diff scheme.               C
C                                                                    C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). See above.                     C
C     MHOBC     : Dimension ( NDCS ). See above.                     C
C                                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     VEC       : Solution vector - dim ( NR*NH )                    C
C                                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHVECF( IH, INARR, MHT, MHL, MHM, MHP, NDCS,
     1                   MHIBC, MHOBC, XARR, VEC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IH, INARR( 3 ), MHT( * ), MHL( * ), MHM( * ), MHP( * ),
     1        NDCS, MHIBC( NDCS ), MHOBC( NDCS )
      DOUBLE PRECISION XARR( * ), VEC( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NR, NH, IR, L, IS, IIBCF, IOBCF, M, ICS, IT, MAXM,
     1        INTPA( 1 ), IND, INDFUN
      DOUBLE PRECISION RAD, COEFS( 4 ), RI, RO, DPRPA( 1 ),
     1        PVNSRF, TVSFRF, PVSFRF, PMFIRF, FAC, PI, F1
      LOGICAL OK
      PARAMETER ( PI=3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      INTPA( 1 ) = 0
      DPRPA( 1 ) = 0.0d0
C
      NR = INARR( 2 )
      NH = INARR( 3 )
      OK = .FALSE.
C  
      RI = XARR( 1  )
      RO = XARR( NR )
C  
      IF ( IH.LT.1 .OR. IH.GT.NH ) THEN
        PRINT *,' Subroutine SHVECF.'
        PRINT *,' IH = ', IH,' NH = ', NH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IT      = MHT( IH )
      IS      = MHP( IH )
      IIBCF   = MHIBC( IS )
      IOBCF   = MHOBC( IS )
C     .
      L = MHL( IH )
C     .
      IF ( MHM( IH ).LT.0 ) THEN
        M   = -MHM( IH )
        ICS = 2
      ELSE
        M   = MHM( IH )
        ICS = 1
      ENDIF
C     .
      IF ( IIBCF.EQ.1 .AND. IOBCF.EQ.1 ) OK = .TRUE.
      IF ( IT.EQ.1 .AND. IIBCF.EQ.4 .AND. IOBCF.EQ.4 ) OK = .TRUE.
      IF ( IT.EQ.1 .AND. IIBCF.EQ.5 .AND. IOBCF.EQ.5 ) OK = .TRUE.
      IF ( IT.EQ.2 .AND. IIBCF.EQ.2 .AND. IOBCF.EQ.2 ) OK = .TRUE.
      IF ( IT.EQ.2 .AND. IIBCF.EQ.6 .AND. IOBCF.EQ.6 ) OK = .TRUE.
      IF ( IT.EQ.3 .AND. IIBCF.EQ.2 .AND. IOBCF.EQ.2 ) OK = .TRUE.
      IF ( IT.EQ.3 .AND. IIBCF.EQ.2 .AND. IOBCF.EQ.3 ) OK = .TRUE.
      IF ( IT.EQ.3 .AND. IIBCF.EQ.3 .AND. IOBCF.EQ.2 ) OK = .TRUE.
      IF ( IT.EQ.4 .AND. IIBCF.EQ.7 .AND. IOBCF.EQ.7 .AND.
     1     L.GE.1 ) OK = .TRUE.
      IF ( IT.EQ.5 .AND. IIBCF.EQ.2 .AND. IOBCF.EQ.2 ) OK = .TRUE.
C     .
      IF ( .NOT. OK ) THEN
        PRINT *,' Subroutine SHVECF.'
        PRINT *,' IH = ', IH,' L = ', L,' M = ', M,' ICS = ', ICS
        PRINT *,' IIBCF = ', IIBCF,' IOBCF = ', IOBCF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . 'Random' function ...
C     .
      IF ( IIBCF.EQ.1 .AND. IOBCF.EQ.1 ) THEN
        DO IR = 1, NR
          IND = INDFUN( IR, IH, INARR )
          RAD = XARR( IR )
          VEC( IND ) = RAD*RAD - RAD
        ENDDO
      ENDIF
C     .
C     . Poloidal velocity - no slip boundaries
C     .
      IF ( IT.EQ.1 .AND. IIBCF.EQ.4 ) THEN
C       .
        MAXM = 2
        COEFS( 1 ) = 1.5d0
        COEFS( 2 ) = -0.8d0
        COEFS( 3 ) = -0.3d0
        COEFS( 4 ) = 0.2d0
C       .
        DO IR = 1, NR
          IND = INDFUN( IR, IH, INARR )
          RAD = XARR( IR )
          VEC( IND ) = PVNSRF( RAD, RI, RO, MAXM, COEFS,
     1                         INTPA, DPRPA )
        ENDDO
      ENDIF
C     .
C     . Poloidal velocity - stress free boundaries
C     . or toroidal velocity with no slip boundaries
C     . or toroidal magnetic field with insulating bdries
C     . or temperature with fixed tm inner and outer
C     .
      IF ( ( IT.EQ.1 .AND. IIBCF.EQ.5 ) .OR.
     1     ( IT.EQ.2 .AND. IIBCF.EQ.2 ) .OR.
     2     ( IT.EQ.5 .AND. IIBCF.EQ.7 ) .OR.
     3     ( IT.EQ.3 .AND. IIBCF.EQ.2 .AND. IOBCF.EQ.2 ) ) THEN
C       .
        MAXM = 4
        COEFS( 1 ) = 1.5d0
        COEFS( 2 ) = -0.8d0
        COEFS( 3 ) = -0.3d0
        COEFS( 4 ) = 0.2d0
C       .
        DO IR = 1, NR
          IND = INDFUN( IR, IH, INARR )
          RAD = XARR( IR )
          VEC( IND ) = PVSFRF( RAD, RI, RO, MAXM, COEFS,
     1                         INTPA, DPRPA )
        ENDDO
      ENDIF
C     .
C     . Toroidal velocity - stress free boundaries
C     .
      IF ( IT.EQ.2 .AND. IIBCF.EQ.6 ) THEN
C       .
        MAXM = 3
        COEFS( 1 ) = 0.0d0
        COEFS( 2 ) = 0.8d0
        COEFS( 3 ) = -0.3d0
        COEFS( 4 ) = 0.2d0
C       .
        DO IR = 1, NR
          IND = INDFUN( IR, IH, INARR )
          RAD = XARR( IR )
          VEC( IND ) = TVSFRF( RAD, RI, RO, MAXM, COEFS,
     1                         INTPA, DPRPA )
        ENDDO
      ENDIF
C     .
C     . Poloidal magnetic field - insulating boundaries
C     .
      IF ( IT.EQ.4 .AND. IIBCF.EQ.7 ) THEN
C       .
        INTPA( 1 ) = L
        MAXM = 4
        COEFS( 1 ) = 1.5d0
        COEFS( 2 ) = -0.8d0
        COEFS( 3 ) = -0.3d0
        COEFS( 4 ) = 0.2d0
C       .
        DO IR = 1, NR
          IND = INDFUN( IR, IH, INARR )
          RAD = XARR( IR )
          VEC( IND ) = PMFIRF( RAD, RI, RO, MAXM, COEFS,
     1                         INTPA, DPRPA )
        ENDDO
      ENDIF
C     .
C     . Temperature with opposite b.c.s
C     .
      IF ( IT.EQ.3 .AND. IIBCF.NE.IOBCF ) THEN
        IF ( IIBCF.EQ.2 ) FAC = 0.0d0
        IF ( IIBCF.EQ.3 ) FAC = 0.5d0*PI
        F1 = 0.5d0*PI/(RO-RI)
        DO IR = 1, NR
          IND = INDFUN( IR, IH, INARR )
          RAD = XARR( IR )
          VEC( IND ) = DSIN( F1*(RAD - RI) + FAC )
        ENDDO
      ENDIF
C     .
      RETURN
      END
C*********************************************************************
