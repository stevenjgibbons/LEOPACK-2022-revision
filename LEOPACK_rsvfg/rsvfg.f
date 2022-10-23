C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C Mon Mar  6 10:14:00 GMT 2000                                       C
C                                                                    C
C Random Solution Vector File Generator.                             C
C                                                                    C
C Reads in from standard input a filename stem, radial spacing info, C
C symmetry info etc. and boundary conditions and fills (using the    C
C routine SHVECF) a very arbitrary non-zero function into the        C
C solution vector.                                                   C
C                                                                    C
C ROOT.ints, ROOT.vecs and ROOT.xarr are output is normal.           C
C                                                                    C
C*********************************************************************
      PROGRAM rsvfg
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, NHMAX, LHMAX, NDCS, ISVMAX
      PARAMETER ( NRMAX = 300, NHMAX = 6000, LHMAX = 64,
     1            NDCS  = LHMAX + 6, ISVMAX = NRMAX*NHMAX )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ), MHP( NHMAX ),
     1        MHIBC( NDCS ), MHOBC( NDCS ), ISYMA( 5 ), 
     2        LHARR( 5 ), MMODES( LHMAX + 1 )
      DOUBLE PRECISION XARR( NRMAX ), VEC( ISVMAX )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      CHARACTER *(80) LINE, ROOT, FNAME
      INTEGER I, ILEN, NR, INSF, IFORMF, INARR( 3 ), IH, LU,
     1        ISM, ISV, LHM, LHV, MINC, MMAX, M, NH, NMODES,
     2        IVELBC, ITHEBC, IFORM
      DOUBLE PRECISION RI, RO, LOW, ZERO, FAC
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( LOW = 1.0d-7, ZERO = 0.0d0 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFORM = 1
      LU    = 91
C
C Fill in integer arrays for boundary conditions
C
      DO I = 1, LHMAX
        MHIBC( I ) = 7
        MHOBC( I ) = 7
      ENDDO
C ----------------------------- f(r) = 0 at both boundaries
      MHIBC( LHMAX + 1 ) = 2
      MHOBC( LHMAX + 1 ) = 2
C ----------------------------- toroidal stress free velocity
      MHIBC( LHMAX + 2 ) = 6
      MHOBC( LHMAX + 2 ) = 6
C ----------------------------- f(r) = 0 at both boundaries
      MHIBC( LHMAX + 3 ) = 4
      MHOBC( LHMAX + 3 ) = 4
C ----------------------------- f(r) = 0 at both boundaries
      MHIBC( LHMAX + 4 ) = 5
      MHOBC( LHMAX + 4 ) = 5
C ----------------------------- f(ri) = 0, f'(ro) = 0
      MHIBC( LHMAX + 5 ) = 2
      MHOBC( LHMAX + 5 ) = 3
C ----------------------------- f'(ri) = 0, f(ro) = 0
      MHIBC( LHMAX + 6 ) = 3
      MHOBC( LHMAX + 6 ) = 2
C----------------------------------------------------------------
C                Enter stem of filename for output vectors      |
C----------------------------------------------------------------
      PRINT *,' Enter name of filename stem.'
 80   FORMAT(A)
 31   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 31
C
      DO I = 1, 80
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 45
        ENDIF
      ENDDO
 45   CONTINUE
      ROOT  = LINE(1:ILEN)
C
      PRINT *,' ----------------------------------------'
      PRINT *,' Enter NR, IFORMF.'
      PRINT *,' nrmax = ', NRMAX
      PRINT *,' iformf = 3 --> ind = ( ir - 1 )*nh + ih '
      PRINT *,' iformf = 4 --> ind = ( ih - 1 )*nr + ir '
      PRINT *,' ----------------------------------------'
 32   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 32
      READ ( LINE, * ) NR, IFORMF
      PRINT *,' NR = ', NR,' IFORMF = ', IFORMF
C
C Check values of NR, IFORMF
C
      IF ( NR.LT.6 .OR. NR.GT.NRMAX ) THEN
        PRINT *,' NR = ', NR,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
        PRINT *,' IFORMF = ', IFORMF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INARR( 1 ) = IFORMF
      INARR( 2 ) = NR
C
C Enter inner and outer radii
C
      PRINT *,' -------------------------------------------'
      PRINT *,' Enter RI, RO, INSF, FAC '
      PRINT *,' insf = 1 --> equally spaced nodes '
      PRINT *,' insf = 2 --> Chebyshev zero spaced nodes. '
      PRINT *,' fac simply scales whole solution vector '
      PRINT *,' -------------------------------------------'
 33   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 33
      READ ( LINE, * ) RI, RO, INSF, FAC
      PRINT *,' RI = ', RI,' RO = ', RO,' INSF = ',INSF
C
C Check values of RI, RO, INSF
C
      IF ( RI.LT.ZERO ) THEN
        PRINT *,' RI = ', RI
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( (RO-RI).LT.LOW ) THEN
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( INSF.NE.1 .AND. INSF.NE.2 ) THEN
        PRINT *,' INSF = ', INSF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Calculate grid nodes
C
      IF ( INSF.EQ.1 ) THEN
        CALL ESNAAS( NR, XARR, RI, RO)
      ELSE
        CALL ZCPAAS( NR, XARR, RI, RO)
      ENDIF
C
C----------------------------------------------------------------
C                Enter dimensions of the new vector.            |
C----------------------------------------------------------------
C
      PRINT *,' Enter LHV, ISV, LHM, ISM, MINC, MMAX '
      PRINT *,' LHV  = maximum degree, l for velocity/temp. '
      PRINT *,' ISV  = 1 --> velocity is eq. symmetric.'
      PRINT *,' ISV  = 2 --> velocity is eq. anti-symmetric.'
      PRINT *,' ISV  = 3 --> velocity has both symmetries.'
      PRINT *,' LHM  = maximum degree, l for magnetic field.'
      PRINT *,' ISM  = 1 --> mag. field is eq. symmetric.'
      PRINT *,' ISM  = 2 --> mag. field is eq. anti-symmetric.'
      PRINT *,' ISM  = 3 --> mag. field has both symmetries.'
      PRINT *,' MINC is increment of wavenumber.'
      PRINT *,' MMAX is maximum desired wavenumber.'
C
 34   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 34
      READ ( LINE, * ) LHV, ISV, LHM, ISM, MINC, MMAX
C
      IF ( LHV.GT.LHMAX ) THEN
        PRINT *,' LHV = ', LHV,' LHMAX = ', LHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LHM.GT.LHMAX ) THEN
        PRINT *,' LHM = ', LHM,' LHMAX = ', LHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      M = 0
      NMODES = 1
      MMODES( 1 ) = 0
 40   CONTINUE
      M = M + 1
      IF ( M.GT.MMAX ) GOTO 50
      IF ( M/MINC*MINC.NE.M ) GOTO 40
      NMODES = NMODES + 1
      MMODES( NMODES ) = M
      GOTO 40
 50   CONTINUE
C
      LHARR( 1 ) = LHV
      LHARR( 2 ) = LHV
      LHARR( 3 ) = LHV
      LHARR( 4 ) = LHM
      LHARR( 5 ) = LHM
C
      ISYMA( 1 ) = ISV
      ISYMA( 2 ) = ISV
      ISYMA( 3 ) = ISV
      ISYMA( 4 ) = ISM
      ISYMA( 5 ) = ISM
C
      CALL HMINDA( LHARR, ISYMA, NMODES, MMODES, NH, NHMAX,
     1             MHT, MHL, MHM, LHMAX )
      PRINT *,' Total harmonics selected = ', NH
      IF ( NH.EQ.0 ) THEN
        PRINT *,' Program aborted. NH = ', NH
        STOP
      ENDIF
C
      PRINT *,' Enter IVELBC, ITHEBC '
      PRINT *,' (ivelbc = 1 --> no slip )'
      PRINT *,' (ivelbc = 2 --> stress free)'
      PRINT *,' (ithebc = 1 --> fixed tm inner and outer)'
      PRINT *,' (ithebc = 2 --> fixed tm inner hf outer)'
      PRINT *,' (ithebc = 3 --> fixed hf inner tm outer)'
C
 35   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 35
      READ ( LINE, * ) IVELBC, ITHEBC
C
      IF ( IVELBC.NE.1 .AND. IVELBC.NE.2 ) THEN
        PRINT *,' IVELBC = ', IVELBC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ITHEBC.NE.1 .AND. ITHEBC.NE.2 .AND. ITHEBC.NE.3 ) THEN
        PRINT *,' ITHEBC = ', ITHEBC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INARR( 3 ) = NH
      DO IH = 1, NH
        IF ( MHT( IH ).EQ.1 ) THEN
          IF ( IVELBC.EQ.1 ) MHP( IH ) = LHMAX + 3
          IF ( IVELBC.EQ.2 ) MHP( IH ) = LHMAX + 4
        ENDIF
        IF ( MHT( IH ).EQ.2 ) THEN
          IF ( IVELBC.EQ.1 ) MHP( IH ) = LHMAX + 1
          IF ( IVELBC.EQ.2 ) MHP( IH ) = LHMAX + 2
        ENDIF
        IF ( MHT( IH ).EQ.3 ) THEN
          IF ( ITHEBC.EQ.1 ) MHP( IH ) = LHMAX + 1
          IF ( ITHEBC.EQ.2 ) MHP( IH ) = LHMAX + 5
          IF ( ITHEBC.EQ.3 ) MHP( IH ) = LHMAX + 6
        ENDIF
        IF ( MHT( IH ).EQ.4 ) THEN
          MHP( IH ) = MHL( IH )
        ENDIF
        IF ( MHT( IH ).EQ.5 ) THEN
          MHP( IH ) = LHMAX + 1
        ENDIF
        CALL SHVECF( IH, INARR, MHT, MHL, MHM, MHP, NDCS,
     1               MHIBC, MHOBC, XARR, VEC )
      ENDDO
C
C Scale vector by FAC
C
      DO I = 1, NR*NH
        VEC( I ) = VEC( I )*FAC
      ENDDO
C
C Finally: write out all of the files
C
      FNAME(1:ILEN) = ROOT(1:ILEN)
      FNAME(ILEN+1:ILEN+5) = '.ints'
      FNAME = FNAME(1:ILEN+5)
C
      CALL HMFWT( NH, MHT, MHL, MHM, MHP, NDCS,
     1            MHIBC, MHOBC, LU, FNAME )
C
      FNAME(ILEN+1:ILEN+5) = '.xarr'
      FNAME = FNAME(1:ILEN+5)
      CALL XARRWT( NR, XARR, LU, FNAME, IFORM )
C
      FNAME(ILEN+1:ILEN+5) = '.vecs'
      FNAME = FNAME(1:ILEN+5)
      CALL SVFWT( INARR, LU, IFORM, VEC, FNAME )
C
      END
C*********************************************************************
C*********************************************************************
C subroutine Equally Spaced Node Abscissa Allocation Subroutine. *****
C            -       -      -    -        -          -           *****
C Steve Gibbons Thu Oct 21 11:35:00 BST 1999                         C
C                                                                    C
C This short routine simply fills the array XARR with xvalues        C
C which are evenly spaced with                                       C
C                                                                    C
C  x_j = r_i + ( j - 1 )*h with h = ( r_o - r_i )/( nr - 1 )         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     RI        : Radius of inner boundary                           C
C     RO        : Radius of outer boundary                           C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ESNAAS( NR, XARR, RI, RO)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR
      DOUBLE PRECISION RI, RO, XARR( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR
      DOUBLE PRECISION H, TOL
      PARAMETER ( TOL = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NR.LT.4 ) THEN
         PRINT *,' Subroutine ESNAAS.'
         PRINT *,' NR = ', NR
         STOP
      ENDIF
C
      IF ( (RO-RI).LT.TOL ) THEN
         PRINT *,' Subroutine ESNAAS.'
         PRINT *,' RI = ', RI
         PRINT *,' RO = ', RO
         STOP
      ENDIF
C
      H = ( RO - RI )/DBLE( NR - 1 )
      DO IR = 1, NR
        XARR( IR ) = RI + DBLE( IR - 1 )*H
      ENDDO
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Zeroes of Chebyshev Polynomial Abscissa Allocation Sub. *
C            -         -         -          -        -          -    *
C Steve Gibbons Tue Sep 21 17:46:55 BST 1999                         C
C This routine is almost entirely derived from ZECHGA written        C
C by D. Funaro of Consiglio Nazionale Delle Ricerche,                C
C Via Abbiategrasso, 209 - 27100 Pavia, Italy.                       C
C However, rather than giving the N zeros of the Chebyshev polynom.  C
C of degree N in the interval (-1,1), it returns XARR( 1 ) = RI,     C
C XARR( NR ) = RO and XARR( i + 1 ) as the i^{th} zero of the        C
C Chebyshev polynomial of degree (NR-2) scaled into the interval     C
C (RI,RO) ie. x := (ro-ri)(x+1.0d0)/2.0d0 + ri                       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     RI        : Radius of inner boundary                           C
C     RO        : Radius of outer boundary                           C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ZCPAAS( NR, XARR, RI, RO)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR
      DOUBLE PRECISION RI, RO, XARR( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER N, N2, I, IN
      DOUBLE PRECISION PH, DN, C, SI, DI, CSX, RIPRO2,
     1                 RIMRO2
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NR.LT.4 ) THEN
         PRINT *,' Subroutine ZCPAAS.'
         PRINT *,' NR = ', NR
         STOP
      ENDIF
C
      IF ( RI.GE.RO ) THEN
         PRINT *,' Subroutine ZCPAAS.'
         PRINT *,' RI = ', RI
         PRINT *,' RO = ', RO
         STOP
      ENDIF
C
      RIPRO2 = 0.5d0*(RO + RI)
      RIMRO2 = 0.5d0*(RO - RI)
C
      XARR(  1 ) = RI
      XARR( NR ) = RO
      N = NR - 2
      XARR(2) = 0.D0
      N2 = N/2
      IN = 1+4*N2-2*N
      PH = 1.57079632679489661923D0
      DN = DFLOAT(N)
      C  = PH/DN
      SI = -1.D0
      DO 10 I = 1, N2
         DI = DFLOAT(I)
         CSX = DCOS(C*(2.D0*DI-1.D0))
         XARR(I+1) = (1.0d0-CSX)*RIMRO2 + RI
         XARR(N-I+2) = (1.0d0+CSX)*RIMRO2 + RI
         SI = -SI
 10   CONTINUE
C
      IF (IN .EQ. 1) RETURN
      XARR(N2+2) = RIPRO2
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine HarMonic INDex Allocation subroutine ********************
C            -  -     ---   -                     ********************
C Steve Gibbons Sat Sep 25 11:52:41 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Fills the arrays MHT, MHL and MHM according to user-specified      C
C conditions on symmetry (equatorial and rotational) dimension and   C
C choice of variables.                                               C
C                                                                    C
C HMINDA will allocate indices to these arrays in ascending order    C
C of modes with symmetries separate, regardless of type of function. C
C                                                                    C
C The maximum degree of spherical harmonic, LH, can be set           C
C differently for each kind of function and they can have different  C
C specified equatorial symmetries. HMINDA does not check on whether  C
C the supplied combination is a physically sensible or meaningful    C
C choice.                                                            C
C                                                                    C
C  For example, set LH( 1 ) = 8      (poloidal velocity)             C
C                   LH( 2 ) = 8      (toroidal velocity)             C
C                   LH( 3 ) = 8      (temperature)                   C
C                   LH( 4 ) = 0      (poloidal magnetic field)       C
C                   LH( 5 ) = 0      (toroidal magnetic field)       C
C                                                                    C
C  Set NMODES = 2 and MMODES( 1 ) = 0                                C
C                     MMODES( 2 ) = 4                                C
C                                                                    C
C  Set ISYM( 1 ) = 3  This will select both equatorially symmetric   C
C      ISYM( 2 ) = 3  and equatorially anti-symmetric                C
C      ISYM( 3 ) = 3    harmonics for first 3 components             C
C      ISYM( 4 ) = 0    There are not to be any magnetic field       C
C      ISYM( 5 ) = 0      harmonics.                                 C
C                                                                    C
C  HMINDA will loop around the wavenumbers (imode)                   C
C  We then loop around the equatorial symmetries first ES then EA    C
C  If either is not wanted, we move straight onto the next.          C
C  It will then loop around ITYPE ( = 1, 2, 3, 4 , 5)                C
C  If LH( itype ).lt.MMODES( imode ) then we move onto the           C
C  next type.                                                        C
C  Finally, if MMODES( imode ).ne.0, we loop around ICS = cos        C
C                                                     to ICS = sin   C
C                                                                    C
C  In the following, ITYPE will always refer to                      C
C                                                                    C
C         ITYPE = 1 for a poloidal velocity harmonic.                C
C         ITYPE = 2 for a toroidal velocity harmonic.                C
C         ITYPE = 3 for a temperature harmonic.                      C
C         ITYPE = 4 for a poloidal magnetic field harmonic.          C
C         ITYPE = 5 for a toroidal magnetic field harmonic.          C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LHARR     : Dim ( 5 ). LHARR( itype ) contains maximum degree  C
C                 of spherical harmonic for type ITYPE.              C
C                  LHARR( itype ) not referred to if                 C
C                   ISYMA( itype ).eq.0.                             C
C                                                                    C
C     ISYMA     : Dim ( 5 ). Specifies equatorial symmetry property  C
C                 of type ITYPE.                                     C
C                                                                    C
C                ISYMA( itype ).eq.0 --> we are simply not including C
C                                        any harmonics of this type. C
C                ISYMA( itype ).eq.1 --> only equatorially symmetric C
C                ISYMA( itype ).eq.2 --> only equatorially anti-sym. C
C                ISYMA( itype ).eq.3 --> all eq. symmetries          C
C                                                                    C
C     NMODES    : Number of wavenumbers ( m ) to be included in the  C
C                   solution vector.                                 C
C                                                                    C
C     MMODES    : Integer array Dim.( NMODES )                       C
C                  NMODES( im ) contains m for this mode             C
C                                                                    C
C     NH        : Output only integer giving the total number of     C
C                  harmonics selected.                               C
C                                                                    C
C     NHMAX     : Maximum number of harmonics permitted. If the      C
C                  specified parameters demand more harmonics than   C
C                   permitted by this bound then HMINDA calculates   C
C                    how many harmonics are necessary and then       C
C                     aborts with an appropriate message.            C
C                                                                    C
C     MHT       : MHT( ih ) contains itype for harmonic 'ih'         C
C                                                                    C
C     MHL       : MHL( ih ) contains degree, l, for harmonic 'ih'    C
C                                                                    C
C     MHM       : MHM( ih ) contains order, m, for harmonic 'ih' if  C
C                  ih has cos m phi dependency and (-m) if ih has    C
C                   sin m phi dependency.                            C
C                                                                    C
C     LHMAX     : Global maximum permitted degree, l, of a spherical C
C                  harmonic.                                         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE HMINDA( LHARR, ISYMA, NMODES, MMODES, NH, NHMAX,
     1                   MHT, MHL, MHM, LHMAX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LHARR( 5 ), ISYMA( 5 ), NMODES, MMODES( NMODES ),
     1        NH, NHMAX, MHT( * ), MHL( * ), MHM( * ), LHMAX
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L, M, IMODE, LMIN, LMAX, ISYM, ICS, M2, ITYPE
      LOGICAL OES( 5 ), OEA( 5 ), OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      DO ITYPE = 1, 5
C
        OES( ITYPE ) = .FALSE.
        OEA( ITYPE ) = .FALSE.
        IF (           ISYMA( ITYPE ).EQ.1 .OR.
     1                 ISYMA( ITYPE ).EQ.3      )
     2                         OES( ITYPE ) = .TRUE.
        IF (           ISYMA( ITYPE ).EQ.2 .OR.
     1                 ISYMA( ITYPE ).EQ.3      )
     2                         OEA( ITYPE ) = .TRUE.
C
        IF ( LHARR( ITYPE ).GT.LHMAX ) THEN
          PRINT *,' Subroutine HMINDA.'
          PRINT *,' LHARR(',ITYPE,') = ',LHARR( ITYPE )
          PRINT *,' LHMAX = ', LHMAX
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
      ENDDO
c
      NH = 0
      OK = .TRUE.
c
      DO IMODE = 1, NMODES
        M = MMODES( IMODE )
C
C isym = 1 will consider the equatorially symmetric
C harmonics, isym = 2 the equatorially antisymmetric
C
        DO ISYM = 1, 2
          DO ITYPE = 1, 5
C
            IF ( ITYPE.EQ.1 ) M2 = 0
            IF ( ITYPE.EQ.2 ) M2 = 1
            IF ( ITYPE.EQ.3 ) M2 = 0
            IF ( ITYPE.EQ.4 ) M2 = 0
            IF ( ITYPE.EQ.5 ) M2 = 1
C
C If an ES function has (l-m) even then M2 is 0
C If an ES function has (l-m) odd then M2 is 1
C
            IF ( ISYMA( ITYPE ).EQ.0 ) GOTO 50
            IF ( .NOT. OES( ITYPE ) .AND. ISYM.EQ.1 ) GOTO 50
            IF ( .NOT. OEA( ITYPE ) .AND. ISYM.EQ.2 ) GOTO 50
C
C set LMIN and LMAX - first checking for the monopole term
C
            LMIN = MMODES( IMODE )
            IF (     LMIN.EQ.0 .AND. 
     1             ( ISYM.NE.1 .OR. ITYPE.NE.3 ) ) LMIN = 1
            LMAX = LHARR( ITYPE )
C
C Loop around L from LMIN, LMAX
C
            DO L = LMIN, LMAX
c             .
c             . ignore if our symmetry is not satisfied
c             .
              IF (     ISYM.EQ.1 .AND.
     1              MOD( (L-M), 2 ).NE.M2 
     2                                       ) GOTO 49
              IF (     ISYM.EQ.2 .AND.
     1              MOD( (L-M), 2 ).EQ.M2 
     2                                       ) GOTO 49
c             .
              DO ICS = 1, 2
                IF ( M.EQ.0 .AND. ICS.EQ.2 ) GOTO 49
c               .
c               . ok this harmonic DOES go in
c               . (provided we have enough room)
c               .
                NH = NH + 1
                IF ( NH.GT.NHMAX ) OK = .FALSE.
                IF ( OK ) MHT( NH ) = ITYPE
                IF ( OK ) MHL( NH ) = L
                IF ( OK .AND. ICS.EQ.1 ) MHM( NH ) = M
                IF ( OK .AND. ICS.EQ.2 ) MHM( NH ) = -M
              ENDDO
c             .
 49         CONTINUE
            ENDDO
C
 50       CONTINUE
          ENDDO
        ENDDO
      ENDDO
c
      IF ( OK ) RETURN
      PRINT *,' Subroutine HMINDA. Your specifications'
      PRINT *,' require ',NH,' harmonics. Maximum was set'
      PRINT *,' at ',NHMAX,'. Program aborted.'
c
      STOP
      END
C*********************************************************************
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
C*********************************************************************
C subroutine HarMonic File WriTe *************************************
C            -  -     -    -  -  *************************************
C Steve Gibbons Fri Nov 12 11:21:17 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Writes out the integer indices of spherical harmonic sets incl.    C
C the appropriate boundary conditions.                               C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Number of vector spherical harmonics.              C
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
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     NDCS      : Number of distinct finite difference schemes.      C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). Governs behaviour at inner     C
C                  boundary for finite diff. scheme ( is )           C
C                                                                    C
C  MHIBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( is ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                        where L = MHL( ih )                         C
C                                                                    C
C     MHOBC     : Dimension ( NDCS ). Governs behaviour at outer     C
C                  boundary for finite diff. scheme ( is )           C
C                                                                    C
C  MHOBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( is ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                        where L = MHL( ih )                         C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE HMFWT( NH, MHT, MHL, MHM, MHP, NDCS, MHIBC, MHOBC,
     1                  LU, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, MHT( * ), MHL( * ), MHM( * ), MHP( * ), NDCS,
     1        MHIBC( NDCS ), MHOBC( NDCS ), LU
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH, IWR, IS, IIBF, IOBF
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for writing
C
      IWR = 3
      CALL FOPEN ( LU, FNAME, IWR )
C
C Write number of spherical harmonics
C
      WRITE ( LU, 40 ) NH
C
C     Now loop around the harmonics and write out their
C     properties
C
      DO IH = 1, NH
        IS   = MHP( IH )
        IIBF = MHIBC( IS )
        IOBF = MHOBC( IS )
        WRITE ( LU, 41 ) MHT( IH ), MHL( IH ), MHM( IH ), IIBF, IOBF
      ENDDO
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 40   FORMAT(I5)
 41   FORMAT(I2,I4,I5,I3,I3)
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine X value ARRay WriTe *************************************
C            -       ---   -  -  *************************************
C Steve Gibbons Fri Nov 12 08:53:38 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Writes out the XARR array of abscissae to a file.                  C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                 Note that NR is not checked for correspondence to  C
C                 any other value - merely for being not greater     C
C                 than NRMAX.                                        C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C     IFORM     : Specifies how the x values are stored on the file. C
C                 Current values are:-                               C
C                                                                    C
C                   IFORM = 1 --> (5(1PD16.7))                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim ( * ) but length atleast NR. Location of       C
C                  radial grid nodes.                                C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XARRWT( NR, XARR, LU, FNAME, IFORM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, LU, IFORM
      CHARACTER *(*) FNAME
      DOUBLE PRECISION XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IWR
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine XARRWT.'
        PRINT *,' IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Open file for writing
C
      IWR = 3
      CALL FOPEN ( LU, FNAME, IWR )
C
C Write number of radial grid nodes
C
      WRITE ( LU, 40 ) NR, IFORM
C
C OK, so write X values ...
C
      IF ( IFORM.EQ.1 ) WRITE ( LU, 41 ) ( XARR( I ), I = 1, NR )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 40   FORMAT(I5,I5)
 41   FORMAT(5(1PD16.7))
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Solution Vector File WriTe ******************************
C            -        -      -    -  -  ******************************
C Steve Gibbons Sat Nov 13 14:30:24 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Writes out a solution vector to a file.                            C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1, 3. INDFUN = ( IR - 1 )*NH + IH       C
C                   IFORMF = 2, 4. INDFUN = ( IH - 1 )*NR + IR       C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C     IFORM     : Specifies how the x values are stored on the file. C
C                 Current values are:-                               C
C                                                                    C
C                   IFORM = 1 --> (5(1PD16.7))                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Dim ( * ) but length atleast NR*NH.                C
C                  Solution vector defined by INARR.                 C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVFWT( INARR, LU, IFORM, SV, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), LU, IFORM
      CHARACTER *(*) FNAME
      DOUBLE PRECISION SV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IWR, ILEN, IFORMF, NR, NH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFORMF = INARR( 1 )
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C
      ILEN   = NR*NH
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine SVFWT.'
        PRINT *,' IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Open file for writing
C
      IWR = 3
      CALL FOPEN ( LU, FNAME, IWR )
C
C Write iformf, nr, nh, iform
C  
       WRITE ( LU, 40 ) IFORMF, NR, NH, IFORM
C
C OK, so write X values ...
C
      IF ( IFORM.EQ.1 ) WRITE ( LU, 41 ) ( SV( I ), I = 1, ILEN )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 40   FORMAT(I5,I5,I5,I5)
 41   FORMAT(5(1PD16.7))
      RETURN
      END
C*********************************************************************

C********************************************************************
C SUBROUTINE File OPEN **********************************************
C            -    ---- **********************************************
C Steve Gibbons 14.4.97 (Adapted from Dan Gordon's Code)            C
C Routine modified 15.3.99
C___________________________________________________________________C
C Opens a file with number LU, name FNAME, access OACCES  and       C
C a flag IRW to indicate whether the file is to be read or written  C
C to. ( IRW=1 ==> read only, IRW=2 ==> write but only if the file   C
C doesn't already exist, IRW=3 ==> write regardless of whether file C
C exists or not.)						    C
C___________________________________________________________________C
C Input Variables :-						    C
C ===============   						    C
C  Integer							    C
C  -------							    C
C     LU	: File number					    C
C     IRW	: Read / Write Flag 				    C
C                  = 1 for read only		                    C
C                  = 2 for write (provided that file doesn't exist. C
C                  = 3 for write (regardless of existence of file.  C
C                  = 4 for append status.                           C
C  Character							    C
C  ---------							    C
C     FNAME	: File name					    C
C___________________________________________________________________C
C Working Variables :-						    C
C =================   						    C
C  Character							    C
C  ---------							    C
C     OACCES 	: Access flag - should be set to 'OLD' for read     C
C                          and 'UNKNOWN' for write                  C
C     CONTYN    : For a yes/no to IWR = 2 option.		    C
C     LABEL     : Null string to pass into FNAMER option.	    C
C  Logical							    C
C  -------							    C
C     LEXIST	: Existence of file. File present <==> LEXIST=.TRUE.C
C___________________________________________________________________C
C Subroutines Used :-                                               C
C ================                                                  C
C     FNAMER	: For the case of IRW = 2; trying to write to an    C
C		   existing file. Used if alternative filename is   C
C                   asked for.					    C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE FOPEN ( LU, FNAME, IRW)
      IMPLICIT NONE
C___________________________________________________________________C
C Variable declarations - Parameters ...............................C
      INTEGER LU, IRW
      CHARACTER *(*) FNAME
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      LOGICAL LEXIST
      CHARACTER *(7) OACCES
      CHARACTER *(1) CONTYN
      CHARACTER *(1) LABEL
C___________________________________________________________________C
      IF ( LU.EQ.0 ) THEN
         PRINT *,' Subroutine FOPEN'
         PRINT *,' I bet you ve forgotten to set LU ...??'
         PRINT *,' Think again and come back when you have'
         PRINT *,' remembered that LU must be a non zero integer!'
         PRINT *,' See you later. Bye for now!!'
         STOP
      ENDIF
C************************
C temporary code : SJG Thu Jun  1 08:07:06 BST 2000
C The Linux compiler will not allow file opening with
C lu.ge.100, so the following line should prevent it.
C 
      IF ( LU.GT.99 ) THEN
         PRINT *,' Subroutine FOPEN'
         PRINT *,' LU = ', LU,' too large.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C************************
C--------------
 600  CONTINUE
      INQUIRE (FILE=FNAME, EXIST=LEXIST)
C Case of read only - LEXIST must be = .TRUE.
      IF ( IRW.EQ.1 ) THEN
         OACCES = 'OLD'
         IF ( .NOT. LEXIST ) THEN
            PRINT *,' Subroutine FOPEN. You are trying to open an'
            PRINT *,' old file which does not exist.'
            PRINT *,' Filename = ',FNAME,' LU= ',LU
            PRINT *,' Program aborted.'
            STOP
         ELSE
            GOTO 500
         ENDIF
      ENDIF
C Case of write to file provided that it doesn't exist
      IF ( IRW.EQ.2 ) THEN
         OACCES = 'UNKNOWN'
         IF ( LEXIST ) THEN
            PRINT *, ' Subroutine FOPEN. You are trying to write'
            PRINT *, ' to an existing file with IRW set to 2.'
            PRINT *,' Filename = ',FNAME,' LU= ',LU
            PRINT *,' Do you wish to give an alternative FNAME?'
            PRINT *,' Type y or n.'
            READ ( 5, 267) CONTYN
 267         FORMAT (A)
            IF (CONTYN.NE.'y'.AND.CONTYN.NE.'Y') THEN
               PRINT *, ' Program Aborted.'
               STOP
            ELSE
               LABEL=' '
               CALL FNAMER ( FNAME, LABEL )
               GOTO 600
            ENDIF
         ELSE
            GOTO 500
         ENDIF
      ENDIF
C Case of write to file regardless of the existence of file.
      IF ( IRW.EQ.3 ) THEN
         OACCES = 'UNKNOWN'
         GOTO 500
      ENDIF
C Treat appendment case
      IF ( IRW.EQ.4 ) THEN
         OACCES = 'UNKNOWN'
         OPEN ( UNIT=LU , FILE=FNAME , STATUS=OACCES,
     1          ACCESS='APPEND', ERR=999 )
         RETURN
      ENDIF
C___________________________________________________________________C
C All the IRW cases as of 14.4.97 have now been covered.
      PRINT *,' Subroutine FOPEN. IRW must be set to 1, 2, 3 or 4.'
      PRINT *,' Program aborted.'
      STOP

 500  CONTINUE
      OPEN ( UNIT=LU , FILE=FNAME , STATUS=OACCES, ERR=999 )
      RETURN

 999  PRINT *,' Subroutine FOPEN. Error in opening file ',FNAME
      STOP

      END
C********************************************************************
C___________________________________________________________________C
C********************************************************************
C subroutine File CLOSE *********************************************
C            -    ----- *********************************************
C Steve Gibbons 14.4.97                                             C
C  ( note that this is essentially the routine of Dan Gordon        C
C___________________________________________________________________C
C Closes file with integer logical unit LU, filename FNAME.         C
C LABEL contains any other information regarding the nature of the  C
C the file.							    C
C___________________________________________________________________C
C Input Variables :-						    C
C ===============						    C
C  Integer							    C
C  -------							    C
C     LU	: Number of file.				    C
C  Character							    C
C  ---------							    C
C     FNAME	: Name of file. Undefined length		    C
C     LABEL	: Any further information. Undefined length         C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE FCLOSE ( LU, FNAME, LABEL )
      IMPLICIT NONE
C___________________________________________________________________C
C Variable Declarations - Parameters ...............................C
      INTEGER LU
      CHARACTER *(*) FNAME
      CHARACTER *(*) LABEL
C___________________________________________________________________C
C START OF PROGRAM *************************************************C
C___________________________________________________________________C
      IF ( LU.EQ.0 ) THEN
         PRINT *,' Subroutine FCLOSE '
         PRINT *,' I bet you ve forgotten to set LU ...??'
         PRINT *,' Think again and come back when you have'
         PRINT *,' remembered that LU must be a non zero integer!'
         PRINT *,' See you later. Bye for now!!'
         STOP
      ENDIF
C----------------------------
      CLOSE (UNIT=LU, STATUS='KEEP', ERR=989 )
      RETURN
C
 989  PRINT *,' Error.  Failed to close ', LABEL, ' file ', FNAME
      STOP
      END
C********************************************************************
C*********************************************************************
C subroutine SiMple ITerative SoLve **********************************
C            - -    --        - -   **********************************
C Steve Gibbons Thu Oct  7 16:26:59 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C If the double precision function FUNC, which has been declared     C
C EXTERNAL in the (sub)program which calls SMITSL, gives f( X )      C
C for a given value of X then SMITSL will iterate linearly towards   C
C a solution of f( X ) = 0.                                          C
C                                                                    C
C FUNC must have the argument list ( X, INTARR, DPRARR )             C
C                                                                    C
C with INTARR and DPRARR being respectively integer and double       C
C precision arrays of undefined length. INTARR and DPRARR are not    C
C referred to by SMITSL other than in the call to FUNC.              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FUNC      : Declared EXTERNALly. Has the form                  C
C                 FUNC( X, INTPAR, DPPAR )                           C
C                 (See above).                                       C
C                                                                    C
C     GUESS     : First value of X to be used.                       C
C     FACTOR    : Second value of X to be used to begin the          C
C                 iteration is FACTOR*GUESS                          C
C                                                                    C
C                 So if FACTOR = 1.1, then X2 = 1.1*X1               C
C                                                                    C
C     ERR       : User specified tolerance of error of the solution. C
C                 If f( X ) has a modulus less than ERR then X is    C
C                 returned as SOL.                                   C
C                                                                    C
C     SOL       : Approximation of solution.                         C
C                                                                    C
C     DPRARR    : Double precision array of parameters for FUNC.     C
C                                                  Dim (*)           C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ITMX      : Maximum number of iterations permitted.            C
C                                                                    C
C     INFO      : Output variable. If solution converges             C
C                 succesfully then INFO returns the number of        C
C                 iterations required to converge - this is          C
C                 an integer between 1 and ITMX.                     C
C                                                                    C
C                 A failiure will result in a negative value         C
C                 of INFO.                                           C
C                                                                    C
C                 INFO = -1 means two consecutive X values were      C
C                                                      identical.    C
C                                                                    C
C                 INFO = -2 means two consecutive f(X) values were   C
C                                                      identical.    C
C                                                                    C
C                 INFO = -3 means the maximum number of iterations   C
C                                          were exceeded.            C
C                                                                    C
C     INTPAR    : Integer array of parameters for FUNC. Dim (*)      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SMITSL( FUNC, GUESS, FACTOR, ITMX, INFO, ERR, SOL,
     1                   DPRARR, INTARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION FUNC, GUESS, FACTOR, ERR, SOL, DPRARR( * )
      INTEGER ITMX, INFO, INTARR( * )
      EXTERNAL FUNC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION X, FX, XOLD, FXOLD, DPM, DPC, FXSUB, XSUB,
     1                 TOL
      PARAMETER ( TOL = 1.0d-8 )
      INTEGER NOIT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Perform initial stage
C
      XOLD = GUESS
      FXOLD = FUNC( XOLD, INTARR, DPRARR )
C
C Now set X as a factor of FACTOR * GUESS
C
      X = FACTOR * GUESS
      FX = FUNC( X, INTARR, DPRARR )
C
      NOIT = 0
C
C Now begin iterative loop
C
 50   CONTINUE
      NOIT = NOIT + 1
      XSUB = X - XOLD
C
C Error trapping. If we have two identical consecutive
C values for x, this will result in a division by zero
C Set INFO = -1
C
      IF ( ABS( XSUB ).LT.TOL ) THEN
        INFO = -1
        RETURN
      ENDIF
C
C Set the gradient of our straight line to FXSUB/XSUB
C
      FXSUB = FX - FXOLD
      DPM = FXSUB/XSUB
C
C Error trapping. If gradient is zero then the function
C has had the same value for the previous iterations.
C Now this might be a (very unfortunate!) coincidence,
C or it could mean that the function is constant and
C no root can be found. In any case, this algorithm
C is unable to get out of this fix ...
C Set INFO = -2
C
      IF ( ABS( DPM ).LT.TOL ) THEN
        INFO = -2
        RETURN
      ENDIF
C
C Ok so gradient is non-zero and we can take the y
C intercept, DPC.
C
      DPC = FX - DPM*X
C
C So calculate our linearly extrapolated solution
C
      SOL = (-1.0d0)*DPC/DPM
C
C Set XOLD and FXOLD to X and FX
C
      XOLD = X
      FXOLD = FX
C
C Set X to SOL and calculate FX
C
      X = SOL
      FX = FUNC( X, INTARR, DPRARR )
C
C Check to see if we have converged
C If so then set INFO to the number of iterations taken
C and return
C
      IF ( ABS( FX ).LT.ERR ) THEN
        INFO = NOIT
        RETURN
      ENDIF
C
C See if we have reached the maximum number of iterations ...
C Set INFO = -3 and return
C
      IF ( NOIT.EQ.ITMX ) THEN
        INFO = -3
        RETURN
      ENDIF
C
C No problems, but no solution either, so go back
C to beginning of loop
C
      GOTO 50
C
      END
C*********************************************************************
C*********************************************************************
C subroutine Poloidal Magnetic Field Alpha and Beta Find *************
C            -        -        -     -         -    -    *************
C Steve Gibbons Mon Oct 11 11:59:12 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  If a poloidal radial function is expanded in the form             C
C                                                                    C
C  P( r ) = sum_i cos ( a_i r - b_i )                                C
C                                                                    C
C  when both inner and outer boundaries are insulating,              C
C  then for L, RI and RO (spherical harmonic degree, radius of       C
C  inner and outer boundary respectively) then a ( ALPHA ) and b     C
C  ( BETA ) must satisfy equations F1 = 0 and F2 = 0 as described    C
C  in Function PMFSEF.                                               C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     ALPHA     : Value a_i. Root to F1 and F2 in PMFSEF.            C
C                 On input is a guess for the i^{th} such root -     C
C                 This means that ALPHA is returned lying between    C
C                 i * pi/dist and ( i - 1 )*pi/dist.                 C
C                  ( dist = ro - ri )                                C
C                                                                    C
C     BETA      : Corresponding b_i value. Not constrained between   C
C                 any limits as is unique only up to periodicity.    C
C                                                                    C
C     RI        : Radius at inner boundary                           C
C     RO        : Radius at outer boundary                           C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     I         : Number of root. See above.                         C
C                                                                    C
C     L         : Spherical harmonic degree.                         C
C                                                                    C
C     ITMX      : Maximum iterations allowed to find a_i and b_i.    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PMFABF( ALPHA, BETA, I, RI, RO, L, ITMX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER I, L, ITMX
      DOUBLE PRECISION ALPHA, BETA, RI, RO
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NEQN
      PARAMETER ( NEQN = 2 )
      DOUBLE PRECISION PMFSEF, PI, ERR, XVEC( NEQN ), WORKV( NEQN ),
     1                 WORKA( NEQN, NEQN ), DPRARR( 2 ), PIN,
     2                 PINM1, DELPI, DIST
      INTEGER INFO, IWORK( NEQN ), INTARR( 1 ), NOIT
      EXTERNAL PMFSEF
      PARAMETER ( PI=3.14159265358979312D0, ERR = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( BETA.EQ.0.0d0 ) BETA = 1.5d0
C
      IF ( ITMX.LT.2 ) THEN
         PRINT *,' Subroutine PMFABF. ITMX = ', ITMX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C Fill arrays for sending to MNEWTR
C
      DPRARR( 1 ) = RI
      DPRARR( 2 ) = RO
      DIST = RO - RI
C
      INTARR( 1 ) = L
C
      PIN    = DBLE( I )*PI/DIST
      PINM1  = (DBLE( I - 1 )*PI + ERR)/DIST
C
C ALPHA and BETA are our initial guesses -
C there is no constraint upon BETA, although ALPHA
C must lie between ( I - 1 )*PI and I*PI
C otherwise, we will ignore this guess and loop around
C lot's of random guesses - this will be at line 50.
C
      IF ( ALPHA.LE.PINM1 .OR. ALPHA.GT.PIN ) GOTO 50
C
C Ok so our guess for ALPHA seems o.k.
C Let's try and solve using this and our guess for BETA
C
      XVEC( 1 ) = ALPHA
      XVEC( 2 ) = BETA
C
      CALL MNEWTR( PMFSEF, XVEC, NEQN, ITMX, ERR, INFO, WORKV,
     1             WORKA, IWORK, INTARR, DPRARR )
C
C Check to see if MNEWTR has found a (possibly a useless)
C solution ...
C
      IF ( INFO.EQ.-1 ) THEN
        PRINT *,' In MNEWTR, the LAPACK routine DGETRF'
        PRINT *,' returned a value of ', IWORK( 1 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( INFO.EQ.-2 ) THEN
        PRINT *,' In MNEWTR, the LAPACK routine DGETRS'
        PRINT *,' returned a value of ', IWORK( 1 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( INFO.EQ.-3 ) THEN
        PRINT *,' In MNEWTR, too many iterations were'
        PRINT *,' needed. ITMX = ', ITMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check our solution is in the right bounds
C If so, return ALPHA and BETA
C
      IF ( XVEC( 1 ).GT.PINM1 .AND. XVEC( 1 ).LE.PIN ) THEN
        ALPHA = XVEC( 1 )
        BETA  = XVEC( 2 )
        RETURN
      ENDIF
C
 50   CONTINUE
C
C O.k. our guesses weren't any good and so we will
C have to try and guess values ...
C
      DELPI = PI/(DBLE( ITMX + 1 )*DIST)
      DO NOIT = 1, ITMX
        XVEC( 1 ) = PIN - DBLE( NOIT )*DELPI
        XVEC( 2 ) = BETA
C
        CALL MNEWTR( PMFSEF, XVEC, NEQN, ITMX, ERR, INFO, WORKV,
     1               WORKA, IWORK, INTARR, DPRARR )
C
C Check to see if MNEWTR has found a (possibly a useless)
C solution ...
C
        IF ( INFO.EQ.-1 ) THEN
          PRINT *,' In MNEWTR, the LAPACK routine DGETRF'
          PRINT *,' returned a value of ', IWORK( 1 )
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
        IF ( INFO.EQ.-2 ) THEN
          PRINT *,' In MNEWTR, the LAPACK routine DGETRS'
          PRINT *,' returned a value of ', IWORK( 1 )
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
        IF ( INFO.EQ.-3 ) THEN
          PRINT *,' In MNEWTR, too many iterations were'
          PRINT *,' needed. ITMX = ', ITMX
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
C o.k. - so A solution was found. Let's see
C if it is in range ...
C
        IF ( XVEC( 1 ).GT.PINM1 .AND. XVEC( 1 ).LE.PIN ) THEN
          ALPHA = XVEC( 1 )
          BETA  = XVEC( 2 )
          RETURN
        ENDIF
C
 60   CONTINUE
      ENDDO
C
      PRINT *,' All you guesses for ALPHA are '
      PRINT *,' now exhausted and your solution '
      PRINT *,' still cannot be found.'
      STOP
      END
C*********************************************************************
C*********************************************************************
C subroutine File NAME giveR *****************************************
C            -    ----     - *****************************************
C Steve Gibbons 14.4.97 (Adapted from Dan Gordon's course.)          C
C____________________________________________________________________C
C Asks the user for a file name FNAME. Pretty simple really ...      C
C____________________________________________________________________C
C Input Variable :-						     C
C ==============   						     C
C  Character							     C
C  ---------							     C
C     LABEL	: Message arbitrary length  			     C
C____________________________________________________________________C
C Output Variable :-						     C
C ===============   						     C
C  Character							     C
C  ---------							     C
C     FNAME	: Filename arbitrary length			     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FNAMER ( FNAME, LABEL )
      IMPLICIT NONE
      CHARACTER *(*) FNAME
      CHARACTER *(*) LABEL

      PRINT *,' Please enter a filename for FNAME.'
      PRINT *, LABEL
      READ (5, 200) FNAME
 200  FORMAT (A)
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Multi-dimensional NEWTon raphson Solve ******************
C            -                 ----           -     ******************
C Steve Gibbons Sat Oct  9 18:11:31 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Say we have NEQN equations and NEQN unknowns.                      C
C For IEQN = 1, NEQN; we seek a solution ( x_1 , x_2, ... , x_NEQN ) C
C to the equations                                                   C
C                                                                    C
C f_IEQN( x_1 , x_2, ... , x_NEQN ) = 0.                             C
C                                                                    C
C FUNC is a double precision function with the calling sequence      C
C                                                                    C
C FUNC( NEQN, IEQN, ID, XVEC, INTARR, DPRARR )                       C
C                                                                    C
C FUNC returns f_IEQN( x_1 , x_2, ... , x_NEQN ) when ID = 0.        C
C If ID = ICMP, FUNC returns the partial derivative of f_IEQN        C
C with respect to x_ICMP.                                            C
C                                                                    C
C XVEC is a vector of length NEQN with X(icmp) = the icmp^{th} var.  C
C                                                                    C
C INTPAR and DPRARR are unlimited arrays of integer and double       C
C precision variables respectively which are not referred to by      C
C MNEWTS other than to supply information to FUNC.                   C
C                                                                    C
C                                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FUNC      : Declared EXTERNALly. Has the form                  C
C                 FUNC( NEQN, IEQN, ID, XVEC, INTARR, DPRARR )       C
C                 (see above).                                       C
C                                                                    C
C     XVEC      : Vector containing initial estimates of the         C
C                 unknowns. On return, XVEC contains the solution    C
C                 provided that convergence is attained.             C
C                 Dim ( NEQN )                                       C
C                                                                    C
C     ERR       : User specified tolerance of error of the solution. C
C                 If the norm of the r.h.s. vector is less than      C
C                 ERR, convergence will be judged to have been       C
C                 achieved.                                          C
C                                                                    C
C     WORKV     : Work array dim ( NEQN )                            C
C     WORKA     : Work array dim ( NEQN, NEQN )                      C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NEQN      : Number of equations and unknowns.                  C
C     ITMX      : Maximum number of iterations permitted.            C
C                                                                    C
C     INFO      : Output variable. If solution converges             C
C                 succesfully then INFO returns the number of        C
C                 iterations required to converge - this is          C
C                 an integer between 1 and ITMX.                     C
C                                                                    C
C                 A failiure will result in a negative value         C
C                 of INFO.                                           C
C                                                                    C
C                 INFO = -1: Jacobian matrix was singular.           C
C                 The LAPACK routine DGETRF was unable to perform    C
C                 an LU decomposition.                               C
C                 The value of 'INFO' returned by DGETRF             C
C                 is returned from MNEWTR in IWORK( 1 )              C
C                                                                    C
C                 INFO = -2:                                         C
C                 The LAPACK routine DGETRS was unable to solve      C
C                 the matrix equation.                               C
C                 The value of 'INFO' returned by DGETRS             C
C                 is returned from MNEWTR in IWORK( 1 )              C
C                                                                    C
C                 INFO = -3: Maximum number of iterations exceeded.  C
C                                                                    C
C     IWORK     : Work array dim ( NEQN ).                           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MNEWTR( FUNC, XVEC, NEQN, ITMX, ERR, INFO, WORKV,
     1                   WORKA, IWORK, INTARR, DPRARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NEQN, ITMX, INFO, IWORK( NEQN ), INTARR( * )
      DOUBLE PRECISION FUNC, XVEC( NEQN ), WORKV( NEQN ), ERR,
     1                 WORKA( NEQN, NEQN ), DPRARR( * )
      EXTERNAL FUNC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION DPNORM
      INTEGER NOIT, ID, IEQN, IERR, NRHS
      CHARACTER *(1) TRANS
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NRHS = 1
      TRANS = 'N'
      NOIT = 0
C
 50   CONTINUE
      NOIT = NOIT + 1
C
C Check maximum iterations haven't been exceeded
C
      IF ( NOIT.GT.ITMX ) THEN
        INFO = -3
        RETURN
      ENDIF
C
C O.k. - we form the right hand side
C Store this in WORKV
C Also calculate the solution norm
C
      ID = 0
      DPNORM = 0.0d0
      DO IEQN = 1, NEQN
        WORKV( IEQN ) =
     1    (-1.0d0)*FUNC( NEQN, IEQN, ID, XVEC, INTARR, DPRARR )
        DPNORM = DPNORM + ABS( WORKV( IEQN ) )
      ENDDO
C
C Return if change in solution is sufficiently small
C
      IF ( DPNORM.LT.ERR ) THEN
        INFO = NOIT
        RETURN
      ENDIF
C
C Now form the Jacobian matrix
C WORKA( i, j ) must contain the derivative of F_i with
C respect to x_j
C
      DO ID = 1, NEQN
        DO IEQN = 1, NEQN
          WORKA( IEQN, ID ) = FUNC( NEQN, IEQN, ID, XVEC,
     1                              INTARR, DPRARR )
        ENDDO
      ENDDO
C
C Now perform an LU decomposition upon the matrix WORKA
C Use the LAPACK routine DGETRF
C
      CALL DGETRF( NEQN, NEQN, WORKA, NEQN, IWORK, IERR )
C
C Return if LU decomposition has been unsuccessful
C
      IF ( IERR.NE.0 ) THEN
        IWORK( 1 ) = IERR
        INFO = -1
        RETURN
      ENDIF
C
C Ok - the LU decomposition seems to have gone
C without a problem. So let's solve the equation
C WORKA . DX = WORKV
C Use LAPACK routine DGETRS
C
      CALL DGETRS( TRANS, NEQN, NRHS, WORKA, NEQN, IWORK,
     1             WORKV, NEQN, IERR )
C
C Return if solution has been unsuccessful
C
      IF ( IERR.NE.0 ) THEN
        IWORK( 1 ) = IERR
        INFO = -2
        RETURN
      ENDIF
C
C Update solution.
C
      DO ID = 1, NEQN
        XVEC( ID ) = XVEC( ID ) + WORKV( ID )
      ENDDO
C
      GOTO 50
C
      END
C*********************************************************************
C*********************************************************************
C function Poloidal Velocity No Slip Radial Function *****************
C          -        -        -  -    -      -        *****************
C Steve Gibbons Fri Oct  8 09:46:45 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C This function is a series of functions as detailed in              C
C Chandrasekhar Appendix V, page 635                                 C
C 'Hydodynamic and Hydromagnetic Stability' 1981, Dover, New York.   C
C                                                                    C
C PVNSRF( r ) = \sum_{m=1}^{m=MAXM} \left(                           C
C                 c_m C_m( x ) + s_m S_m( x ) \right)                C
C                                                                    C
C where c_m = COEFS( 2m - 1 ), s_m = COEFS( 2m ),                    C
C                                                                    C
C x = ( r - ri )/(ro - ri ) - 0.5d0,                                 C
C                                                                    C
C            COSH( lambda_m*x )       COS( lambda_m*x )              C
C C_m( x ) = ----------------     -    ----------------              C
C             COSH( 0.5 * x )           COS( 0.5 * x )               C
C                                                                    C
C             SINH( myu_m*x )           SIN( myu_m*x )               C
C S_m( x ) = ----------------     -    ----------------              C
C             SINH( 0.5 * x )           SIN( 0.5 * x )               C
C                                                                    C
C  lambda_m and myu_m are respectively the m^{th} non-zero roots of  C
C                                                                    C
C tanh ( 0.5*lambda ) + tan( 0.5*lambda ) = 0 and                    C
C                                                                    C
C coth ( 0.5*myu ) - cot( 0.5*myu ) = 0                              C
C                                                                    C
C which are calculated by calls to SMITSL with different options     C
C on the EXTERNAL function CHORCH.                                   C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     R         : Radius                                             C
C     RI        : Radius at inner boundary                           C
C     RO        : Radius at outer boundary                           C
C     COEFS     : c_m in COEFS( 2m - 1 )                             C
C                 s_m in COEFS( 2m )                                 C
C     DPRPA     : D.p. array( * ). Not referred to                   C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     MMAX      : Highest value of M                                 C
C     INTPA     : Int array( * ). Not referred to                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PVNSRF( R, RI, RO, MAXM, COEFS, INTPA, DPRPA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION PVNSRF, R, RI, RO, COEFS( * ), DPRPA( * )
      INTEGER MAXM, INTPA( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER M, IT, IND, INTARR( 1 ), ITMX, INFO
      DOUBLE PRECISION SOL, ERR, X, CHORCH, CHORCF, CHORSF, CCF, SCF,
     1                 GUESS, FACTOR, PI
      PARAMETER (PI=3.14159265358979312D0, ERR = 1.0d-9, ITMX = 40 )
      EXTERNAL CHORCH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      FACTOR = 1.02d0
      IF ( ABS(RO-RI).LT.ERR ) THEN
        PRINT *,' Function PVNSRF.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      IF ( R.LT.RI .OR. R.GT.RO ) THEN
        PRINT *,' Function PVNSRF. R = ', R
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Convert R (between RI and RO) to X
C     . (between -0.5 and 0.5 )
C     .
      X = (R-RI)/(RO-RI) - 0.5d0
C     .
      PVNSRF = 0.0d0
C     .
C     . Loop around values of M
C     .
      DO M = 1, MAXM
C       .
C       . it = 1  --> C_m( x )  (see CHORCF)
C       .
        IT = 1
        INTARR( 1 ) = IT
        IND = 2*M - 1
        CCF = COEFS( IND )
C
C Try to locate the M^th zero of 
C tanh( 0.5*x ) + tan( 0.5*x )
C initial guess for this is (2m-0.5)*pi
C
        GUESS = (2.0d0*DBLE( M ) - 0.5d0 )*PI
        CALL SMITSL( CHORCH, GUESS, FACTOR, ITMX, INFO, ERR,
     1               SOL, DPRPA, INTARR )
C
C Check that we have successfully located a root
C
        IF ( INFO.LT.0 ) THEN
          PRINT *,' Function PVNSRF.'
          PRINT *,' SMITSL has returned INFO = ', INFO
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
C So SOL now contains lambda_m; so add contribution
C to PVNSRF
C
        PVNSRF = PVNSRF + CCF*CHORCF( X, SOL )
C       .
C       . it = 2  --> S_m( x )  (see CHORSF)
C       .
        IT = 2
        INTARR( 1 ) = IT
        IND = 2*M
        SCF = COEFS( IND )
C       .
C Try to locate the M^th zero of 
C coth( 0.5*x ) - cot( 0.5*x )
C initial guess for this is (2m+0.5)*pi
C
        GUESS = (2.0d0*DBLE( M ) + 0.5d0 )*PI
        CALL SMITSL( CHORCH, GUESS, FACTOR, ITMX, INFO, ERR,
     1               SOL, DPRPA, INTARR )
C
C Check that we have successfully located a root
C
        IF ( INFO.LT.0 ) THEN
          PRINT *,' Function PVNSRF.'
          PRINT *,' SMITSL has returned INFO = ', INFO
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C
C So SOL now contains lambda_m; so add contribution
C to PVNSRF
C
        PVNSRF = PVNSRF + SCF*CHORSF( X, SOL )
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function Toroidal Velocity Stress Free Radial Function *************
C          -        -        -      -    -      -        *************
C Steve Gibbons Fri Oct  8 16:57:14 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Fills a function f( r ) = \sum_{m = 0}^{MAXM}                      C
C                 r cos( m pi (r - ri)/(ro - ri) )                   C
C                                                                    C
C This satisfies the boundary conditions                             C
C                                                                    C
C       d [ f ]                                                      C
C      -- [---] = 0     at    r = ri    and    r = ro                C
C      dr [ r ]                                                      C
C                                                                    C
C as required by the toroidal velocity condition with stress free    C
C boundary conditions.                                               C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     R         : Radius                                             C
C     RI        : Radius at inner boundary                           C
C     RO        : Radius at outer boundary                           C
C     COEFS     : c_m in COEFS( m + 1 )                              C
C     DPRPA     : D.p. array( * ). Not referred to                   C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     MMAX      : Highest value of M                                 C
C     INTPA     : Int array( * ). Not referred to                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION TVSFRF( R, RI, RO, MAXM, COEFS, INTPA, DPRPA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION TVSFRF, R, RI, RO, COEFS( * ), DPRPA( * )
      INTEGER MAXM, INTPA( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER M
      DOUBLE PRECISION PI, FAC, ERR, F2
      PARAMETER (PI=3.14159265358979312D0, ERR = 1.0d-9)
      EXTERNAL CHORCH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( ABS(RO-RI).LT.ERR ) THEN
        PRINT *,' Function TVSFRF.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      IF ( R.LT.RI .OR. R.GT.RO ) THEN
        PRINT *,' Function TVSFRF. R = ', R
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      FAC = PI*(R - RI)/(RO-RI)
C     .
      TVSFRF = 0.0d0
C     .
C     . Loop around values of M
C     .
      DO M = 0, MAXM
C       .
        F2 = FAC*DBLE( M )
        TVSFRF = TVSFRF + R*DCOS( F2 )*COEFS( M + 1)
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function Poloidal Velocity Stress Free Radial Function *************
C          -        -        -      -    -      -        *************
C Steve Gibbons Fri Oct  8 15:24:51 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Fills a function f( r ) = \sum_{m = 0}^{MAXM}                      C
C                   sin( m pi (r - ri)/(ro - ri) )                   C
C                                                                    C
C This satisfies the boundary conditions                             C
C                                                                    C
C p( ri ) = p( ro ) = p''( ri ) = p''( ro ) = 0                      C
C                                                                    C
C as required by the poloidal velocity condition with stress free    C
C boundary conditions.                                               C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     R         : Radius                                             C
C     RI        : Radius at inner boundary                           C
C     RO        : Radius at outer boundary                           C
C     COEFS     : c_m in COEFS( m )                                  C
C     DPRPA     : D.p. array( * ). Not referred to                   C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     MMAX      : Highest value of M                                 C
C     INTPA     : Int array( * ). Not referred to                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PVSFRF( R, RI, RO, MAXM, COEFS, INTPA, DPRPA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION PVSFRF, R, RI, RO, COEFS( * ), DPRPA( * )
      INTEGER MAXM, INTPA( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER M
      DOUBLE PRECISION PI, FAC, ERR, F2
      PARAMETER (PI=3.14159265358979312D0, ERR = 1.0d-9)
      EXTERNAL CHORCH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( ABS(RO-RI).LT.ERR ) THEN
        PRINT *,' Function PVSFRF.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      IF ( R.LT.RI .OR. R.GT.RO ) THEN
        PRINT *,' Function PVSFRF. R = ', R
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      FAC = PI*(R - RI)/(RO-RI)
C     .
      PVSFRF = 0.0d0
C     .
C     . Loop around values of M
C     .
      DO M = 1, MAXM
C       .
        F2 = FAC*DBLE( M )
        PVSFRF = PVSFRF + DSIN( F2 )*COEFS( M )
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function Poloidal Magnetic Field (Insulated) Radial Function *******
C          -        -        -      -          -      -        *******
C Steve Gibbons Mon Oct 11 11:59:12 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  If a poloidal radial function is expanded in the form             C
C                                                                    C
C  P( r ) = sum_i c_i cos ( a_i r - b_i )                            C
C                                                                    C
C  when both inner and outer boundaries are insulating,              C
C  then for L, RI and RO (spherical harmonic degree, radius of       C
C  inner and outer boundary respectively) then PMFIRF returns the    C
C  value of P( r ) where c_i is given by COEFS( i ) and the a_i and  C
C  b_i are determined by the subroutine PMFABF.                      C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     R         : Radius                                             C
C     RI        : Radius at inner boundary                           C
C     RO        : Radius at outer boundary                           C
C     COEFS     : c_i in COEFS( i )                                  C
C     DPRPA     : D.p. array( * ). Not referred to                   C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IMAX      : Highest value of I                                 C
C     INTPA     : Int array( * ). INTPA( 1 ) contains sph. harm      C
C                                              degree, L.            C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMFIRF( R, RI, RO, IMAX, COEFS, INTPA, DPRPA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION PMFIRF, R, RI, RO, COEFS( * ), DPRPA( * )
      INTEGER          IMAX, INTPA( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, L, ITMX
      DOUBLE PRECISION COEF, ALPHA, BETA, PI, ERR, BM1, BM2, DIST
      PARAMETER ( PI=3.14159265358979312D0, ERR = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( ABS(RO-RI).LT.ERR ) THEN
        PRINT *,' Function PMFIRF.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( R.LT.RI .OR. R.GT.RO ) THEN
        PRINT *,' Function PMFIRF. R = ', R
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DIST = RO - RI
C
      L      = INTPA( 1 )
      ITMX   = 60
      PMFIRF = 0.0d0
C
C Initialise guesses
C
      ALPHA = PI*0.75d0/DIST
      BETA  = PI*0.5d0
C
C Loop around the i
C
      DO I = 1, IMAX
C
        COEF = COEFS( I )
C
        IF ( I.GT.2 ) BETA = 2.0d0*BM1 - BM2
C
C       Calculate the constants a_i and b_i
C
        CALL PMFABF( ALPHA, BETA, I, RI, RO, L, ITMX )
C
C       Now add onto our total for PMFIRF
C
        PMFIRF = PMFIRF + COEF*COS( ALPHA*R - BETA )
C
C Update guesses
C
        BM2 = BM1
        BM1 = BETA
        ALPHA = ALPHA + PI/DIST
C
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function CHandrasekhar's ORthogonal CHaracteristic function ********
C          --              --         --                      ********
C Steve Gibbons Thu Oct  7 15:00:53 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C If IFLAG = INTARR( 1 ) = 1 then                                    C
C              CHORCH = TANH( 0.5*X ) + TAN( 0.5*X )                 C
C                                                                    C
C If IFLAG = INTARR( 1 ) = 2 then                                    C
C              CHORCH = 1.0d0/TANH( 0.5*X ) - 1.0d0/TAN( 0.5*X )     C
C                                                                    C
C DPRARR is not referred to by CHORCH. It is merely there so that    C
C roots to this equation may be found by SMITSL.                     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION CHORCH( X, INTARR, DPRARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION CHORCH, DPRARR( * ), X
      INTEGER INTARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFLAG
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFLAG = INTARR( 1 )
C     .
C     . Check IFLAG is valid
C     .
      IF ( IFLAG.NE.1 .AND. IFLAG.NE.2 ) THEN
        PRINT *,' Function CHORCH.'
        PRINT *,' INTARR( 1 ) = ', IFLAG
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IFLAG.EQ.1 ) CHORCH = TANH( 0.5*X ) + TAN( 0.5*X )
      IF ( IFLAG.EQ.2 ) CHORCH = 1.0d0/TANH( 0.5*X ) -
     1                 1.0d0/TAN( 0.5*X )
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function CHandrasekhar's ORthogonal Cosine Function ****************
C          --              --         -      -        ****************
C Steve Gibbons Thu Oct  7 15:00:53 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C For a value of lambda and x, CHORCF returns the value of           C
C                                                                    C
C          COSH( lambda*x )       COS( lambda*x )                    C
C C( x ) = ----------------   -   ----------------                   C
C          COSH( 0.5 * x )        COS( 0.5 * x )                     C
C                                                                    C
C as defined in Appendix V (page 635) of S. Chandrasekhar,           C
C 'Hydodynamic and Hydromagnetic Stability' 1981, Dover, New York.   C
C                                                                    C
C C( x ) = 0 at x = -0.5 and x = 0.5 and                             C
C dC/dx = 0 at x = -0.5 and x = 0.5 if lambda is the root of         C
C the characteristic equation                                        C
C                                                                    C
C    tanh ( 0.5*lambda ) + tan( 0.5*lambda ) = 0                     C
C                                                                    C
C Function and both input variables are double precision scalars.    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION CHORCF( X, LAMBDA)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION CHORCF, X, LAMBDA
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TOL, Q1, Q2, Q3, Q4, LAMX, HL
      PARAMETER (TOL=1.0d-10)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      LAMX = LAMBDA*X
      HL   = 0.5d0*LAMBDA
      Q1   = COSH( LAMX )
      Q2   = COSH( HL )
      Q3   = DCOS( LAMX )
      Q4   = DCOS( HL )
C
      IF ( ABS( Q2 ).LT.TOL ) THEN
        PRINT *,' Function CHORCF.'
        PRINT *,' COSH(',HL,') = ', Q2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ABS( Q4 ).LT.TOL ) THEN
        PRINT *,' Function CHORCF.'
        PRINT *,' COS(',HL,') = ', Q4
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      CHORCF = Q1/Q2 - Q3/Q4
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function CHandrasekhar's ORthogonal Sine Function ******************
C          --              --         -    -        ******************
C Steve Gibbons Thu Oct  7 15:00:53 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C For a value of myu and x, CHORSF returns the value of              C
C                                                                    C
C          SINH( myu*x )           SIN( myu*x )                      C
C C( x ) = ----------------   -   ----------------                   C
C          SINH( 0.5 * x )        SIN( 0.5 * x )                     C
C                                                                    C
C as defined in Appendix V (page 635) of S. Chandrasekhar,           C
C 'Hydodynamic and Hydromagnetic Stability' 1981, Dover, New York.   C
C                                                                    C
C C( x ) = 0 at x = -0.5 and x = 0.5 and                             C
C dC/dx = 0 at x = -0.5 and x = 0.5 if myu is the root of            C
C the characteristic equation                                        C
C                                                                    C
C    coth ( 0.5*myu ) - cot( 0.5*myu ) = 0                           C
C                                                                    C
C Function and both input variables are double precision scalars.    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION CHORSF( X, MYU)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION CHORSF, X, MYU
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TOL, Q1, Q2, Q3, Q4, MYUX, HL
      PARAMETER (TOL=1.0d-10)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      MYUX = MYU*X
      HL   = 0.5d0*MYU
      Q1   = SINH( MYUX )
      Q2   = SINH( HL )
      Q3   = DSIN( MYUX )
      Q4   = DSIN( HL )
C
      IF ( ABS( Q2 ).LT.TOL ) THEN
        PRINT *,' Function CHORSF.'
        PRINT *,' SINH(',HL,') = ', Q2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ABS( Q4 ).LT.TOL ) THEN
        PRINT *,' Function CHORSF.'
        PRINT *,' SIN(',HL,') = ', Q4
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      CHORSF = Q1/Q2 - Q3/Q4
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function Poloidal Magnetic Field Spectral Expansion Function *******
C          -        -        -     -        -         -        *******
C Steve Gibbons Sun Oct 10 15:54:27 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Auxiliary function to MNEWTR to solve for alpha and beta           C
C when F1( alpha, beta) = 0  and F2( alpha, beta) = 0.               C
C                                                                    C
C Here F1( alpha, beta) = -ri alpha sin ( alpha ri - beta )          C
C                         - l cos ( alpha ri - beta )                C
C                                                                    C
C Here F2( alpha, beta) = -ro alpha sin ( alpha ro - beta )          C
C                         + (l + 1 ) cos ( alpha ri - beta )         C
C                                                                    C
C Also gives partial derivatives of F1 and F2 with respect to        C
C alpha and beta.                                                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     XVEC      : XVEC( 1 ) contains alpha.                          C
C                 XVEC( 2 ) contains beta.                           C
C                                                                    C
C     DPRARR    : DPRARR( 1 ) contains RI.                           C
C                 DPRARR( 2 ) contains RO. ( Array of dim ( * ) )    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NEQN      : Number of equations and unknowns. Must be 2.       C
C                                                                    C
C     IEQN      : Number of current equation. Must be either 1 or 2. C
C                 (Corresponds to F1 and F2 above resp.)             C
C                                                                    C
C     ID        : ID = 0 --> return value of F( IEQN )               C
C                 ID = 1 --> return value of dF( IEQN )/d alpha      C
C                 ID = 2 --> return value of dF( IEQN )/d beta       C
C                                                                    C
C     INTARR    : INTARR( 1 ) contains L. Array of dim ( * )         C
C                                                                    C
C____________________________________________________________________C
C Output :-                                                          C
C ======                                                             C
C  Double Precision                                                  C
C  ----------------                                                  C
C     PMFSEF    : Value of function or partial derivative.           C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMFSEF( NEQN, IEQN, ID, XVEC, INTARR, DPRARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NEQN, IEQN, ID, INTARR( * )
      DOUBLE PRECISION PMFSEF, XVEC( * ), DPRARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L
      DOUBLE PRECISION FACINN, FACOUT, RI, RO, ALPHA, BETA
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First, trap any bad values of NEQN, IEQN and ID.
C
      IF ( NEQN.NE.2 ) THEN
         PRINT *,' Function PMFSEF.'
         PRINT *,' NEQN = ', NEQN,' and must be 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IEQN.NE.1 .AND. IEQN.NE.2 ) THEN
         PRINT *,' Function PMFSEF.'
         PRINT *,' IEQN = ', IEQN,' and must be 1 or 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( ID.NE.0 .AND. ID.NE.1 .AND. ID.NE.2 ) THEN
         PRINT *,' Function PMFSEF.'
         PRINT *,' ID = ', ID,' and must be 0, 1 or 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C ok - our command is valid
C
      L      = INTARR( 1 )
      RI     = DPRARR( 1 )
      RO     = DPRARR( 2 )
      ALPHA  = XVEC( 1 )
      BETA   = XVEC( 2 )
C
      FACINN = ALPHA*RI - BETA
      FACOUT = ALPHA*RO - BETA
C
C First do the functions
C
      IF ( ID.EQ.0 ) THEN
        IF ( IEQN.EQ.1 ) THEN
          PMFSEF = (-1.0d0)*RI*ALPHA*SIN( FACINN ) -
     1             DBLE( L )*COS( FACINN )
        ENDIF
        IF ( IEQN.EQ.2 ) THEN
          PMFSEF = (-1.0d0)*RO*ALPHA*SIN( FACOUT ) +
     1             DBLE( L + 1 )*COS( FACOUT )
        ENDIF
      ENDIF
C
C Now do the derivatives w.r.t. alpha
C
      IF ( ID.EQ.1 ) THEN
        IF ( IEQN.EQ.1 ) THEN
          PMFSEF = (-1.0d0)*RI*RI*ALPHA*COS( FACINN ) -
     1             RI*SIN( FACINN ) +
     2             DBLE( L )*RI*SIN( FACINN )
        ENDIF
        IF ( IEQN.EQ.2 ) THEN
          PMFSEF = (-1.0d0)*RO*RO*ALPHA*COS( FACOUT ) -
     1             RO*SIN( FACOUT ) -
     2             DBLE( L + 1 )*RO*SIN( FACOUT )
        ENDIF
      ENDIF
C
C Now do the derivatives w.r.t. beta
C
      IF ( ID.EQ.2 ) THEN
        IF ( IEQN.EQ.1 ) THEN
          PMFSEF = RI*ALPHA*COS( FACINN ) -
     1             DBLE( L )*SIN( FACINN )
        ENDIF
        IF ( IEQN.EQ.2 ) THEN
          PMFSEF = RO*ALPHA*COS( FACOUT ) +
     1             DBLE( L + 1 )*SIN( FACOUT )
        ENDIF
      ENDIF
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C integer function INDex FUNction ************************************
C                  ---   ---      ************************************
C Steve Gibbons Thu Sep 16 11:13:38 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Returns the location in the solution vector of the value of the    C
C IR^{th} grid node of the IH^{th} harmonic.                         C
C                                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Number of radial grid node.                        C
C     IH        : Number of harmonic.                                C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1/3. INDFUN = ( IR - 1 )*NH + IH        C
C                   IFORMF = 2/4. INDFUN = ( IH - 1 )*NR + IR        C
C                                                                    C
C  Here, NH is the TOTAL number of harmonics in the solution         C
C  vector - or atleast the part of which is visible to that part     C
C  of the program. NR is the number of radial grid nodes             C
C  corresponding to that harmonic (which for cases IFORMF = 1 and    C
C  IFORMF = 2 is identical for all harmonics - more complicated      C
C  options (for example magnetic fields which have to be resolved    C
C  beyond the region of fluid flow) may be added later and this      C
C  routine should be flexible to all possibilities with extra        C
C  constraints being added in other elements of INARR.               C
C                                                                    C
C  IFORMF = 1/3 is the option likely to be used for the purpose of   C
C  solution as it allows the banded formation of a matrix.           C
C                                                                    C
C  IFORMF = 2/4 is the option likely to be used for the purpose of   C
C  displaying solutions as it stores adjacent nodes for each         C
C  harmonic together.                                                C
C                                                                    C
C                 INARR( 2 ) = NR. See above.                        C
C                 INARR( 3 ) = NH. See above.                        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION INDFUN ( IR, IH, INARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INDFUN, IR, IH, INARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFORMF, NR, NH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFORMF = INARR( 1 )
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C
      IF ( IR.LT.1 .OR. IR.GT.NR ) THEN
        PRINT *,' Function INDFUN. IR = ', IR
        PRINT *,' NR = ', NR,'. Program aborted.'
        STOP
      ENDIF
C
      IF ( IH.LT.1 .OR. IH.GT.NH ) THEN
        PRINT *,' Function INDFUN. IH = ', IH
        PRINT *,' NH = ', NH,'. Program aborted.'
        STOP
      ENDIF
C
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.3 ) THEN
        INDFUN = ( IR - 1 )*NH + IH
        RETURN
      ENDIF
C
      IF ( IFORMF.EQ.2 .OR. IFORMF.EQ.4 ) THEN
        INDFUN = ( IH - 1 )*NR + IR
        RETURN
      ENDIF
C
      PRINT *,' Function INDFUN. IFORMF = ', IFORMF
      PRINT *,' Not current option. Program aborted.'
      STOP
      END
C*********************************************************************
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*inv( A' )*B.
*
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*inv( A' ).
*
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DGER  .
*
      END
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A' * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, JP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IP, IX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSWAP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.EQ.0 )
     $   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
*     End of ILAENV
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
