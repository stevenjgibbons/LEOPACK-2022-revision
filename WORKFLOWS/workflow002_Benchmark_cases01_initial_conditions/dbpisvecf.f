C*********************************************************************
C                                                                    C
C Steve Gibbons - Dynamo Benchmark Project Initial State             C
C                 VECtor Form                                        C
C                                                                    C
C Fri Feb  4 10:39:33 GMT 2000                                       C
C                                                                    C
C                                                                    C
C*********************************************************************
      PROGRAM dbpisvecf
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, NHMAX, LHMAX, NDCS
      PARAMETER ( NRMAX = 200, NHMAX = 1200, LHMAX = 100,
     1            NDCS = LHMAX + 2 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ), MHP( NHMAX ),
     1        MMODES( LHMAX + 1 ), INARR( 3 ), ISYMA( 5 ), LHARR( 5 ),
     2        MHBCS( NDCS )
      DOUBLE PRECISION XARR( NRMAX ), VEC( NRMAX*NHMAX )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER ISP, NR, MMAX, MINC, ISYMV, ISYMM, LHM, LHV, M, NH,
     1        IFLAG, NMODES, I, LU, ILEN, IFORM
      DOUBLE PRECISION ASP, RI, RO
      CHARACTER *(80) FNAME,ROOT
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IFORM = 1 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      LU = 91
C
      DO I = 1, 80
        FNAME (I:I) = ' '
        ROOT  (I:I) = ' '
      ENDDO
C
      PRINT *,' Radial set up. Enter NR, ASP, ISP '
      PRINT *,' NR = Number of grid nodes: Maximum = ', NRMAX
      PRINT *,' ASP = ratio of ri/ro '
      PRINT *,' ISP = 1: equally spaced nodes. '
      PRINT *,' ISP = 2: Checbyshev polynomial nodes. '
      PRINT *,' ------------------------------------- '
      READ (5,*) NR, ASP, ISP
C
      IF ( ASP.LT.0.0d0 .OR. ASP.GE.1.0d0 ) THEN
        PRINT *,' ASP = ',ASP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NR.LT.10 .OR. NR.GT.NRMAX ) THEN
        PRINT *,' NR = ', NR,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ISP.NE.1 .AND. ISP.NE.2 ) THEN
        PRINT *,' ISP = ', ISP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RO = 1.0d0/( 1.0d0 - ASP )
      RI = RO - 1.0d0
C
      IF ( ISP.EQ.1 ) THEN
        CALL ESNAAS( NR, XARR, RI, RO )
      ELSE
        CALL ZCPAAS( NR, XARR, RI, RO )
      ENDIF
C
      PRINT *,' Lateral resolution. Enter LHV, LHM, MINC, MMAX '
      PRINT *,' LHV  - maximum degree for velocity/temp. expansion.'
      PRINT *,' LHM  - maximum degree for mag. field expansion.'
      PRINT *,' MINC - increment for wavenumber '
      PRINT *,' MMAX - maximum wavenumber '
      PRINT *,' e.g. MINC = 2, MMAX = 6 will select m = 0, 2, 4 6.'
      PRINT *,' --------------------------------------------------'
      READ (5, * ) LHV, LHM, MINC, MMAX
C
      IF ( MINC.GT.4 .OR. MMAX.LT.4 .OR. MINC.EQ.3 .OR.
     1     MINC.LT.1 ) THEN
        PRINT *,' MINC = ', MINC,' MMAX = ', MMAX
        PRINT *,' Wavenumbers must include M = 4!!'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      LHARR( 1 ) = LHV
      LHARR( 2 ) = LHV
      LHARR( 3 ) = LHV
      LHARR( 4 ) = LHM
      LHARR( 5 ) = LHM
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
      PRINT *,' Symmetry considerations. Enter ISYMV, ISYMM.'
      PRINT *,' ISYMV = 1 for ES vel/temp, 2 for EA, 3 for both.'
      PRINT *,' ISYMM = 1 for ES mag, f.  2 for EA, 3 for both.'
      PRINT *,' e.g. a symmetric velocity field and dipole magnetic'
      PRINT *,' field would be given by ISYMV = 1, ISYMM = 2.'
      PRINT *,' -------------------------------------------------- '
      READ ( 5, * ) ISYMV, ISYMM
C
      ISYMA( 1 ) = ISYMV
      ISYMA( 2 ) = ISYMV
      ISYMA( 3 ) = ISYMV
      ISYMA( 4 ) = ISYMM
      ISYMA( 5 ) = ISYMM
C
      CALL HMINDA( LHARR, ISYMA, NMODES, MMODES, NH, NHMAX,
     1                   MHT, MHL, MHM, LHMAX )
      PRINT *,' Total harmonics selected = ', NH
      IF ( NH.EQ.0 ) THEN
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( LHM.EQ.0 ) THEN
        IFLAG = 2
      ELSE
        IFLAG = 1
      ENDIF
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C
      DO I = 1, LHM
        MHBCS( I ) = 7
      ENDDO
      MHBCS( LHM + 1 ) = 2
      MHBCS( LHM + 2 ) = 4
C
      DO I = 1, NH
        IF ( MHT( I ).EQ.1 ) MHP( I ) = LHM + 2
        IF ( MHT( I ).EQ.2 ) MHP( I ) = LHM + 1
        IF ( MHT( I ).EQ.3 ) MHP( I ) = LHM + 1
        IF ( MHT( I ).EQ.4 ) MHP( I ) = MHL( I )
        IF ( MHT( I ).EQ.5 ) MHP( I ) = LHM + 1
      ENDDO
C
      CALL DBISVF( INARR, MHT, MHL, MHM, IFLAG, VEC, XARR )
C
      PRINT *,' Enter ROOT for output files.'
      READ ( 5, 80 ) ROOT
C
      DO I = 1, 80
        IF ( ROOT(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 100
        ENDIF
      ENDDO
 100  CONTINUE
      FNAME(1:ILEN) = ROOT(1:ILEN)
      FNAME(ILEN+1:ILEN+5) = '.ints'
      FNAME = FNAME(1:ILEN+5)
C
      CALL HMFWT( NH, MHT, MHL, MHM, MHP, NDCS, MHBCS, MHBCS,
     1                  LU, FNAME )
C
      FNAME(ILEN+1:ILEN+5) = '.xarr'
      FNAME = FNAME(1:ILEN+5)
      CALL XARRWT( NR, XARR, LU, FNAME, IFORM )
C
      FNAME(ILEN+1:ILEN+5) = '.vecs'
      FNAME = FNAME(1:ILEN+5)
      CALL SVFWT( INARR, LU, IFORM, VEC, FNAME )
C
 80   FORMAT(A)
      STOP
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
C subroutine Dynamo Benchmark Initial State Vector Form **************
C            -      -         -       -     -      -    **************
C Steve Gibbons Fri Feb  4 08:36:23 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Fills in vector with the initial state for dynamo benchmark        C
C as suggested by U. Christensen. Values are given by the function   C
C DBPICS.                                                            C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Indexes the array VEC. Dim (3). See INDFUN.        C
C                                                                    C
C     MHT       : Dim (*). MHT(ih) = 1 --> poloidal velocity.        C
C                          MHT(ih) = 2 --> toroidal velocity.        C
C                          MHT(ih) = 3 --> temperature.              C
C                          MHT(ih) = 4 --> poloidal magnetic field   C
C                          MHT(ih) = 5 --> toroidal magnetic field   C
C                                                                    C
C     MHL       : Dim (*). MHL(ih) = l, Spherical harmonic degree.   C
C                                                                    C
C     MHM       : Dim (*). MHM(ih) = m,  for cos m phi dependence.   C
C                          MHM(ih) = -m, for sin m phi dependence.   C
C                                                                    C
C     IFLAG     : 1 --> full system including magnetic field.        C
C                 2 --> velocity and temperature only.               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Dim ( * ). Vector containing initial solution.     C
C                                                                    C
C     XARR      : Dim ( * ). Vector containing grid node radii.      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DBISVF( INARR, MHT, MHL, MHM, IFLAG, VEC, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( 3 ), MHT( * ), MHL( * ), MHM( * ), IFLAG
      DOUBLE PRECISION VEC( * ), XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER LH, NTH, NPH, MMAX, NPHPTS
      PARAMETER ( LH = 4, NTH = 6, NPH = 16 )
      DOUBLE PRECISION SHC( LH*( LH + 2) ), ZCOEF, X1, X2,
     1                 SF( NPH, NTH ), GAUW( NTH ), GAUX( NTH )
      DOUBLE PRECISION PA ( ( LH + 1 )*( LH + 2 )/2 , NTH ),
     1                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTH )
      DOUBLE PRECISION FTF1( 2*NPH ), FTF2( 2*NPH ), FTF3( 2*NPH ),
     1                 VF( 2, NTH, 3)
      DOUBLE PRECISION QST( LH*( LH + 2), 3 ), RI, RO, RAD, LOW,
     1                 THE, PHI, COSTH, DPHI, PI, ZERO, DBPISC, COEF,
     2                 SQRLL1
      PARAMETER ( PI=3.14159265358979312D0, ZERO = 0.0d0,
     1            X1 = -1.0d0, X2 =  1.0d0, LOW = 1.0d-7 )
      INTEGER IR, NR, NH, IOP, ILEN, ITH, IPH, ICOMP, IHARM, NHARM,
     1        L, M, ICS, IH, MM, IND, INDFUN, ILN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C
      RI     = XARR(  1 )
      RO     = XARR( NR )
      NHARM  = LH*(LH + 2)
C
      IF ( DABS( RI ).LT.LOW ) THEN
        ILN = 2
      ELSE
        ILN = 1
      ENDIF
C
C Zero vector
C
      IOP  = 0
      ILEN = NR*NH
      CALL VECOP( VEC, ZERO, ILEN, IOP )
C
C Calculate Gauss points, Legendre functions etc.
C
      CALL GAUWTS ( X1, X2, GAUX, GAUW, NTH )
      CALL SCHNLA ( PA, DPA, GAUX, LH, NTH )
C
C Do initial temperature field
C
      ICOMP  = 4
      MMAX   = 4
      NPHPTS = NPH
      DPHI   = 2.0d0*PI/DBLE( NPHPTS )
      DO  IR = ILN, NR
        RAD  = XARR( IR )
C       .
C       . Fill scalar function array.
C       .
        DO ITH = 1, NTH
          COSTH = GAUX( ITH )
          THE   = ACOS( COSTH )
          DO IPH = 1, NPHPTS
            PHI  = DBLE( IPH - 1 )*DPHI
            SF( IPH, ITH ) = DBPISC( ICOMP, RAD, THE, PHI, RI, RO )
          ENDDO
        ENDDO
C       .
C       . SF is filled. Convert to spectral coefficients
C       .
        CALL FORSST ( SHC, SF, GAUW, PA, FTF1, LH, MMAX,
     1                NTH, NPHPTS, ZCOEF )
C       .
        DO IHARM = 1, NHARM
          COEF   = SHC( IHARM )
          IF ( DABS( COEF ).GT.LOW ) THEN
            CALL LMFIND( IHARM, L, M, ICS )
            IF ( ICS.EQ.1 ) MM = M
            IF ( ICS.EQ.2 ) MM = -M
            DO IH = 1, NH
              IF ( MHT( IH ).EQ.3 .AND. MHL( IH ).EQ.L
     1             .AND. MHM( IH ).EQ.MM ) THEN
                IND = INDFUN( IR, IH, INARR )
                VEC( IND ) = COEF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
C       .
      ENDDO
C
      IF ( IFLAG.EQ.2 ) RETURN
C
C Do initial magnetic field.
C
      MMAX   = 0
      NPHPTS = 2
      DO  IR = ILN, NR
        RAD  = XARR( IR )
C       .
C       . Fill vector function array.
C       .
        DO ITH = 1, NTH
          COSTH = GAUX( ITH )
          THE   = ACOS( COSTH )
          DO ICOMP = 1, 3
            VF( 1, ITH, ICOMP ) =
     1                DBPISC( ICOMP, RAD, THE, PHI, RI, RO )
            VF( 2, ITH, ICOMP ) = VF( 1, ITH, ICOMP )
          ENDDO
        ENDDO
C       .
C       . VF is filled. Convert to spectral coefficients
C       .
        CALL VF2QST ( QST, VF, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1                FTF3, ZCOEF, LH, NTH, NPHPTS, MMAX )
C       .
        DO IHARM = 1, NHARM
C         .
C         . First do poloidal parts
C         .
          COEF   = QST( IHARM, 1 )
          IF ( DABS( COEF ).GT.LOW ) THEN
            CALL LMFIND( IHARM, L, M, ICS )
            IF ( ICS.EQ.1 ) MM = M
            IF ( ICS.EQ.2 ) MM = -M
            DO IH = 1, NH
              IF ( MHT( IH ).EQ.4 .AND. MHL( IH ).EQ.L
     1             .AND. MHM( IH ).EQ.MM ) THEN
                IND = INDFUN( IR, IH, INARR )
                VEC( IND ) = COEF*RAD/DBLE( L*L + L )
              ENDIF
            ENDDO
          ENDIF
C         .
C         . Now do toroidal parts
C         .
          COEF   = QST( IHARM, 3 )
          IF ( DABS( COEF ).GT.LOW ) THEN
            CALL LMFIND( IHARM, L, M, ICS )
            IF ( ICS.EQ.1 ) MM = M
            IF ( ICS.EQ.2 ) MM = -M
            DO IH = 1, NH
              IF ( MHT( IH ).EQ.5 .AND. MHL( IH ).EQ.L
     1             .AND. MHM( IH ).EQ.MM ) THEN
                IND = INDFUN( IR, IH, INARR )
                VEC( IND ) = (-1.0d0)*COEF/SQRLL1( L )
              ENDIF
            ENDDO
          ENDIF
C         .
        ENDDO
C       .
      ENDDO
C
      IF ( IFLAG.EQ.1 ) RETURN
C
      PRINT *,' Subroutine DBPSVF.'
      PRINT *,' IFLAG is neither 1 nor 2.'
      PRINT *,' Program aborted.'
      STOP
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
C subroutine SCHmidt Normalised Legendre function Array **************
C            ---     -          -                 -     **************
C Steve Gibbons 22.4.97                                              C
C____________________________________________________________________C
C Does the same as SCHNLF except that instead of a single valued X   C
C for one theta point, it fills arrays PA and DPA with the           C
C legendre Functions etc. for each of the NTHPTS values of cos(theta)C
C in the array GAUX.                                                 C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH	: Highest degree, l, of spherical harmonic.          C
C     NTHPTS	: Number of theta points.                            C
C  Double Precision                                                  C
C  ----------------                                                  C
C     PA	: Schmidt Normalised Legendre Functions Dimension.   C
C		   {  ( LH + 1 )*( LH + 2 )/2 , NTHPTS }             C
C		   P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA	: Derivatives of the above.                          C
C     GAUX	: Array of cosines to the NTHPTS angles.             C
C                  Dimension ( NTHPTS ).                             C
C____________________________________________________________________C
C Functions ... Calling Proceedures :-                               C
C  Double Precision                                                  C
C  ----------------                                                  C
C PMM ( M, S )					                     C
C DPMM ( M, C, S )				                     C
C PMM1 ( M, X, PMM0 )                                                C
C PLM ( L, M, X, PLMIN1, PLMIN2 )				     C
C DPMM1 ( M , X , S, PMM , DPMM)                                     C
C DPLM ( L, M , X , S, PMM1 , DPMM1, DPMM2 )                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SCHNLA ( PA, DPA, GAUX, LH, NTHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS
      DOUBLE PRECISION GAUX( NTHPTS ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER INDEX,L,M,IOLD1,IOLD2,NTHETA
      DOUBLE PRECISION SINE,PMIN1,PMIN2,TOL,DPMIN1,
     1                 DPMIN2,X
      PARAMETER (TOL=1.0d-6)
C____________________________________________________________________C
C Variable declarations - Functions called ..........................C
      DOUBLE PRECISION PMM,PMM1,PLM,DPMM,DPMM1,DPLM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check validity of arguments ....
C
C ........................... now loop around theta points
      DO NTHETA = 1, NTHPTS
C
        X = GAUX ( NTHETA )
        SINE = X*X
        IF ( SINE.GT.1.0d0 ) THEN
           PRINT *,' Subroutine SCHNLA.'
           PRINT *,' Illegal Cos(theta) has been entered.'
           PRINT *,' ( For NTHETA = ',NTHETA,' )'
           PRINT *,' Program stopped.'
           STOP
        ENDIF
C
C Set SINE (theta) in terms of X
        SINE = DSQRT ( (1.0d0 + X)*(1.0d0 - X) )
        IF ( SINE.LT.TOL ) THEN
           PRINT *,' Subroutine SCHNLA.'
           PRINT *,' SINE is too small. Division by zero imminent.'
           PRINT *,' ( For NTHETA = ',NTHETA,' )'
           PRINT *,' Program stopped.'
           STOP
        ENDIF
C..................... first calculate the P_l^m ..........
        DO M = 0, LH - 2
C                        ............. Calculate P_M^M .....
           L = M
           INDEX = L*(L+1)/2+M+1
           PA ( INDEX , NTHETA) = PMM ( M , SINE )
           DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
C                         ............. Calculate P_(M+1)^M .
           PMIN1 = PA ( INDEX , NTHETA)
           DPMIN1 = DPA ( INDEX , NTHETA)
           IOLD2 = INDEX
           L = L + 1
           INDEX = L*(L+1)/2+M+1
           PA (INDEX , NTHETA) = PMM1 ( M , X , PMIN1 )
           DPA (INDEX , NTHETA) = DPMM1 (M,X,SINE,PMIN1, DPMIN1)
           IOLD1 = INDEX
C                         ......... Calculate P_L^M general .
           DO L = M + 2, LH
              PMIN2 = PA ( IOLD2 , NTHETA)
              PMIN1 = PA ( IOLD1 , NTHETA)
              DPMIN2 = DPA ( IOLD2 , NTHETA)
              DPMIN1 = DPA ( IOLD1 , NTHETA)
              INDEX = L*(L+1)/2+M+1
              PA ( INDEX , NTHETA) = PLM ( L,M,X,PMIN1,PMIN2 )
              DPA ( INDEX , NTHETA) = DPLM (L,M,X,SINE , PMIN1, 
     1                              DPMIN1, DPMIN2 )
              IOLD2 = IOLD1
              IOLD1 = INDEX
           ENDDO
        ENDDO
        M = LH - 1
        L = M
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM ( M , SINE )
        DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
        PMIN1 = PA( INDEX , NTHETA)
        DPMIN1 = DPA( INDEX , NTHETA)
        L = LH
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM1 ( M , X , PMIN1 )
        DPA(INDEX ,NTHETA) = DPMM1 (M ,X ,SINE,PMIN1,DPMIN1 )
        M = LH
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM ( M , SINE )
        DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
C......................finished calculating P_l^m .........
      ENDDO
C......................finished looping around theta points

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
C*********************************************************************
C subroutine VECtor OPeration ****************************************
C Steve Gibbons 22.4.97 Fills vector with a constant, multiplies a   C
C                       vector by a constant or adds a constant.     C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IOP	: Type of operation required.                        C
C                  IOP=0  -->  Each element of the vector = CONST    C
C                  IOP=1  -->  Each el. is multiplied by CONST       C
C                  IOP=2  -->  Each el. is added to CONST            C
C     N		: Length of the vector.                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VEC	: Vector - dimension ( N )                           C
C     CONST     : Double precision constant.                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VECOP ( VEC, CONST, N, IOP )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N, IOP
      DOUBLE PRECISION VEC( N ), CONST
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First do the case of making an element constant.
      IF ( IOP.EQ.0 ) THEN
         DO I = 1, N
            VEC ( I ) = CONST
         ENDDO
         RETURN
      ENDIF
C Now do multiplying a vector
      IF ( IOP.EQ.1 ) THEN
         DO I = 1, N
            VEC ( I ) = VEC( I )*CONST
         ENDDO
         RETURN
      ENDIF
C Now do adding a vector
      IF ( IOP.EQ.2 ) THEN
         DO I = 1, N
            VEC ( I ) = VEC( I ) + CONST
         ENDDO
         RETURN
      ENDIF
C____________________________________________________________________C

      PRINT *,' Subroutine VECOP. IOP must be 0, 1 or 2.'
      PRINT *,'Program aborted.'
      STOP
      END
C*********************************************************************


C*********************************************************************
C subroutine GAUWTS **************************************************
C Adapted 22.4.97 from Numerical Recipes routine GAULEG              C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHPTS	: Number of theta points.                            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X1	: Starting value for integration.                    C
C     X2	: Ending value for integration.                      C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX	: Array containing abscissae of the Gauss-Legendre   C
C                  NTHPTS-points quadrature formula.                 C
C     GAUW      : Array containing the weights for the above points. C
C                                                                    C
C ( Both GAUX and GAUW have dimension NTHPTS ).                      C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE GAUWTS ( X1, X2, GAUX, GAUW, NTHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHPTS
      DOUBLE PRECISION X1, X2, GAUX( NTHPTS ), GAUW( NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER M,I,J
      DOUBLE PRECISION XM,XL,P1,P2,P3,EPS,PP,Z,Z1
      PARAMETER (EPS=1.0d-13)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C 
      IF ( ABS( X1 - X2 ).LT.EPS ) THEN
        PRINT *,' Subroutine GAUWTS,'
        PRINT *,' X1 = ', X1
        PRINT *,' X2 = ', X2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C ....... the roots are symmetric in the interval so need only find
C half of them.
      M = ( NTHPTS + 1 )/2
      XM = 0.5d0 * ( X2 + X1 )
      XL = 0.5d0 * ( X2 - X1 )
C ........start looping over the desired roots .......................
      DO I = 1, M
         Z = DCOS( PI*( DBLE(I) - 0.25d0)/( DBLE(NTHPTS) + 0.5d0 ))
C           ..... starting with this approximation to the Ith root, we
C                enter the main loop of refinement by Newton's method.
 100     CONTINUE
            P1 = 1.0D0
            P2 = 0.0D0
C           ........... Loop up the recurrence relation to get the
C                      legendre Polynomial evaluated at Z.
            DO J = 1, NTHPTS
               P3 = P2
               P2 = P1
               P1 = ((2.0d0*J-1.0d0)*Z*P2 - (J-1.0d0)*P3)/DBLE( J )
            ENDDO
C           ..................... finish recurrence relation loop ...
C ... P1 is now the desired Legendre Polynomial. We now compute PP,
C    its derivative by a standard relation involving also P2, the 
C    polynomial of one order lower.
            PP = NTHPTS*(Z*P1-P2)/(Z*Z-1.0d0)
            Z1 = Z
            Z = Z1 - P1/PP
         IF ( ABS(Z-Z1).GT.EPS ) GOTO 100
C ...........scale the root to the desired interval .................
         GAUX( I ) = XM - XL*Z
C ...........and add its symmetric counterpart ......................
         GAUX( NTHPTS+1-I ) = XM + XL*Z
C ...........calculate the weight ...................................
         GAUW( I ) = 2.0d0*XL/((1.0d0-Z*Z)*PP*PP)
C ...........and add its symmetric counterpart ......................
         GAUW( NTHPTS + 1 - I ) = GAUW( I )
      ENDDO
C ......... end looping over the desired roots .......................
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
C*********************************************************************
C subroutine Vector Function TO QST coefficients *********************
C            -      -        -- ---              *********************
C Steve Gibbons 26.4.97                                              C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum spherical harmonic order, m.               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     GAUW      : Corresponding weights. As above.                   C
C     VF        : Vector Function. An array of dimensions            C
C                  ( NPHPTS, NTHPTS, 3) which contain the R, THETA   C
C                  and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IPHI, ITHETA, 2 ) is the Theta      C
C                  compontent of the vector at (iphi, itheta).       C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     QST       : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector c.f. eqn (38).            C
C                  Has dimensions (  LH*(LH+2) , 3).                 C
C               QST (l*l+2m,1) = q_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,1) = q_l^mc(r_i)                         C
C               QST (l*l+2m,2) = s_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,2) = s_l^mc(r_i)                         C
C               QST (l*l+2m,3) = t_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,3) = t_l^mc(r_i)                         C
C                                                                    C
C     ZCOEF     : Coefficient of the monopole term in the R compon.  C
C                 This is not necessarily zero in general.           C
C                                                                    C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VF2QST ( QST, VF, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1                    FTF3, ZCOEF, LH, NTHPTS, NPHPTS, MMAX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS, NPHPTS, MMAX
      DOUBLE PRECISION QST(  LH*(LH+2) , 3), ZCOEF,
     1                 VF( NPHPTS, NTHPTS, 3),
     2                 GAUX( NTHPTS ), GAUW( NTHPTS )
      DOUBLE PRECISION FTF1( 2*NPHPTS ),
     1                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ), 
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, ISIGN, ITHREE, ITHETA, 
     1        INDCOS, INDSIN, IP, IND1, L, M, IPHI
      DOUBLE PRECISION ZERO,X,SINE,TERM, WEIGHT, W1, W2
      PARAMETER ( ZERO = 0.0d0 , ITHREE=3 )
C____________________________________________________________________C
C Functions used :-
      DOUBLE PRECISION SQRLL1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the validity of arguments .....
C No need to check NPHPTS is power of 2 - this is done
C by FFTRLV ...
C ......... need to have MMAX < NPHPTS/2
C
      IF ( NTHPTS.LE.LH ) THEN
         PRINT *,' Subroutine VF2QST.'
         PRINT *,' NTHPTS must be atleast ', LH + 1
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPTS.LE.(2*MMAX+1) ) THEN
         PRINT *,' Subroutine VF2QST.'
         PRINT *,' NPHPTS must be atleast ',(2*MMAX+1)
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
C_____________________________________________________________________
C
C ............ set array QST to zero .................................
      IOP = 0
      IND1 = LH*(LH+2)
      CALL MATOP ( QST, ZERO, IND1, ITHREE, IOP )
      ZCOEF = ZERO
C ....................................................................
C ...... now start to loop around theta points .......................
      DO ITHETA = 1, NTHPTS
C
C******** Evaluate lengendre functions at given theta ****************
C
         X = GAUX( ITHETA )
         SINE = DSIN ( ACOS ( X ) )
         WEIGHT = GAUW( ITHETA )
C        .
C        . Zero the arrays ftf1, ftf2 and ftf3
C        .
         IOP = 0
         IND1 = 2*NPHPTS
         CALL DVECZ( FTF1, IND1 )
         CALL DVECZ( FTF2, IND1 )
         CALL DVECZ( FTF3, IND1 )
C
C .................. firstly enter the VF information into FTF1, FTF2,
C .................. FTF3 for r, theta, phi components respectively of
C .................. VF.
C
         DO IPHI = 1, NPHPTS
            FTF1( 2*IPHI - 1 ) = VF( IPHI, ITHETA, 1 )
            FTF2( 2*IPHI - 1 ) = VF( IPHI, ITHETA, 2 )
            FTF3( 2*IPHI - 1 ) = VF( IPHI, ITHETA, 3 )
            ZCOEF = ZCOEF + WEIGHT*VF( IPHI, ITHETA, 1 )/4.0d0
         ENDDO
C .................. now perform Forward Discreet Fourier Transforms
C .................. on FTF1, FTF2, FTF3
         ISIGN = 1
         CALL FFTRLV ( FTF1, NPHPTS, ISIGN )
         CALL FFTRLV ( FTF2, NPHPTS, ISIGN )
         CALL FFTRLV ( FTF3, NPHPTS, ISIGN )
C ...................................................................
C .................. Now let's loop around the Harmonics. ...........
C ...................................................................
         DO L = 1, LH
            W1 = WEIGHT*(2.0d0*L + 1.0d0)/4.0d0
            W2 = W1/SQRLL1( L )
            INDCOS = L*L
C .................. Now let's loop around the order, but considering
C .................. the Case M = 0 separately. 
C ___ CASE M = 0 ___ ************************************************
            M = 0
            IP = L*(L+1)/2+M+1
            QST( INDCOS , 1 ) = QST( INDCOS , 1) +
     1                          W1*PA( IP ,ITHETA )*FTF1( 1 )
            QST( INDCOS , 2 ) = QST( INDCOS , 2) +
     1                          W2*DPA( IP ,ITHETA )*FTF2( 1 )
            QST( INDCOS , 3 ) = QST( INDCOS , 3) +
     1                          W2*DPA( IP ,ITHETA )*FTF3( 1 )
C
            INDCOS = INDCOS - 1
C
C ___ Looping around from 1 to L ___ ********************************
C
C
            DO M = 1, MIN( MMAX, L )
               IP = IP + 1
               TERM = W2*DBLE(M)*PA( IP ,ITHETA )/SINE
               INDCOS = INDCOS + 2
               INDSIN = INDCOS + 1
C                                                   ... eqn 194
C                                                       Qlm cos
C***.Equation for Qlm COS ***************************
               QST( INDCOS , 1 ) = QST( INDCOS , 1) +
     1          W1*PA( IP ,ITHETA )*FTF1( 2*M + 1)
C ..............contribution from B(rad)cos .........
C                                                   ... eqn 195
C                                                       Qlm sin
C***.Equation for Qlm SIN ***************************
               QST( INDSIN , 1 ) = QST( INDSIN , 1) +
     1          W1*PA( IP ,ITHETA )*FTF1( 2*M + 2)
C ..............contribution from B(rad)sin .........
C                                                   ... eqn 203
C                                                       Slm cos
C***.Equation for Slm COS ***************************
               QST( INDCOS , 2 ) = QST( INDCOS , 2) +
     1          W2*DPA( IP ,ITHETA )*FTF2( 2*M + 1 ) -
     2          TERM*FTF3( 2*M + 2 )
C ..............contribution from B(theta)cos .......
C ..............contribution from B(phi)sin .........
C                                                   ... eqn 204
C                                                       Slm sin
C***.Equation for Slm SIN ***************************
               QST( INDSIN , 2 ) = QST( INDSIN , 2) +
     1          W2*DPA( IP ,ITHETA )*FTF2( 2*M + 2 ) +
     2          TERM*FTF3( 2*M + 1 )
C ..............contribution from B(theta)sin .......
C ..............contribution from B(phi)cos .........
C                                                   ... eqn 206
C                                                       Tlm cos
C***.Equation for Tlm COS ***************************
               QST( INDCOS , 3 ) = QST( INDCOS , 3) +
     1          W2*DPA( IP ,ITHETA )*FTF3( 2*M + 1 ) +
     2          TERM*FTF2( 2*M + 2 )
C ..............contribution from B(phi)cos .........
C ..............contribution from B(theta)sin .......
C                                                   ... eqn 207
C                                                       Tlm sin
C***.Equation for Tlm SIN ***************************
               QST( INDSIN , 3 ) = QST( INDSIN , 3) -
     1          TERM*FTF2( 2*M + 1 ) +
     2          W2*DPA( IP ,ITHETA )*FTF3( 2*M + 2 )
C ..............contribution from B(theta)cos .......
C ..............contribution from B(phi)sin .........
C
            ENDDO
C
C
C Done looping around from 1 to L ___ *******************************
C
C
         ENDDO
C ...................................................................
C .................. Ended looping around the harmonics .............
C ...................................................................
      ENDDO
C ...... ended looping around theta points ...........................
C ....................................................................
      ZCOEF = 2.0d0*ZCOEF/DBLE(NPHPTS)
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine FORward Scalar Spherical Transform (array version) ******
C            ---     -      -         -                              C
C Steve Gibbons 23.4.97                                              C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SF	: Scalar Function. This is an array containing a     C
C                  function over a set of theta and phi points.      C
C                    Dimensions are                                  C
C                      ( NPHPTS , NTHPTS )                           C
C                                                                    C
C                 SF( i, j ) contains the value a scalar function    C
C                 at theta = acos( x_j ) where x_j is the j^{th}     C
C                 element of GAUX returned by GAUWTS and             C
C                 phi = 2 pi ( i - 1 ) / nphpts                      C
C                                                                    C
C     GAUW      : Gauss weights as calculated by GAUWTS. ( NTHPTS )  C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     ZCOEF     : Coefficient of the l=0 spherical harmonic.         C
C                                                                    C
C  Integer                                                           C
C  -------							     C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum spherical harmonic order, m.               C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SHC	: Spherical Harmonic Coefficients.                   C
C                  Has the usual ordering.                           C
C                   SHC( l*l ) = P_l^0                               C
C                   SHC( l*l + 2*m - 1 ) = P_l^mc ( m non zero )     C
C                   SHC( l*l + 2*m ) = P_l^ms ( m non zero )         C
C                  Dimension is ( LH*( LH + 2) )                     C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF       : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FORSST ( SHC, SF, GAUW, PA, FTF, LH, MMAX,
     1                    NTHPTS, NPHPTS, ZCOEF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          LH, NTHPTS, NPHPTS, MMAX
      DOUBLE PRECISION SHC( LH*( LH + 2) ),
     1                 SF( NPHPTS , NTHPTS ), GAUW(NTHPTS)
      DOUBLE PRECISION FTF( 2*NPHPTS ), ZCOEF,
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IOP,ITHETA,ILEV,INDCOS,INDSIN,LENV,ISIGN,I,M, IP
      DOUBLE PRECISION ZERO,WEIGHT
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C    
C Check the validity of arguments .....
C (Don't bother checking NPHPTS is a power of 2
C as this will be done by FFTRLV)
C
      IF ( NTHPTS.LE.LH ) THEN
         PRINT *,' Subroutine FORSST.'
         PRINT *,' NTHPTS must be atleast ', LH + 1
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPTS.LE.(2*MMAX+1) ) THEN
         PRINT *,' Subroutine FORSST.'
         PRINT *,' NPHPTS must be atleast ',2*MMAX+1
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C_____________________________________________________________________
C
      ISIGN = 1
C ... zero the array of spherical harm coeffs --> IOP = 0
      IOP = 0
      LENV = LH*(LH+2)
      CALL VECOP ( SHC, ZERO, LENV, IOP )
      ZCOEF = 0.0d0
 
C ............... begin looping around theta points ..................
      DO ITHETA = 1, NTHPTS
         LENV = 2*NPHPTS
         CALL VECOP ( FTF, ZERO, LENV, IOP )
         WEIGHT = GAUW( ITHETA )
C ................ now (considering equation we have w_j,
C all of the P_l^m(cos(theta_j) ) so now need the Fourier transform
C of the SF for (theta_j) ......
         DO I = 1, NPHPTS
            FTF ( 2*I - 1) = SF( I, ITHETA )
            ZCOEF = ZCOEF + WEIGHT * SF( I, ITHETA )/DBLE(NPHPTS)
         ENDDO
C ................ o.k. now let's Fourier Transform .................
C Array ftf currently contains value of function at phi point
C (i) in element ftf( 2i - 1 ) ...
C
         CALL FFTRLV ( FTF, NPHPTS, ISIGN )
C
C ftf now contains c_m in ftf(2m+1) and s_m in ftf(2m+2)
C ................ let's loop around the harmonics ..................
         DO ILEV = 1, LH         
            INDCOS = ILEV*ILEV
C ................................... first let's do the M = 0 case .
            M = 0
            IP = ILEV*(ILEV+1)/2+M+1
            SHC ( INDCOS ) = SHC( INDCOS ) + PA( IP , ITHETA)*
     1                        WEIGHT*FTF( 1 )
            INDCOS = INDCOS + 1
            INDSIN = INDCOS + 1
C ............ now let's loop around M from 1 to ILEV ..............
            DO M = 1, MIN( ILEV, MMAX )
               IP = ILEV*(ILEV+1)/2+M+1
               SHC(INDCOS) = SHC(INDCOS) + PA(IP, ITHETA)*
     1                        WEIGHT*FTF(2*M+1)
               SHC(INDSIN) = SHC(INDSIN) + PA(IP, ITHETA)*
     1                        WEIGHT*FTF(2*M+2)
               INDCOS = INDCOS + 2
               INDSIN = INDSIN + 2
            ENDDO
         ENDDO
C .................. end loop around the harmonics ..................
      ENDDO
C ................. end looping around theta points ..................
C Now just loop around the SHC's to multiply by the schmidt normalisation
C and scale factor of fft.
      ZCOEF = ZCOEF / 2.0d0
      DO ILEV = 1, LH
         WEIGHT = DBLE( 2*ILEV + 1 )/4.0d0
         DO M = 0 , 2*MIN( ILEV, MMAX )
            SHC( ILEV*ILEV + M ) = SHC ( ILEV*ILEV + M ) * WEIGHT
         ENDDO
      ENDDO
      RETURN
      END
C*********************************************************************
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
C d. p. function Dynamo Benchmark Project Initial State Components ***
C                -      -         -       -       -     -          ***
C Steve Gibbons Thu Feb  3 17:52:23 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Supplies values of components given for initial state for          C
C dynamo benchmark as suggested by U. Christensen.                   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ICOMP     : Selects the component to be calculated.            C
C                                                                    C
C                 (1) --> B_r (radial magnetic field)                C
C                 (2) --> B_theta (theta magnetic field)             C
C                 (3) --> B_phi (phi magnetic field)                 C
C                 (4) --> Theta scalar function.                     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Distance from centre of sphere.                    C
C     THE       : Co-latitude in radians.                            C
C     PHI       : Longitude in radians.                              C
C     RI        : Radius of inner boundary                           C
C     RO        : Radius of outer boundary                           C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DBPISC( ICOMP, RAD, THE, PHI, RI, RO )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION DBPISC, RAD, THE, PHI, RI, RO
      INTEGER ICOMP
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION PI, LOW, COSTH, RAD3, RI4, PIRMRI, TWOTHE,
     1                 X, X2, X4, X6, SINTH, SINTH4, FOURPH, BRAC,
     2                 A, FRAC
      PARAMETER ( PI = 3.14159265358979312D0, LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( RAD.LT.LOW ) THEN
        PRINT *,' Function DBPISC.'
        PRINT *,' RAD = ', RAD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ICOMP.NE.1 .AND. ICOMP.NE.2 .AND.
     1     ICOMP.NE.3 .AND. ICOMP.NE.4 ) THEN
        PRINT *,' Function DBPISC.'
        PRINT *,' ICOMP = ', ICOMP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Case of B_r
C
      IF ( ICOMP.EQ.1 ) THEN
        RI4    = RI*RI*RI*RI
        RAD3   = RAD*RAD*RAD
        COSTH  = COS( THE )
        BRAC   = 6.0d0*RAD - 8.0d0*RO + 2.0d0*RI4/RAD3
        DBPISC = -5.0d0*BRAC*COSTH/8.0d0
        RETURN
      ENDIF
C
C Case of B_theta
C
      IF ( ICOMP.EQ.2 ) THEN
        RI4    = RI*RI*RI*RI
        RAD3   = RAD*RAD*RAD
        SINTH  = SIN( THE )
        BRAC   = -9.0d0*RAD + 8.0d0*RO + RI4/RAD3
        DBPISC = -5.0d0*BRAC*SINTH/8.0d0
        RETURN
      ENDIF
C
C Case of B_phi
C
      IF ( ICOMP.EQ.3 ) THEN
        PIRMRI = PI*( RAD - RI )
        TWOTHE = 2.0d0*THE
        DBPISC = 5.0d0*SIN( PIRMRI )*SIN( TWOTHE )
        RETURN
      ENDIF
C
C Case of Theta
C
      IF ( ICOMP.EQ.4 ) THEN
        X      = 2.0d0*RAD - RI - RO
        X2     = X*X
        X4     = X2*X2
        X6     = X4*X2
        BRAC   = ( 1.0d0 - 3.0d0*X2 + 3.0d0*X4 - X6 )
        SINTH  = SIN( THE )
        SINTH4 = SINTH*SINTH*SINTH*SINTH
        FOURPH = 4.0d0*PHI
        A      = 0.1d0
        FRAC   = 210.0d0/SQRT( 17920.0d0*PI )
        DBPISC = A*FRAC*BRAC*SINTH4*COS( FOURPH )
        RETURN
      ENDIF
C
      END
C*********************************************************************
C*********************************************************************
C function SQuare Root of L*(L+1) ************************************
C          --     -       -  - -  ************************************
C Steve Gibbons 25.4.97                                              C
C____________________________________________________________________C
C
      FUNCTION SQRLL1 ( L )
      IMPLICIT NONE
      DOUBLE PRECISION SQRLL1,Q
      INTEGER L
   
      IF ( L.LT.0 ) THEN
         PRINT *,' Function SQRLL1. L less than 0. Program aborted.'
         STOP
      ENDIF
      Q = DBLE( L )
      SQRLL1 = DSQRT( Q*Q + Q )
      RETURN
      END
C*********************************************************************
C********************************************************************
C subroutine L and M FIND *******************************************
C            -     - ---- *******************************************
C Steve Gibbons 12.4.97                                             C
C Modified 12.6.97 to have IT = 1 for cos and 2 for sine            C
C___________________________________________________________________C
C Given a number of harmonic, N; LMFIND will return the level, L,   C
C the order, M, and IT - which is equal to 2 for sin                C
C spherical harmonics and 1 for cosine ones.                        C
C All the above are integers - no point in a variable list .....    C
C     N = L*L for M = 0, IT = 1                                     C
C     N = L*L + 2*M - 1 for non-zero M and IT = 1                   C
C     N = L*L + 2*M for non-zero M and IT = 2                       C
C If N = 0, we have a monopole: L = 0, M = 0 and IT = 1.            C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE LMFIND ( N, L, M, IT)
      IMPLICIT NONE
C___________________________________________________________________C
C Variable Declarations - Parameters ...............................C
      INTEGER N,L,M,IT
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      INTEGER LL,IDIFF,N2,ITWO
      PARAMETER (ITWO=2)
C___________________________________________________________________C
C First put N into N2 so that N is not altered
      N2=N
      IF ( N2.EQ.0 ) THEN
         L=0
         M=0
         IT=1
         RETURN
      ENDIF
      IF ( N2.LT.1 ) THEN 
         WRITE (6,989)
         STOP
 989  FORMAT (' Subroutine LMFIND. N<1 - Program stopped.')
      ENDIF
      IF ( N2.EQ.1 ) THEN
         L=1
         M=0
         IT=1
         RETURN
      ENDIF
      L=1
 500  CONTINUE
      LL=L*L
      IDIFF = N2 - LL
      IF ( IDIFF.GT.0 ) THEN
         L=L+1
         GOTO 500
      ENDIF
      IF ( IDIFF.EQ.0 ) THEN
         M=0
         IT=1
         RETURN
      ENDIF
      L=L-1
      LL=L*L
      N2=N2-LL
C so we know now that N2 is equal to either 2*M or 2*M-1
C corresponding to IT=1 and IT=2 respectively
      IDIFF = MOD ( N2, ITWO )
      IF ( IDIFF.EQ.1) THEN
         IT = 1
         M = (N2+1)/ITWO
      ELSE
         IT = 2
         M = N2/ITWO
      ENDIF
      RETURN
      END
C********************************************************************
C*********************************************************************
C function DPLM *******************************************************
C Steve Gibbons 18.4.97                                              C
C____________________________________________________________________C
C                                                                    C
C Calculates general P_l^m derivative from eqn 187 in my notes       C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C     L         : Well it's L isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     S         : Sin(theta)                                         C
C     PLMIN1    : P_(l-1)^m ( X )                                    C
C     DPLMN1   : P_(l-1)^m ( X ) derivative                         C
C     DPLMN2   : P_(l-2)^m ( X ) derivative                         C
C____________________________________________________________________C
C Local Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     RL        : DBLE  ( L )                                        C
C     RM        : DBLE  ( M )                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPLM ( L, M, X, S, PLMIN1, DPLMN1, DPLMN2 )
      IMPLICIT NONE
      DOUBLE PRECISION DPLM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M,L
      DOUBLE PRECISION X, S, PLMIN1, DPLMN1, DPLMN2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM,RL
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function DPLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function DPLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' DPLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = DBLE( M )
      RL = DBLE( L )
      DPLM = ( 2.0d0*RL - 1.0d0 )*(X*DPLMN1-S*PLMIN1)
      DPLM = DPLM - DPLMN2*DSQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      DPLM = DPLM/DSQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function PMM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Gives the Schmidt Normalised Legendre Function P_m^m ( X )         C
C from eqn. 175 in my notes ie                                       C
C                                                                    C
C                   ( 2m - 1)!! (1- XX)^(m/2) * SQRT (2.0d0 )        C
C   P_m^m( X ) =  ---------------------------------------------      C
C                        SQRT ( (2m)! )                              C
C                                                                    C
C       for m non-zero and                                           C
C                                                                    C
C   P_0^0( X ) = 1.0d0                                               C
C                                                                    C
C N.B. The double factorial sign means the product of all ODD        C
C integers between 1 and ( 2m - 1 ).                                 C
C Best to use the form in eq 184 to calculate PMM.                   C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SINE      : Sin (theta)                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMM ( M, S )
      IMPLICIT NONE
      DOUBLE PRECISION PMM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION S
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
      DOUBLE PRECISION TWOI
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM. M < 0 error. Stopped.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         PMM = 1.0d0
         RETURN
      ENDIF
C ................................. so M is greater than 0
      PMM = DSQRT ( 2.0d0 )
      DO I = 1, M
         TWOI = 2.0d0*DBLE(I)
         PMM = PMM * S * DSQRT ( (TWOI - 1.0d0 )/TWOI )
      ENDDO
      RETURN
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
C function DPMM ******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Calculates the derivative of P_m^m(theta) according to equation    C
C 185 in my notes i.e.                                               C
C                                                                    C
C d P_m^m(theta)/ d(theta) = sqrt(2.0)*M*cos(theta)/sin(theta) * A   C
C  with A = product_(i=1)^m ( SQRT( (2i-1)/2i ) * sin(theta) )       C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SINE      : Sin (theta)                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPMM ( M , C , S )
      IMPLICIT NONE
      DOUBLE PRECISION DPMM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION C, S
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
      DOUBLE PRECISION TWOI
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function DPMM. M < 0 error. Stopped.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         DPMM = 0.0d0
         RETURN
      ENDIF
C ................................. so M is greater than 0
      DPMM = DSQRT ( 2.0d0 )*DBLE(M)*C/S
      DO I = 1, M
         TWOI = 2.0d0*DBLE(I)
         DPMM = DPMM * S * DSQRT ( (TWOI - 1.0d0 )/TWOI )
      ENDDO
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine MATrix OPeration ****************************************
C Steve Gibbons 23.4.97 Does operation on a two-dimensional array.   C
C                                                Can set equal to a  C
C                       constant; multiply by a constant or have a   C
C                       constant added to it.                        C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IOP	: Type of operation to be done.                      C
C                    WHOLE MATRIX OPERATIONS                         C
C                  IOP=0  -->  Each element of the matrix = CONST    C
C                  IOP=1  -->  Each el. is multiplied by CONST       C
C                  IOP=2  -->  Each el. is added to CONST            C
C     NDIM1     : First dimension of the matrix.	             C
C     NDIM2     : Second dimension of the matrix.	             C
C  Double Precision                                                  C
C  ----------------                                                  C
C     MAT	: Matrix with dimension ( NDIM1, NDIM2 )             C
C     CONST     : Double precision constant.                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MATOP ( MAT, CONST, NDIM1, NDIM2, IOP)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM1, NDIM2, IOP
      DOUBLE PRECISION MAT ( NDIM1, NDIM2 ), CONST
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, J
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First do case IOP=0 
      IF ( IOP.EQ.0 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               MAT ( I, J) = CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
C Now do case IOP=1
      IF ( IOP.EQ.1 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               MAT ( I, J) = MAT ( I, J)*CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
C Now do case IOP=2
      IF ( IOP.EQ.2 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               MAT ( I, J) = MAT ( I, J) + CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
      PRINT *,' Subroutine MATOP. IOP must be 0,1 or 2.'
      STOP
      END
C*********************************************************************
C*********************************************************************
C function PMM1 ******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Evaluates the Schmidt Normalised Legendre Function P_(m+1)^m (X)   C
C according to equation 179 in my notes ; i.e.                       C
C                                                                    C
C    P_(m+1)^m (X) = SQRT( 2m+1 ).X.P^m_m(X)                         C
C                                                                    C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     PMM0      : P_m^m(X) as evaluated by Function PMM.             C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMM1 ( M, X, PMM0 )
      IMPLICIT NONE
      DOUBLE PRECISION PMM1
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION X, PMM0
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = DBLE( M )
      PMM1 = X*PMM0*DSQRT( 2.0d0*RM+1.0d0 )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Fast Fourier Transform ReaL Version *********************
C            -    -       -         -  - -       *********************
C Steve Gibbons 22.4.97 (Based on FOUR1 from Numerical Recipes       C
C      I've made slight alterations for double precision, error      C
C      checking at the start and rescaling the coefficients after    C
C      the forward transform.)                                       C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NN	: Number of data to be transformed.		     C
C                  NN MUST be a power of two. This is checked for by C
C                  the subroutine POWTWO.                            C
C                                                                    C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       C
C                                                                    C
C    The real function f( PHI ) must be periodic, satisfying         C
C    f( PHI ) = f( PHI + 2*PI )                                      C
C                                                                    C
C    f( PHI ) has the formulation                                    C
C                                                                    C
C    f( PHI ) = sum_{m=0}^M [ c_m cos(m PHI) + s_m sin( m PHI) ]     C
C                                                                    C
C    FFTRLV will transform between the coefficients {c_m, s_m}       C
C    and discretised values of f at NN equally spaced points PHI_i   C
C    with                                                            C
C                                                                    C
C    PHI_i = (i-1)*DELTA_phi          and                            C
C                                                                    C
C    DELTA_phi = 2*PI/dble( NN )                                     C
C                                                                    C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       C
C                                                                    C
C     ISIGN	:                                                    C
C                                                                    C
C    If ISIGN = 1, on entry DATA( 2*I - 1 ) must contain the value   C
C    of f at the point PHI = PHI_i as defined above.                 C
C                                                                    C
C    On exit, DATA( 2*m + 1 ) will contain the coeff. c_m and        C
C             DATA( 2*m + 2 ) will contain the coeff. s_m and        C
C                                                                    C
C    for m = 0, M.                                                   C
C      For accuracy, M **MUST** be less than NN/2                    C
C    (This is the Nyquist criterion - see Numerical Recipes ).       C
C                                                                    C
C    If ISIGN = -1, DATA must be zero except for the values          C
C    c_m contained in DATA( 2*m + 1 )  and                           C
C    s_m contained in DATA( 2*m + 2 ) on input.                      C
C                                                                    C
C    On output, f( PHI_i ) will be contained in DATA( 2*I - 1 ).     C
C                                                                    C
C    Note that this routine SCALES the coefficients after            C
C    transforming, which the Numerical Recipes version doesn't.      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     DATA	: Array of dimension (2*NN). See above ...           C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FFTRLV ( DATA, NN, ISIGN )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NN, ISIGN
      DOUBLE PRECISION DATA( 2 * NN )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TEMPR,TEMPI,THETA,WPR,WPI,WR,WI,WTEMP,FAC
      INTEGER I,J,N,M,MMAX,ISTEP
c     LOGICAL POWT
      DOUBLE PRECISION PI, HTHETA
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C 
C  SJG Mon Feb  5 18:18:33 WET 2001
C  I have commented out the power of two checking
C  for speed: reinstate if necessary:
C 
C First check that the number NN supplied is infact a power of 2.
c     CALL POWTWO ( NN, POWT)
c     IF ( .NOT. POWT) THEN
c        PRINT *,' Subroutine FFTRLV. NN is not a power of 2.'
c        PRINT *,'Program aborted.'
c        STOP
c     ENDIF
C Now check that ISIGN = 1 or -1.
      IF ( ISIGN*ISIGN .NE. 1) THEN
         PRINT *,' Subroutine FFTRLV. ISIGN must be 1 or -1.'
         PRINT *,'Program aborted.'
         STOP
      ENDIF
C____________________________________________________________________C
C Here we will zero the imaginary part of the function
C just incase it contains non-zero data which may harm the
C outcome.
C
      IF ( ISIGN.EQ.1 ) THEN
        DO I = 1, NN
          DATA( 2*I ) = 0.0d0
        ENDDO
      ENDIF
C
C____________________________________________________________________C
C
      N = 2 * NN
      J = 1
      DO I = 1, N, 2
         IF ( J.GT.I ) THEN
            TEMPR = DATA( J )
            TEMPI = DATA( J+1 )
            DATA( J ) = DATA( I )
            DATA( J+1 ) = DATA( I+1 )
            DATA( I ) = TEMPR
            DATA( I+1 ) = TEMPI
         ENDIF
C ... now calculate inverse bit map for next value of I.
         M = N/2
 500     IF ( M.GE.2 .AND. J.GT.M ) THEN
            J = J-M
            M = M/2
            GOTO 500
         ENDIF
         J = J + M
      ENDDO
C ............. now do the real part of transform.
      MMAX = 2
 502  IF ( N.GT.MMAX ) THEN
         ISTEP = 2*MMAX
         HTHETA = PI/DBLE(ISIGN*MMAX)
c        THETA = 2.0d0*PI/DBLE(ISIGN*MMAX)
C set THETA temporarily
         THETA = DSIN(HTHETA)
         WPR = -2.0d0*THETA*THETA
c        WPR = -2.0d0*DSIN(HTHETA)**2
         THETA = 2.0d0*HTHETA
         WPI = DSIN(THETA)
         WR = 1.0d0
         WI = 0.0d0
         DO M = 1, MMAX, 2
            DO I = M, N, ISTEP
               J = I + MMAX
               TEMPR = WR*DATA( J ) - WI*DATA( J + 1 )
               TEMPI = WR*DATA( J + 1 ) + WI*DATA( J )
               DATA( J ) = DATA( I ) - TEMPR
               DATA( J + 1 ) = DATA( I + 1 ) - TEMPI
               DATA( I ) = DATA( I ) + TEMPR
               DATA( I + 1 ) = DATA( I + 1 ) + TEMPI
            ENDDO
            WTEMP = WR
            WR = WR*WPR - WI*WPI + WR
            WI = WI*WPR + WTEMP*WPI + WI
         ENDDO
         MMAX = ISTEP
         GOTO 502
      ENDIF
C
      IF ( ISIGN.EQ.1 ) THEN
        FAC = 2.0d0/DBLE( NN )
        DO I = 1, NN
          DATA( I ) = FAC*DATA( I )
        ENDDO
        DO I = NN+1, 2*NN
          DATA( I ) = 0.0d0
        ENDDO
      ENDIF
C
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine Double precision VECtor Zero ****************************
C Steve Gibbons Wed Feb 14 09:01:11 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Sets to zero a double precision vector, VEC, of length N.          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N		: Length of the vector.                              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VEC	: Vector - dimension ( N )                           C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DVECZ( VEC, N )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          N
      DOUBLE PRECISION VEC( N )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          I
      DOUBLE PRECISION DZERO
      PARAMETER      ( DZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO I = 1, N
        VEC( I ) = DZERO
      ENDDO
      RETURN
      END
C*********************************************************************


C*********************************************************************
C function DPMM1 *****************************************************
C Steve Gibbons 18.4.97                                              C
C____________________________________________________________________C
C Calculates d(P_(m+1)^m/d(theta) by equation 185 in my notes        C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     S         : Sin(theta)                                         C
C     PMM0      : P_m^m(X) as evaluated by Function PMM.             C
C     DPMM0     : d(P_m^m(X))/d(theta) as evaluated by Function DPMM C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPMM1 ( M, X, S, PMM0, DPMM0)
      IMPLICIT NONE
      DOUBLE PRECISION DPMM1
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION X, S, PMM0, DPMM0
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function DPMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = DBLE( M )
      DPMM1 = DSQRT( 2.0d0*RM+1.0d0 )*( X*DPMM0 - S*PMM0 )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function PLM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Calculates the Schmidt Normalised Legendre Function P_l^m (x)      C
C given P_(l-1)^m and P_(l-2)^m according to equation 183 in my notesC
C i.e.                                                               C
C   P_l^m( X ) = { A * P_(l-1)^m - B * P_(l-2)^m }/C                 C
C                                                                    C
C where A = (2*l - 1)*X ,                                            C
C                                                                    C
C B = SQRT( (L+M-1)*(L-M-1) ) and C = SQRT( (L+M)*(L-M) )            C
C                                                                    C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C     L         : Well it's L isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     PLMIN1    : P_(l-1)^m ( X )                                    C
C     PLMIN2    : P_(l-2)^m ( X )                                    C
C____________________________________________________________________C
C Local Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     RL        : DBLE  ( L )                                        C
C     RM        : DBLE  ( M )                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PLM ( L, M, X, PLMIN1, PLMIN2 )
      IMPLICIT NONE
      DOUBLE PRECISION PLM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M,L
      DOUBLE PRECISION X, PLMIN1, PLMIN2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM,RL
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function PLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function PLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' PLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = DBLE( M )
      RL = DBLE( L )
      PLM = ( 2.0d0*RL - 1.0d0 )*X*PLMIN1
      PLM = PLM - PLMIN2*DSQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      PLM = PLM/DSQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
