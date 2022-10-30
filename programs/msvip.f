C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C Mon Mar  6 10:14:00 GMT 2000                                       C
C                                                                    C
C Multiple Solution Vector Ingestion Program.                        C
C                                                                    C
C Reads in from standard input a filename stem, radial spacing info, C
C symmetry info etc. and boundary conditions and fills (using the    C
C routine SHVECF) a very arbitrary non-zero function into the        C
C solution vector.                                                   C
C                                                                    C
C It then loops around a list of file names for                      C
C .ints, .xarr and .vecs files (and FAC1, FAC2 coefficients)         C
C adding FAC2 of the new vector to FAC1 * the old vector.            C
C                                                                    C
C When the list is completely read,                                  C
C ROOT.ints, ROOT.vecs and ROOT.xarr are output is normal.           C
C                                                                    C
C*********************************************************************
      PROGRAM msvip
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, NHMAX, LHMAX, NDCS, ISVMAX, NNDM
      PARAMETER ( NRMAX = 300, NHMAX = 6000, LHMAX = 104,
     1            NDCS  = LHMAX + 6, ISVMAX = NRMAX*NHMAX,
     2            NNDM = 6 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ), MHP( NHMAX ),
     1        MHIBC( NDCS ), MHOBC( NDCS ), ISYMA( 5 ), 
     2        LHARR( 5 ), MMODES( LHMAX + 1 ), IWORK( NNDM ),
     3        MHT2( NHMAX ), MHL2( NHMAX ), MHM2( NHMAX )
C
      DOUBLE PRECISION XARR( NRMAX ), VEC( ISVMAX ),
     1                 XARR2( NRMAX ), VEC2( ISVMAX ), WORKA( NNDM ),
     2                 WORKB( NNDM ), WORKC( NNDM, NNDM )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      CHARACTER *(80) LINE, ROOT, FNAME, FNNEW
      INTEGER I, ILEN, NR, INSF, IFORMF, INARR( 3 ), IH, LU,
     1        ISM, ISV, LHM, LHV, MINC, MMAX, M, NH, NMODES,
     2        IVELBC, ITHEBC, IFORM, NR1, NR2, INARR2( 3 ),
     3        NNDS, ILEN2, NH2
      DOUBLE PRECISION RI, RO, LOW, ZERO, FAC, FAC1, FAC2
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
      PRINT *,' Enter NR, IFORMF, NNDS.'
      PRINT *,' nrmax = ', NRMAX
      PRINT *,' iformf = 3 --> ind = ( ir - 1 )*nh + ih '
      PRINT *,' iformf = 4 --> ind = ( ih - 1 )*nr + ir '
      PRINT *,' nnds: nodes for interpolation.'
      PRINT *,' nndm = ', NNDM
      PRINT *,' ----------------------------------------'
 32   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 32
      READ ( LINE, * ) NR, IFORMF, NNDS
      PRINT *,' NR = ', NR,' IFORMF = ', IFORMF,' NNDS = ', NNDS
C
C Check values of NR, IFORMF, NNDS
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
      IF ( NNDS.GT.NNDM ) THEN
        PRINT *,' NNDS = ', NNDS,' NNDM = ', NNDM
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
C Loop around reading in groups of 4 lines.
C The first must contain 2 numbers: FAC1 FAC2
C The second must contain the harmonics (.ints) file
C The third must contain the vector (.vecs) file
C The fourth must contain the spacings (.xarr) file
C
 36   CONTINUE
      READ ( 5, 80, END=301 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 36
      READ ( LINE, * ) FAC1, FAC2
C
C Now read in integer harmonics file
C
 37   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 37
C
      DO I = 1, 80
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN2 = I-1
          GOTO 46
        ENDIF
      ENDDO
 46   CONTINUE
      FNNEW = LINE(1:ILEN2)
C
C Read in harmonics file
C
      CALL BIHFRD( NH2, NHMAX, MHT2, MHL2, MHM2, LU, FNNEW )
      INARR2( 3 ) = NH2
C
C Now read in vectors file
C
 38   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 38
C
      DO I = 1, 80
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN2 = I-1
          GOTO 47
        ENDIF
      ENDDO
 47   CONTINUE
      FNNEW = LINE(1:ILEN2)
C
      CALL SVFRD( INARR2, LU, NRMAX, VEC2, FNNEW )
      NR1 = INARR2( 2 )
C
C Now read in radial spacings filename
C
 39   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 39
C
      DO I = 1, 80
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN2 = I-1
          GOTO 48
        ENDIF
      ENDDO
 48   CONTINUE
      FNNEW = LINE(1:ILEN2)
C
      CALL XARRRD( NR2, NRMAX, XARR2, LU, FNNEW )
C
      IF ( NR1.NE.NR2 ) THEN
        PRINT *,' Solution vector and radial node '
        PRINT *,' file claim differing numbers of '
        PRINT *,' grid nodes. Program aborted.'
        STOP
      ENDIF
C
C Now incorporate new vector into old.
C
      CALL ASSVAR( INARR, MHT, MHL, MHM, INARR2, MHT2, MHL2,
     1             MHM2, NNDS, IWORK, VEC, XARR, VEC2, XARR2,
     2             WORKA, WORKB, WORKC, FAC1, FAC2 )     
C
      GOTO 36
C
 301  CONTINUE
C
C**********************************************************C
C  Now write out the accumulated file.                     C
C**********************************************************C
C
      FNAME(1:ILEN) = ROOT(1:ILEN)
      FNAME(ILEN+1:ILEN+5) = '.ints'
      FNAME = FNAME(1:ILEN+5)
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
