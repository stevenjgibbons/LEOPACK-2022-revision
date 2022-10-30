C*********************************************************************
C                                                                    C
C                  Conducting Inner Core                             C
C                  -          -     -                                C
C Steve Gibbons -  Solution Vector: Perutrbation and New Spatial     C
C                  -        -       -                -   -           C
C Mesh Adaption Prog.                                                C
C -    -        -                                                    C
C Mon Jan  8 12:26:30 WET 2001                                       C
C                                                                    C
C Reads in a solution vector and transforms it onto a spectral/      C
C finite difference mesh which may have more (or fewer) spherical    C
C harmonics, more (or fewer) or differently spaced radial grid nodes C
C                                                                    C
C In SVNSMAP, the new radial functions are zero. This program is     C
C an adaption such that the new radial functions are filled with     C
C very small perturbations (size controlled by parameter DMAG).      C
C                                                                    C
C*********************************************************************
      PROGRAM cicsvpnsmap
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRVMAX, NRMMAX, NHVMAX, NHMMAX, IVVMAX, IVMMAX, LHMAX,
     1        NDCS, NNDM
      PARAMETER ( NRVMAX = 300, NRMMAX = 300, NHVMAX = 6000,
     1            NHMMAX = 6000, LHMAX = 200,
     2            NDCS = LHMAX+4, IVVMAX = NRVMAX*NHVMAX,
     3            IVMMAX = NRMMAX*NHMMAX, NNDM = 6 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          MHIBCI( NDCS ), MHOBCI( NDCS ),
     1                 MHIBCO( NDCS ), MHOBCO( NDCS ),
     2                 LARR( NDCS ), IWORK( NNDM ),
     3                 ISYMA( 5 ), LHARRV( 5 ), LHARRM( 5 ),
     4                 MMODES( LHMAX + 1 )
      DOUBLE PRECISION WORK1( NNDM ), WORK2( NNDM ),
     1                 WORKM( NNDM, NNDM )
      INTEGER          MHTVIN( NHVMAX ), MHLVIN( NHVMAX ), 
     1                 MHMVIN( NHVMAX ), MHPVIN( NHVMAX )
      INTEGER          MHTMIN( NHMMAX ), MHLMIN( NHMMAX ), 
     1                 MHMMIN( NHMMAX ), MHPMIN( NHMMAX )
      INTEGER          MHTVUT( NHVMAX ), MHLVUT( NHVMAX ), 
     1                 MHMVUT( NHVMAX )
      INTEGER          MHTMUT( NHMMAX ), MHLMUT( NHMMAX ), 
     1                 MHMMUT( NHMMAX )
      DOUBLE PRECISION VECVIN( IVVMAX ), VECVUT( IVVMAX ),
     1                 XARRVN( NRVMAX ), XARRVT( NRVMAX )
      DOUBLE PRECISION VECMIN( IVMMAX ), VECMUT( IVMMAX ),
     1                 XARRMN( NRMMAX ), XARRMT( NRMMAX )
      CHARACTER *(80)  FNAME, LINE, ROOT
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER I, ILEN, NCUDS, NRV, NRM, NR1, ISPV,
     1        IFORMF, LU, NRVNW, ISV, M, NNDS, NRICNW, ISPIC,
     2        ISM, MINC, MMAX, NHVUT, NHMUT, NMODES, IS, IFORM,
     3        IHI, IHO, IR, IND, INDFUN, LH1, LH2, LH3, LH4, LH5,
     4        NHV, NHM, NR, IRAD, NRMNW
      INTEGER INARVI( 3 ), INARVO( 3 ),
     1        INARMI( 3 ), INARMO( 3 )
      DOUBLE PRECISION RIV, ROV, RIM, ROM, RAD, RFTOT, DZERO, DMAG,
     1                 DLOW
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IFORM = 1, DZERO = 0.0d0, DLOW = 1.0d-6 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      LU      = 91
 80   FORMAT(A)
C----------------------------------------------------------------
C                Enter filename of velocity harmonics file.     |
C----------------------------------------------------------------
      PRINT *,' Enter name of velocity harmonics file.'
 21   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 21
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 41
        ENDIF
      ENDDO
 41   CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in velocity harmonics file
C
      NCUDS = 0
C     (ncuds - number of diff. schemes already in use).
      CALL HMFRD( NHV, NHVMAX, MHTVIN, MHLVIN, MHMVIN, MHPVIN,
     1            NCUDS, NDCS, MHIBCI, MHOBCI, LARR, LU, FNAME )
      INARVI( 3 ) = NHV
C----------------------------------------------------------------
C                Enter filename of velocity solution vector.    |
C----------------------------------------------------------------
      PRINT *,' Enter name of velocity solution vector file.'
 22   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 22
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 42
        ENDIF
      ENDDO
 42   CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in solution vector file
C
      CALL SVFRD( INARVI, LU, NRVMAX, VECVIN, FNAME )
C
      NR1 = INARVI( 2 )
C----------------------------------------------------------------
C                Enter filename of radial spacing vector.       |
C----------------------------------------------------------------
      PRINT *,' Enter name of radial spacing file.'
 23   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 23
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 43
        ENDIF
      ENDDO
 43   CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in radial spacing file
C
      CALL XARRRD( NRV, NRVMAX, XARRVN, LU, FNAME )
C
      IF ( NR1.NE.NRV ) THEN
        PRINT *,' Solution vector and radial node '
        PRINT *,' file claim differing numbers of '
        PRINT *,' grid nodes. Program aborted.'
        STOP
      ENDIF
C
      RIV = XARRVN(  1  )
      ROV = XARRVN( NRV )
C----------------------------------------------------------------
C                Enter filename of magnetic harmonics file.     |
C----------------------------------------------------------------
      PRINT *,' Enter name of magnetic harmonics file.'
 121  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 121
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 141
        ENDIF
      ENDDO
 141  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in magnetic harmonics file
C
      CALL HMFRD( NHM, NHMMAX, MHTMIN, MHLMIN, MHMMIN, MHPMIN,
     1            NCUDS, NDCS, MHIBCI, MHOBCI, LARR, LU, FNAME )
      INARMI( 3 ) = NHM
C----------------------------------------------------------------
C                Enter filename of magnetic solution vector.    |
C----------------------------------------------------------------
      PRINT *,' Enter name of magnetic solution vector file.'
 122  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 122
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 142
        ENDIF
      ENDDO
 142  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in magnetic solution vector file
C
      CALL SVFRD( INARMI, LU, NRMMAX, VECMIN, FNAME )
C
      NR1 = INARMI( 2 )
C----------------------------------------------------------------
C        Enter filename of magnetic radial spacing vector.      |
C----------------------------------------------------------------
      PRINT *,' Enter name of magnetic radial spacing file.'
 123  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 123
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 143
        ENDIF
      ENDDO
 143  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in magnetic radial spacing file
C
      CALL XARRRD( NRM, NRMMAX, XARRMN, LU, FNAME )
C
      IF ( NR1.NE.NRM ) THEN
        PRINT *,' Solution vector and radial node '
        PRINT *,' file claim differing numbers of '
        PRINT *,' grid nodes. Program aborted.'
        STOP
      ENDIF
C
      RIM = XARRMN(  1  )
      ROM = XARRMN( NRM )
C
      IF ( DABS( ROV - ROM ).GT.DLOW ) THEN
        PRINT *,' ROV = ', ROV,' ROM = ', ROM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C----------------------------------------------------------------
C                Enter stem of filename for output vectors      |
C----------------------------------------------------------------
      PRINT *,' Enter name of filename stem.'
 31   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 31
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 45
        ENDIF
      ENDDO
 45   CONTINUE
      ROOT  = LINE(1:ILEN)
C----------------------------------------------------------------
C                Enter stem of filename for output vectors      |
C----------------------------------------------------------------
      PRINT *,' Enter NRVNW, ISPV, NRICNW, ISPIC, IFORMF, NNDS, DMAG '
      PRINT *,' nrvnw is the new number of grid nodes for outer core.'
      PRINT *,' ispv = 1 --> equally spaced nodes.'
      PRINT *,' ispv = 2 --> Chebyshev nodes.'
      PRINT *,' nricnw the new number of grid nodes for inner core.'
      PRINT *,' ispic = 1 --> equally spaced nodes.'
      PRINT *,' ispic = 2 --> Chebyshev nodes over (-1,1).'
      PRINT *,' ispic = 3 --> Chebyshev nodes over  (0,1).'
      PRINT *,' iformf = 3 --> ind = ( IR - 1 )*NH + IH '
      PRINT *,' iformf = 4 --> ind = ( IH - 1 )*NR + IR '
      PRINT *,' nnds = nodes required for interpolation.'
      PRINT *,' 3 is generally sufficient.'
      PRINT *,' no more than ', NNDM
      PRINT *,' DMAG = amplitude of random functions'
 32   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 32
      READ ( LINE, * ) NRVNW, ISPV, NRICNW, ISPIC,
     1                  IFORMF, NNDS, DMAG
C
      IF ( NRVNW.GT.NRVMAX ) THEN
        PRINT *,' NRVNW = ', NRVNW,' NRVMAX = ', NRVMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ISPV.NE.1 .AND. ISPV.NE.2 ) THEN
        PRINT *,' ISPV = ', ISPV
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NRMNW = NRICNW+NRVNW
      IF ( NRMNW.GT.NRMMAX ) THEN
        PRINT *,' NRVNW  = ', NRVNW
        PRINT *,' NRICNW = ', NRICNW
        PRINT *,' NRMMAX = ', NRMMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ISPIC.NE.1 .AND. ISPIC.NE.2 .AND. ISPIC.NE.3 ) THEN
        PRINT *,' ISPIC = ', ISPIC
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
      IF ( NNDS.LT.2 .OR. NNDS.GT.NNDM ) THEN
        PRINT *,' NNDS = ', NNDS,' NNDM = ', NNDM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C----------------------------------------------------------------
C                Enter dimensions of the new vector.            |
C----------------------------------------------------------------
C
      PRINT *,' Enter LH1, LH2, LH3, LH4, LH5.'
      PRINT *,' LH1  = max. degree, l, for poloidal vel.'
      PRINT *,' LH2  = max. degree, l, for toroidal vel.'
      PRINT *,' LH3  = max. degree, l, for temperature.'
      PRINT *,' LH4  = max. degree, l, for poloidal mag. field.'
      PRINT *,' LH5  = max. degree, l, for toroidal mag. field.'
C
 33   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 33
      READ ( LINE, * ) LH1, LH2, LH3, LH4, LH5
C
C----------------------------------------------------------------
C                Enter dimensions of the new vector.            |
C----------------------------------------------------------------
C
      PRINT *,' Enter ISV, ISM, MINC, MMAX '
      PRINT *,' ISV  = 1 --> velocity is eq. symmetric.'
      PRINT *,' ISV  = 2 --> velocity is eq. anti-symmetric.'
      PRINT *,' ISV  = 3 --> velocity has both symmetries.'
      PRINT *,' ISM  = 1 --> mag. field is eq. symmetric.'
      PRINT *,' ISM  = 2 --> mag. field is eq. anti-symmetric.'
      PRINT *,' ISM  = 3 --> mag. field has both symmetries.'
      PRINT *,' MINC is increment of wavenumber.'
      PRINT *,' MMAX is maximum desired wavenumber.'
C
 34   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 34
      READ ( LINE, * ) ISV, ISM, MINC, MMAX
C
      LHARRV( 1 ) = LH1
      LHARRV( 2 ) = LH2
      LHARRV( 3 ) = LH3
      LHARRV( 4 ) = 0
      LHARRV( 5 ) = 0
C
      LHARRM( 1 ) = 0
      LHARRM( 2 ) = 0
      LHARRM( 3 ) = -1
      LHARRM( 4 ) = LH4
      LHARRM( 5 ) = LH5
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
      ISYMA( 1 ) = ISV
      ISYMA( 2 ) = ISV
      ISYMA( 3 ) = ISV
      ISYMA( 4 ) = ISM
      ISYMA( 5 ) = ISM
C
      CALL HMINDA( LHARRV, ISYMA, NMODES, MMODES, NHVUT, NHVMAX,
     1                   MHTVUT, MHLVUT, MHMVUT, LHMAX )
      PRINT *,' Total velocity harmonics selected = ', NHVUT
      CALL HMINDA( LHARRM, ISYMA, NMODES, MMODES, NHMUT, NHMMAX,
     1                   MHTMUT, MHLMUT, MHMMUT, LHMAX )
      PRINT *,' Total magnetic harmonics selected = ', NHMUT
      IF ( NHMUT.EQ.0 .OR. NHVUT.EQ.0 ) THEN
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now fill the boundary condition arrays:
C
      DO I = 1, NHV
C
        IF ( MHTVIN( I ).EQ.1 ) THEN
          IS   = MHPVIN( I )
          MHIBCO( 1 ) = MHIBCI( IS )
          MHOBCO( 1 ) = MHOBCI( IS )
        ENDIF
C
        IF ( MHTVIN( I ).EQ.2 ) THEN
          IS   = MHPVIN( I )
          MHIBCO( 2 ) = MHIBCI( IS )
          MHOBCO( 2 ) = MHOBCI( IS )
        ENDIF
C
        IF ( MHTVIN( I ).EQ.3 ) THEN
          IS   = MHPVIN( I )
          MHIBCO( 3 ) = MHIBCI( IS )
          MHOBCO( 3 ) = MHOBCI( IS )
        ENDIF
C
      ENDDO
      DO I = 1, NHM
C
        IF ( MHTMIN( I ).EQ.4 ) THEN
          IS   = MHPMIN( I )
          MHIBCO( 4 ) = MHIBCI( IS )
          MHOBCO( 4 ) = MHOBCI( IS )
        ENDIF
C
        IF ( MHTMIN( I ).EQ.5 ) THEN
          IS   = MHPMIN( I )
          MHIBCO( 5 ) = MHIBCI( IS )
          MHOBCO( 5 ) = MHOBCI( IS )
        ENDIF
C
      ENDDO
C
C Fill radial grid node array for velocity
C
      IF ( ISPV.EQ.1 ) THEN
        CALL ESNAAS( NRVNW, XARRVT, RIV, ROV )
      ELSE
        CALL ZCPAAS( NRVNW, XARRVT, RIV, ROV )
      ENDIF
C
      INARVO( 1 ) = IFORMF
      INARVO( 2 ) = NRVNW
      INARVO( 3 ) = NHVUT
C
      INARMO( 1 ) = IFORMF
      INARMO( 2 ) = NRMNW
      INARMO( 3 ) = NHMUT
C
C Fill radial grid node array for inner core
C
      NR = NRICNW + 1
      IF ( ISPIC.EQ.1 ) CALL ESNAAS( NR, XARRMT, RIM, RIV )
      IF ( ISPIC.EQ.2 ) CALL ZCPAAS( NR, XARRMT, RIM, RIV )
      IF ( ISPIC.EQ.3 ) CALL ZCPAA2( NR, XARRMT, RIM, RIV )
C
C Copy the outer-core nodes into XARRMT
C
      DO IR = 1, NRVNW
        IRAD = IR + NRICNW
        XARRMT( IRAD ) = XARRVT( IR )
      ENDDO
C
C Zero VECVUT array
C
      DO I = 1, NRVNW*NHVUT
        VECVUT( I ) = 0.0d0
      ENDDO
C
C Zero VECMUT array
C
      DO I = 1, NRMNW*NHMUT
        VECMUT( I ) = 0.0d0
      ENDDO
C
C Now interpolate velocity
C
      DO IHI = 1, NHV
        DO IHO = 1, NHVUT
          IF (  MHTVIN( IHI ).EQ.MHTVUT( IHO ) .AND.
     1          MHLVIN( IHI ).EQ.MHLVUT( IHO ) .AND.
     2          MHMVIN( IHI ).EQ.MHMVUT( IHO )   ) THEN
            DO IR = 1, NRVNW
              RAD = XARRVT( IR )
              CALL SVRINT( RAD, VECVIN, XARRVN, INARVI, IHI, NNDS,
     1                     WORK1, IWORK, WORK2, WORKM )
              IND = INDFUN( IR, IHO, INARVO )
              VECVUT( IND ) = WORK1( 1 )
            ENDDO
          ENDIF
        ENDDO
      ENDDO
C
C Now interpolate magnetic field
C
      DO IHI = 1, NHM
        DO IHO = 1, NHMUT
          IF (  MHTMIN( IHI ).EQ.MHTMUT( IHO ) .AND.
     1          MHLMIN( IHI ).EQ.MHLMUT( IHO ) .AND.
     2          MHMMIN( IHI ).EQ.MHMMUT( IHO )   ) THEN
            DO IR = 1, NRMNW
              RAD = XARRMT( IR )
              CALL SVRINT( RAD, VECMIN, XARRMN, INARMI, IHI, NNDS,
     1                     WORK1, IWORK, WORK2, WORKM )
              IND = INDFUN( IR, IHO, INARMO )
              VECMUT( IND ) = WORK1( 1 )
            ENDDO
          ENDIF
        ENDDO
      ENDDO
C
C Ok - now we need to put in the 'new' radial for velocity
C functions. First loop around all the NHVUT
C harmonics and decide whether or not they are
C zero.
C
      DO IHO = 1, NHVUT
        RFTOT = 0.0d0
        DO IR = 1, NRVNW
          IND = INDFUN( IR, IHO, INARVO )
            RFTOT = RFTOT +  VECVUT( IND )*VECVUT( IND )
        ENDDO
C
C If rftot is zero, then we need to put a perturbation
C function into the radial function.
C
        IF ( RFTOT.EQ.DZERO ) THEN
C         .
          DO IR = 1, NRVNW
            RAD = XARRVT( IR )
            IND = INDFUN( IR, IHO, INARVO )
            VECVUT( IND ) = DMAG*DSIN( RAD )
          ENDDO
C         .
        ENDIF
      ENDDO
C
C Ok - now we need to put in the 'new' radial for magnetic
C functions. First loop around all the NHVUT
C harmonics and decide whether or not they are
C zero.
C
      DO IHO = 1, NHMUT
        RFTOT = 0.0d0
        DO IR = 1, NRMNW
          IND = INDFUN( IR, IHO, INARMO )
            RFTOT = RFTOT +  VECMUT( IND )*VECMUT( IND )
        ENDDO
C
C If rftot is zero, then we need to put a perturbation
C function into the radial function.
C
        IF ( RFTOT.EQ.DZERO ) THEN
C         .
          DO IR = 1, NRMNW
            RAD = XARRMT( IR )
            IND = INDFUN( IR, IHO, INARMO )
            VECMUT( IND ) = DMAG*DSIN( RAD )
          ENDDO
C         .
        ENDIF
      ENDDO
C
C Finally: write out all of the files
C
      FNAME(1:ILEN) = ROOT(1:ILEN)
      FNAME(ILEN+1:ILEN+6) = '.intsv'
      FNAME = FNAME(1:ILEN+6)
C
      CALL HMFWT( NHVUT, MHTVUT, MHLVUT, MHMVUT, MHTVUT, NDCS,
     1            MHIBCO, MHOBCO, LU, FNAME )
C
      FNAME(ILEN+1:ILEN+6) = '.xarrv'
      FNAME = FNAME(1:ILEN+6)
      CALL XARRWT( NRVNW, XARRVT, LU, FNAME, IFORM )
C
      FNAME(ILEN+1:ILEN+6) = '.vecsv'
      FNAME = FNAME(1:ILEN+6)
      CALL SVFWT( INARVO, LU, IFORM, VECVUT, FNAME )
C
      FNAME(ILEN+1:ILEN+6) = '.intsm'
      FNAME = FNAME(1:ILEN+6)
C
      CALL HMFWT( NHMUT, MHTMUT, MHLMUT, MHMMUT, MHTMUT, NDCS,
     1            MHIBCO, MHOBCO, LU, FNAME )
C
      FNAME(ILEN+1:ILEN+6) = '.xarrm'
      FNAME = FNAME(1:ILEN+6)
      CALL XARRWT( NRMNW, XARRMT, LU, FNAME, IFORM )
C
      FNAME(ILEN+1:ILEN+6) = '.vecsm'
      FNAME = FNAME(1:ILEN+6)
      CALL SVFWT( INARMO, LU, IFORM, VECMUT, FNAME )
      STOP
      END
C*********************************************************************
