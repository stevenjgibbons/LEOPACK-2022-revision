C*********************************************************************
C                                                                    C
C Steve Gibbons -  Solution Vector: Perutrbation and New Spatial     C
C                  -        -       -                -   -           C
C Mesh Adaption Prog.                                                C
C -    -        -                                                    C
C Mon May 15 08:03:32 BST 2000                                       C
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
      PROGRAM svpnsmap
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, NHMAX, ISVMAX, LHMAX, NDCS, NNDM
      PARAMETER ( NRMAX = 300, NHMAX = 6000, LHMAX = 100,
     1            NDCS = LHMAX+4, ISVMAX = NRMAX*NHMAX, NNDM = 6 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          MHIBCI( NDCS ), MHOBCI( NDCS ),
     1                 MHIBCO( NDCS ), MHOBCO( NDCS ),
     2                 LARR( NDCS ), IWORK( NNDM ),
     3                 ISYMA( 5 ), LHARR( 5 ), MMODES( LHMAX + 1 )
      DOUBLE PRECISION WORK1( NNDM ), WORK2( NNDM ),
     1                 WORKM( NNDM, NNDM )
      INTEGER          MHTIN( NHMAX ), MHLIN( NHMAX ), 
     1                 MHMIN( NHMAX ), MHPIN( NHMAX )
      INTEGER          MHTOUT( NHMAX ), MHLOUT( NHMAX ), 
     1                 MHMOUT( NHMAX )
      DOUBLE PRECISION VECIN( ISVMAX ), VECOUT( ISVMAX ),
     1                 XARRIN( NRMAX ), XARROT( NRMAX )
      CHARACTER *(80)  FNAME, LINE, ROOT
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER I, ILEN, NCUDS, NH, NR, ISP, INARRI( 3 ), NR1,
     1        INARRO( 3 ), IFORMF, LU, NRNEW, LHV, ISV, M, NNDS,
     2        LHM, ISM, MINC, MMAX, NHOUT, NMODES, IS, IFORM,
     3        IHI, IHO, IR, IND, INDFUN
      DOUBLE PRECISION RI, RO, RAD, RFTOT, DMAG, DZERO
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IFORM = 1, DZERO = 0.0d0 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      LU      = 91
 80   FORMAT(A)
C----------------------------------------------------------------
C                Enter filename of harmonics file.              |
C----------------------------------------------------------------
      PRINT *,' Enter name of harmonics file.'
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
C Read in harmonics file
C
      NCUDS = 0
C     (ncuds - number of diff. schemes already in use).
      CALL HMFRD( NH, NHMAX, MHTIN, MHLIN, MHMIN, MHPIN, NCUDS, NDCS,
     1            MHIBCI, MHOBCI, LARR, LU, FNAME )
      INARRI( 3 ) = NH
C----------------------------------------------------------------
C                Enter filename of solution vector.             |
C----------------------------------------------------------------
      PRINT *,' Enter name of solution vector file.'
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
      CALL SVFRD( INARRI, LU, NRMAX, VECIN, FNAME )
C
      NR1 = INARRI( 2 )
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
      CALL XARRRD( NR, NRMAX, XARRIN, LU, FNAME )
C
      IF ( NR1.NE.NR ) THEN
        PRINT *,' Solution vector and radial node '
        PRINT *,' file claim differing numbers of '
        PRINT *,' grid nodes. Program aborted.'
        STOP
      ENDIF
C
      RI = XARRIN(  1 )
      RO = XARRIN( NR )
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
      PRINT *,' Enter NRNEW, ISP, IFORMF, NNDS, DMAG '
      PRINT *,' nrnew is the new number of radial grid nodes.'
      PRINT *,' isp = 1 --> equally spaced nodes.'
      PRINT *,' isp = 2 --> Chebyshev nodes.'
      PRINT *,' iformf = 3 --> ind = ( IR - 1 )*NH + IH '
      PRINT *,' iformf = 4 --> ind = ( IH - 1 )*NR + IR '
      PRINT *,' nnds = nodes required for interpolation.'
      PRINT *,' 3 is generally sufficient.'
      PRINT *,' no more than ', NNDM
 32   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 32
      READ ( LINE, * ) NRNEW, ISP, IFORMF, NNDS, DMAG
C
      IF ( NRNEW.GT.NRMAX ) THEN
        PRINT *,' NRNEW = ', NRNEW,' NRMAX = ', NRMAX
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
 33   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 33
      READ ( LINE, * ) LHV, ISV, LHM, ISM, MINC, MMAX
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
      ISYMA( 1 ) = ISV
      ISYMA( 2 ) = ISV
      ISYMA( 3 ) = ISV
      ISYMA( 4 ) = ISM
      ISYMA( 5 ) = ISM
C
      CALL HMINDA( LHARR, ISYMA, NMODES, MMODES, NHOUT, NHMAX,
     1                   MHTOUT, MHLOUT, MHMOUT, LHMAX )
      PRINT *,' Total harmonics selected = ', NHOUT
      IF ( NHOUT.EQ.0 ) THEN
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now fill the boundary condition arrays:
C (First set some given defaults, incase
C there were no harmonics of that type in
C the original solution vector).
C
C First, poloidal velocity: no slip boundary
C
      MHIBCO( 1 ) = 4
      MHOBCO( 1 ) = 4
C
C Now, toroidal velocity: no slip boundary
C
      MHIBCO( 2 ) = 2
      MHOBCO( 2 ) = 2
C
C Now, temperature: constant at boundary
C
      MHIBCO( 3 ) = 2
      MHOBCO( 3 ) = 2
C
C Now, poloidal magnetic field: potential field
C
      MHIBCO( 4 ) = 7
      MHOBCO( 4 ) = 7
C
C Now, toroidal magnetic field: potential field
C
      MHIBCO( 5 ) = 2
      MHOBCO( 5 ) = 2
C
      DO I = 1, NH
C
        IF ( MHTIN( I ).EQ.1 ) THEN
          IS   = MHPIN( I )
          MHIBCO( 1 ) = MHIBCI( IS )
          MHOBCO( 1 ) = MHOBCI( IS )
        ENDIF
C
        IF ( MHTIN( I ).EQ.2 ) THEN
          IS   = MHPIN( I )
          MHIBCO( 2 ) = MHIBCI( IS )
          MHOBCO( 2 ) = MHOBCI( IS )
        ENDIF
C
        IF ( MHTIN( I ).EQ.3 ) THEN
          IS   = MHPIN( I )
          MHIBCO( 3 ) = MHIBCI( IS )
          MHOBCO( 3 ) = MHOBCI( IS )
        ENDIF
C
        IF ( MHTIN( I ).EQ.4 ) THEN
          IS   = MHPIN( I )
          MHIBCO( 4 ) = MHIBCI( IS )
          MHOBCO( 4 ) = MHOBCI( IS )
        ENDIF
C
        IF ( MHTIN( I ).EQ.5 ) THEN
          IS   = MHPIN( I )
          MHIBCO( 5 ) = MHIBCI( IS )
          MHOBCO( 5 ) = MHOBCI( IS )
        ENDIF
C
      ENDDO
C
C Fill radial grid node array
C
      IF ( ISP.EQ.1 ) THEN
        CALL ESNAAS( NRNEW, XARROT, RI, RO )
      ELSE
        CALL ZCPAAS( NRNEW, XARROT, RI, RO )
      ENDIF
C
      INARRO( 1 ) = IFORMF
      INARRO( 2 ) = NRNEW
      INARRO( 3 ) = NHOUT
C
C Zero VECOUT array
C
      DO I = 1, NRNEW*NHOUT
        VECOUT( I ) = 0.0d0
      ENDDO
C
C Now interpolate
C
      DO IHI = 1, NH
        DO IHO = 1, NHOUT
          IF (  MHTIN( IHI ).EQ.MHTOUT( IHO ) .AND.
     1          MHLIN( IHI ).EQ.MHLOUT( IHO ) .AND.
     2          MHMIN( IHI ).EQ.MHMOUT( IHO )   ) THEN
            DO IR = 1, NRNEW
              RAD = XARROT( IR )
              CALL SVRINT( RAD, VECIN, XARRIN, INARRI, IHI, NNDS,
     1                     WORK1, IWORK, WORK2, WORKM )
              IND = INDFUN( IR, IHO, INARRO )
              VECOUT( IND ) = WORK1( 1 )
            ENDDO
          ENDIF
        ENDDO
      ENDDO
C
C Ok - now we need to put in the 'new' radial
C functions. First loop around all the NHOUT
C harmonics and decide whether or not they are
C zero.
C
      DO IHO = 1, NHOUT
        RFTOT = 0.0d0
        DO IR = 1, NRNEW
          IND = INDFUN( IR, IHO, INARRO )
            RFTOT = RFTOT +  VECOUT( IND )*VECOUT( IND )
        ENDDO
C
C If rftot is zero, then we need to put a perturbation
C function into the radial function.
C
        IF ( RFTOT.EQ.DZERO ) THEN
C         .
          DO IR = 1, NRNEW
            RAD = XARROT( IR )
            IND = INDFUN( IR, IHO, INARRO )
            VECOUT( IND ) = VECOUT( IND ) + DMAG*DSIN( RAD )
          ENDDO
C         .
        ENDIF
      ENDDO
C
C Finally: write out all of the files
C
      FNAME(1:ILEN) = ROOT(1:ILEN)
      FNAME(ILEN+1:ILEN+5) = '.ints'
      FNAME = FNAME(1:ILEN+5)
C
      CALL HMFWT( NHOUT, MHTOUT, MHLOUT, MHMOUT, MHTOUT, NDCS,
     1            MHIBCO, MHOBCO, LU, FNAME )
C
      FNAME(ILEN+1:ILEN+5) = '.xarr'
      FNAME = FNAME(1:ILEN+5)
      CALL XARRWT( NRNEW, XARROT, LU, FNAME, IFORM )
C
      FNAME(ILEN+1:ILEN+5) = '.vecs'
      FNAME = FNAME(1:ILEN+5)
      CALL SVFWT( INARRO, LU, IFORM, VECOUT, FNAME )
C
      STOP
      END
C*********************************************************************
