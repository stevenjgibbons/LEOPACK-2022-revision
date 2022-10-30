C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C Sat Apr  1 13:17:45 BST 2000                                       C
C                                                                    C
C Insulating Inner Core 2 Conducting Inner Core Solution Convert.    C
C -          -     -    - -          -     -    -        -           C
C                                                                    C
C Reads in a .vec, .ints and .xarr array from a solution from the    C
C insulating inner core codes and splits it into two separate        C
C solution vectors (one for the velocity and temperature, one for    C
C the magnetic field) to be used by the conducting inner core code.  C
C                                                                    C
C The user is asked for a number of nodes and format for the inner   C
C core solution.                                                     C
C                                                                    C
C                                                                    C
C*********************************************************************
      PROGRAM iic2cicsc
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, NHMAX, ISVMAX, LHMAX, NDCS
      PARAMETER ( NRMAX = 300, NHMAX = 2500, LHMAX = 100,
     1            NDCS = LHMAX+4, ISVMAX = NRMAX*NHMAX )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          MHIBCI( NDCS ), MHOBCI( NDCS )
      INTEGER          INARRI( 3 ), INARRV( 3 ), INARRM( 3 ),
     1                 LARR( NDCS )
      INTEGER       MHTIN( NHMAX ), MHLIN( NHMAX ), MHMIN( NHMAX ),
     1              MHPIN( NHMAX ), MHTVT( NHMAX ), MHLVT( NHMAX ),
     2              MHMVT( NHMAX ), MHPVT( NHMAX ), MHTMF( NHMAX ),
     3              MHLMF( NHMAX ), MHMMF( NHMAX ), MHPMF( NHMAX )
      DOUBLE PRECISION XARRIN( NRMAX ), XARRMF( NRMAX ),
     1                 VECIN( ISVMAX ), VECV( ISVMAX ),
     2                 VECM( ISVMAX )
      CHARACTER *(80)  FNAME, LINE, ROOT
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER LU, ILEN, I, NCUDS, NH, NR, NR1, INDFUN, NHV, NHM,
     1        NRIC, ISPIC, IR, IH, IOP, NRMF, INDIN, INDOUT, IHV,
     2        IHM, IRMF, IFORM
      DOUBLE PRECISION RI, DLOW, ZERO, RAD, VAL
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( DLOW = 1.0d-8, ZERO = 0.0d0, IFORM = 1 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Zero vector arrays
C
      CALL VECOP( VECIN, ZERO, ISVMAX, IOP )
      CALL VECOP( VECV , ZERO, ISVMAX, IOP )
      CALL VECOP( VECM , ZERO, ISVMAX, IOP )
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
C
      IF ( RI.LT.DLOW ) THEN
        PRINT *,' RI = ', RI
        PRINT *,' How can you include an inner core!!'
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
C
C OK - input number of nodes for inner core,
C and get format: ISPIC = 1 --> uniform grid nodes
C                 ISPIC = 2 --> Chebyshev nodes (-1,1)
C                 ISPIC = 3 --> Chebyshev nodes ( 0,1)
C
      PRINT *,' Enter NRIC, ISPIC '
      PRINT *,' NRIC = number of radial nodes for inner core.'
      PRINT *,' ISPIC = 1 --> uniform grid nodes.'
      PRINT *,' ISPIC = 2 --> Chebyshev nodes.'
      READ ( 5, * ) NRIC, ISPIC
C
      NRMF = NRIC+NR
      IF ( NRMF.GT.NRMAX ) THEN
        PRINT *,' Total number of nodes required = ', NRMF
        PRINT *,' Maximum nodes allowed = ', NRMAX
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
C Fill the XARRMF array
C
      IF ( NRIC.EQ.0 ) GOTO 70
      IF ( ISPIC.EQ.1 ) THEN
        IH = NRIC+1
        CALL ESNAAS( IH, XARRMF, ZERO, RI)
      ENDIF
      IF ( ISPIC.EQ.2 ) THEN
        IH = NRIC+1
        CALL ZCPAAS( IH, XARRMF, ZERO, RI)
      ENDIF
      IF ( ISPIC.EQ.3 ) THEN
        IH = NRIC+1
        CALL ZCPAA2( IH, XARRMF, ZERO, RI)
      ENDIF
C
 70   CONTINUE
      DO IR = 1, NR
        XARRMF( IR + NRIC ) = XARRIN( IR )
      ENDDO
C
C Loop around harmonics, counting up those of different type:
C Filling new arrays
C
      NHV = 0
      NHM = 0
      DO IH = 1, NH
C ------*
        IF (      MHTIN( IH ).EQ.1      .OR.
     1            MHTIN( IH ).EQ.2      .OR.
     2            MHTIN( IH ).EQ.3    )      THEN
C       ------------------------------------------
C       This is a velocity or temperature harmonic
C       ------------------------------------------
          NHV          = NHV + 1
          MHTVT( NHV ) = MHTIN( IH )
          MHLVT( NHV ) = MHLIN( IH )
          MHMVT( NHV ) = MHMIN( IH )
          MHPVT( NHV ) = MHPIN( IH )
        ENDIF
C ------*
        IF (      MHTIN( IH ).EQ.4      .OR.
     1            MHTIN( IH ).EQ.5    )      THEN
C       ------------------------------------------
C       This is a magnetic field harmonic
C       ------------------------------------------
          NHM          = NHM + 1
          MHTMF( NHM ) = MHTIN( IH )
          MHLMF( NHM ) = MHLIN( IH )
          MHMMF( NHM ) = MHMIN( IH )
          MHPMF( NHM ) = MHPIN( IH )
        ENDIF
C ------*
      ENDDO
C
C Just check that the numbers add up.
C
      IF ( NH.NE.(NHV+NHM) ) THEN
        PRINT *,' Number of harmonics does not add up.'
        PRINT *,' NH = ', NH,' NHV = ', NHV,' NHM = ', NHM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NHM.EQ.0 ) THEN
        PRINT *,' NHM = ', NHM
        PRINT *,' There is no point running this program'
        PRINT *,' with no magnetic field!'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now we can fill INARRV and INARRM
C
      INARRV( 1 ) = INARRI( 1 )
      INARRM( 1 ) = INARRI( 1 )
C
      INARRV( 2 ) = INARRI( 2 )
      INARRM( 2 ) = NRMF
C
      INARRV( 3 ) = NHV
      INARRM( 3 ) = NHM
C
C Now fill our output vectors ...
C
      IHV = 0
      IHM = 0
      DO IH = 1, NH
C ------*
        IF (      MHTIN( IH ).EQ.1      .OR.
     1            MHTIN( IH ).EQ.2      .OR.
     2            MHTIN( IH ).EQ.3    )      THEN
C       ------------------------------------------
C       This is a velocity or temperature harmonic
C       ------------------------------------------
          IHV          = IHV + 1
          DO IR = 1, NR
            INDIN  = INDFUN( IR,  IH, INARRI )
            INDOUT = INDFUN( IR, IHV, INARRV )
            VECV( INDOUT ) = VECIN( INDIN )
          ENDDO
        ENDIF
C ------*
        IF (      MHTIN( IH ).EQ.4      .OR.
     1            MHTIN( IH ).EQ.5    )      THEN
C       ------------------------------------------
C       This is a magnetic field harmonic
C       ------------------------------------------
          IHM          = IHM + 1
          DO IR = 1, NR
            IRMF = IR + NRIC
            INDIN  = INDFUN( IR,    IH, INARRI )
            INDOUT = INDFUN( IRMF, IHM, INARRM )
            VECM( INDOUT ) = VECIN( INDIN )
          ENDDO
C         .
C         . Now we do a simplistic linear
C         . interpolation of radial function
C         . throughout the inner core.
C         .
          IF ( NRIC.EQ.0 ) GOTO 71
          IR    = 1
          INDIN = INDFUN( IR, IH, INARRI )
          VAL   = VECIN( INDIN )
          DO IRMF  = 1, NRIC
            RAD    = XARRMF( IRMF )
            INDOUT = INDFUN( IRMF, IHM, INARRM )
            VECM( INDOUT ) = RAD*VAL/RI
          ENDDO
 71       CONTINUE
C         .
        ENDIF
C ------*        
      ENDDO
C
C Now output all the files ...
C
      FNAME(1:ILEN) = ROOT(1:ILEN)
      FNAME(ILEN+1:ILEN+6) = '.intsv'
      FNAME = FNAME(1:ILEN+6)
C
      CALL HMFWT( NHV, MHTVT, MHLVT, MHMVT, MHPVT, NDCS,
     1            MHIBCI, MHOBCI, LU, FNAME )
C
      FNAME(ILEN+1:ILEN+6) = '.xarrv'
      FNAME = FNAME(1:ILEN+6)
      CALL XARRWT( NR, XARRIN, LU, FNAME, IFORM )
C
      FNAME(ILEN+1:ILEN+6) = '.vecsv'
      FNAME = FNAME(1:ILEN+6)
      CALL SVFWT( INARRV, LU, IFORM, VECV, FNAME )
C
      FNAME(1:ILEN) = ROOT(1:ILEN)
      FNAME(ILEN+1:ILEN+6) = '.intsm'
      FNAME = FNAME(1:ILEN+6)
C
      CALL HMFWT( NHM, MHTMF, MHLMF, MHMMF, MHPMF, NDCS,
     1            MHIBCI, MHOBCI, LU, FNAME )
C
      FNAME(ILEN+1:ILEN+6) = '.xarrm'
      FNAME = FNAME(1:ILEN+6)
      CALL XARRWT( NRMF, XARRMF, LU, FNAME, IFORM )
C
      do ir = 1, nrmf
        write (6,*) ir, xarrmf( ir )
      enddo
C
      FNAME(ILEN+1:ILEN+6) = '.vecsm'
      FNAME = FNAME(1:ILEN+6)
      CALL SVFWT( INARRM, LU, IFORM, VECM, FNAME )
C
      STOP
      END
C*********************************************************************
