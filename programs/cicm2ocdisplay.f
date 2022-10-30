C*********************************************************************
C                                                                    C
C                 Conducting Inner Core and Mantle                   C
C Steve Gibbons - 2 Outer Core DISPLAY solution                      C
C Sat Feb 22 12:07:03 MET 2003                                       C
C                                                                    C
C This program reads in the same files as would be read in           C
C for cicmibcdts2.exe ( XXX.xarrm, XXX.vecsm, XXX.intsm, XXX.xarrv,  C
C XXX.vecsv, XXX.intsv)                                              C
C                                                                    C
C It then reads in SCAL, the strength of an inhomogeneous            C
C temperature. If SCAL is non-zero, it proceeds to read in the       C
C ( .ints, .vecs and .xarr files used to store the temperature).     C
C                                                                    C
C This program simply writes out a single .ints, .vecs and .xarr     C
C file corresponding to the outer core region (with inhomogeneous    C
C temperature superimposed if necessary).                            C
C                                                                    C
C Note that the output from this program is STRICTLY for display     C
C purposes only as the boundary conditions will not be meaningful.   C
C                                                                    C
C                                                                    C
C*********************************************************************
      PROGRAM cicm2ocdisplay
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          NRVMAX, NRMMAX, NH1MAX, NH2MAX, NH3MAX,
     1                 NH4MAX, NH5MAX, NHVMAX, NHMMAX, NIVMAX, NIMMAX,
     2                 NHMAX, NIMAX, NNDS,
     3                 NIVCMX, NHCMAX, NRMAXC, NDCS
      PARAMETER      ( NRVMAX = 40, NRMMAX = 60, NH1MAX = 200,
     1                 NH2MAX = 200, NH3MAX = 200, NH4MAX = 200,
     2                 NH5MAX = 200, NHVMAX = NH1MAX+NH2MAX+NH3MAX,
     3                 NHMMAX = NH4MAX+NH5MAX, NNDS = 4,
     4                 NIVMAX = NHVMAX*NRVMAX, NHMAX = NHVMAX+NHMMAX,
     5                 NIMMAX = NHMMAX*NRMMAX, NIMAX = NIVMAX+NIMMAX )
      PARAMETER      ( NRMAXC = 100,  NDCS = 1,
     4                 NHCMAX = 300, NIVCMX = NHCMAX*NRMAXC )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER          MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     1                 MHP( NHMAX ), MHIBC( NDCS ), MHOBC( NDCS )
      INTEGER          MHTV( NHVMAX ), MHLV( NHVMAX ), MHMV( NHVMAX )
      INTEGER          MHTM( NHMMAX ), MHLM( NHMMAX ), MHMM( NHMMAX )
      INTEGER          MHTC( NHCMAX ), MHLC( NHCMAX ), MHMC( NHCMAX )
C
      DOUBLE PRECISION SVC( NIVCMX ), XARRC( NRMAXC )
      DOUBLE PRECISION SVV( NIVMAX ), XARRV( NRVMAX )
      DOUBLE PRECISION SVM( NIMMAX ), XARRM( NRMMAX )
      DOUBLE PRECISION SV( NIMAX )
C
      INTEGER          IWORK( NNDS )
      DOUBLE PRECISION WORK1( NNDS ), WORK2( NNDS ),
     1                 WORKM( NNDS, NNDS )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER           NHV, NHM, NHC, NH, ILEN, ILENRT, LU,
     1                  I, NRV, NRM, NRC, NR, INARRV( 3 ),
     2                  INARRM( 3 ), INARRC( 3 ), INARR( 3 ),
     3                  NRCM, NRIC, NR1, IR, IRV, IRM, INDFUN,
     4                  IND, INDV, INDM, IFORM
      INTEGER           IHARM, IH, IHC, IHM, IHV
      DOUBLE PRECISION  DLOW, SCAL, RAD, RI, RO, RIC, ROC
      CHARACTER *(120)  LINE, FNXARR, FNINTS, FNVECS, ROOT, FNAME
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      EXTERNAL          INDFUN
      PARAMETER      ( DLOW = 1.0d-8 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      LU      = 11
 80   FORMAT(A)
C
      DO I = 1, 120
        LINE(I:I) = ' '
      ENDDO
C
      PRINT *,' Enter stem for output files.'
 200  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 200
C
C Extract ROOT from LINE
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILENRT = I - 1
          GOTO 199
        ENDIF
      ENDDO
 199  CONTINUE
      ROOT(1:ILENRT)   = LINE(1:ILENRT)
      FNXARR(1:ILENRT) = ROOT(1:ILENRT)
      FNXARR(ILENRT+1:ILENRT+5) = '.xarr'
      FNXARR = FNXARR(1:ILENRT+5)
C
      FNINTS(1:ILENRT) = ROOT(1:ILENRT)
      FNINTS(ILENRT+1:ILENRT+5) = '.ints'
      FNINTS = FNINTS(1:ILENRT+5)
C
      FNVECS(1:ILENRT) = ROOT(1:ILENRT)
      FNVECS(ILENRT+1:ILENRT+5) = '.vecs'
      FNVECS = FNVECS(1:ILENRT+5)
C
      PRINT *,' Enter name of velocity harmonics file.'
 201  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 201
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 202
        ENDIF
      ENDDO
 202  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in velocity harmonics file
C
      CALL BIHFRD( NHV, NHVMAX, MHTV, MHLV, MHMV, LU, FNAME )
      INARRV( 3 ) = NHV
C
      PRINT *,' Enter name of velocity vector file.'
 204  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 204
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 205
        ENDIF
      ENDDO
 205  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in velocity vector file
C
      CALL SVFRD( INARRV, LU, NRVMAX, SVV, FNAME )
      NR1 = INARRV( 2 )
C
      PRINT *,' Enter name of velocity radial spacing file.'
 206  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 206
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 207
        ENDIF
      ENDDO
 207  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in velocity radial spacing file
C
      CALL XARRRD( NRV, NRVMAX, XARRV, LU, FNAME )
C
      IF ( NR1.NE.NRV ) THEN
        WRITE ( 6, * ) 'Vel. soln. vector and radial node '
        WRITE ( 6, * ) 'file claim differing numbers of '
        WRITE ( 6, * ) 'grid nodes. Program aborted.'
        STOP
      ENDIF
C
      PRINT *,' Enter name of mag. harmonics file.'
 401  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 401
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 402
        ENDIF
      ENDDO
 402  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in magnetic harmonics file
C
      CALL BIHFRD( NHM, NHMMAX, MHTM, MHLM, MHMM, LU, FNAME )
      INARRM( 3 ) = NHM
C
      PRINT *,' Enter name of magnetic field vector file.'
 404  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 404
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I - 1
          GOTO 405
        ENDIF
      ENDDO
 405  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in magnetic field vector file
C
      CALL SVFRD( INARRM, LU, NRMMAX, SVM, FNAME )
      NR1 = INARRM( 2 )
C
      PRINT *,' Enter name of magnetic field radial spacing file.'
 406  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 406
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 407
        ENDIF
      ENDDO
 407  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in magnetic field radial spacing file
C
      CALL XARRRD( NRM, NRMMAX, XARRM, LU, FNAME )
C
      IF ( NR1.NE.NRM ) THEN
        WRITE ( 6, * ) 'Mag. soln. vector and radial node '
        WRITE ( 6, * ) 'file claim differing numbers of '
        WRITE ( 6, * ) 'grid nodes. Program aborted.'
        STOP
      ENDIF
C
C Ensure that our radial nodes for velocity
C and magnetic field are compatible
C
      CALL XARRC2( NRV, XARRV, NRM, XARRM, NRIC )
      NRCM = NRM - NRV - NRIC
C
      WRITE ( 6, * ) ' There are ',NRV,' nodes for outer core.'
      WRITE ( 6, * ) ' There are ',NRM,' nodes for magnetic field.'
      WRITE ( 6, * ) ' There are ',NRIC,' nodes for inner core.'
      WRITE ( 6, * ) ' There are ',NRCM,' nodes for the mantle.'
C     .
      PRINT *,' Enter scaling factor for temperature inhomog.'
 412  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 412
C
      READ ( LINE, * ) SCAL
      WRITE ( 6, * ) 'SCAL  = ', SCAL
      WRITE ( 6, * )
      IF ( DABS( SCAL ).LT.DLOW ) GOTO 600
C
C Skip the next part if our inhomogeneous part is zero
C
      PRINT *,' Enter name of inhomog. temp. harmonics file.'
 506  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 506
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 507
        ENDIF
      ENDDO
 507  CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in harmonics file (ignoring boundary conditions!)
C
      CALL BIHFRD( NHC, NHCMAX, MHTC, MHLC, MHMC, LU, FNAME )
C
      PRINT *,' Enter name of inhomog. temp. vector file.'
 408  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 408
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 409
        ENDIF
      ENDDO
 409  CONTINUE
      FNAME = LINE(1:ILEN)
C
      INARRC( 3 ) = NHC
      CALL SVFRD( INARRC, LU, NRMAXC, SVC, FNAME )
C
      PRINT *,' Enter name of inhomog. radial spacings file.'
 410  CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 410
C
      DO I = 1, 120
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 411
        ENDIF
      ENDDO
 411  CONTINUE
      FNAME = LINE(1:ILEN)
C
      CALL XARRRD( NRC, NRMAXC, XARRC, LU, FNAME )
C
      IF ( INARRC( 2 ).NE.NRC ) THEN
        PRINT *,' Radial vector file and spacings file'
        PRINT *,' are incompatible. Program aborted.'
        STOP
      ENDIF
C
C All files appear to be read in now
C Fill in the index arrays for the new solution vector.
C
 600  CONTINUE
      NH     = NHM + NHV
      NR     = NRV
      INARR( 1 ) = INARRV( 1 )
      INARR( 2 ) = NR
      INARR( 3 ) = NH
      PRINT *,' NHM = ', NHM
      PRINT *,' NHV = ', NHV
      DO IHV = 1, NHV
        IH        = IHV
        MHL( IH ) = MHLV( IHV )
        MHM( IH ) = MHMV( IHV )
        MHT( IH ) = MHTV( IHV )
        MHP( IH ) = 1
        DO IRV = 1, NRV
          IR      = IRV
          INDV    = INDFUN(  IRV, IHV, INARRV )
          IND     = INDFUN(  IR , IH , INARR  )
          SV( IND ) = SVV( INDV )
        ENDDO
      ENDDO
      DO IHM = 1, NHM
        IH        = IHM + NHV
        MHL( IH ) = MHLM( IHM )
        MHM( IH ) = MHMM( IHM )
        MHT( IH ) = MHTM( IHM )
        MHP( IH ) = 1
        DO IRV = 1, NRV
          IRM     = IRV + NRIC
          IR      = IRV
          INDM    = INDFUN(  IRM, IHM, INARRM )
          IND     = INDFUN(  IR , IH , INARR  )
          SV( IND ) = SVM( INDM )
        ENDDO
      ENDDO
C
C Now we need to add the inhomogeneous temperature
C functions if they are required.
C
      IF ( DABS( SCAL ).GT.DLOW ) THEN
C       .
C       . Check the inner and outer core radii
C       .
        RI       = XARRV(   1 )
        RO       = XARRV(  NR )
C       .
        RIC      = XARRC(   1 )
        ROC      = XARRC( NRC )
C       .
        IF ( DABS( RI-RIC ).GT.DLOW  .OR.
     1       DABS( RO-ROC ).GT.DLOW      ) THEN
          PRINT *,' RI = ', RI,' RIC = ', RIC
          PRINT *,' RO = ', RO,' ROC = ', ROC
          PRINT *,' Incompatible radial functions.'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
        DO IHC = 1, NHC
          IF ( MHTC( IHC ).NE.3 ) GOTO 50
C         .
C         . OK so IHC is a temperature harmonic -
C         . so let's try to find the corresponding harmonic.
C         .
          IHARM = 0
          DO IH = 1, NH
            IF (     MHT( IH ).EQ.3             .AND.
     1               MHL( IH ).EQ.MHLC( IHC )   .AND.
     2               MHM( IH ).EQ.MHMC( IHC )        ) IHARM = IH
C           .
          ENDDO
          IF ( IHARM.EQ.0 ) THEN
            PRINT *,' Solution vector contains no harmonic, ih,'
            PRINT *,' with MHL( ih ) = ', MHLC( IHC )
            PRINT *,' and MHM( ih ) = ', MHMC( IHC )
            PRINT *,' Program aborted.'
            STOP
          ENDIF
C         .
          DO IR = 1, NR
            IND    = INDFUN(  IR, IHARM, INARR )
            RAD    = XARRC( IR )
            CALL SVRINT( RAD, SVC, XARRC, INARRC, IHC, NNDS, WORK1,
     1                   IWORK, WORK2, WORKM )
C           .
            SV( IND ) = SV( IND ) + SCAL*WORK1( 1 )
C           .
          ENDDO
 50     CONTINUE
        ENDDO
C       .
      ENDIF
C
C Now ready to output all of these files
C
      MHIBC( 1 ) = 1
      MHOBC( 1 ) = 1
      IFORM      = 1
      CALL HMFWT( NH, MHT, MHL, MHM, MHP, NDCS,
     1            MHIBC, MHOBC, LU, FNINTS )
      CALL XARRWT( NR, XARRV, LU, FNXARR, IFORM )
      CALL SVFWT( INARR, LU, IFORM, SV, FNVECS )
C
      STOP
      END
C*********************************************************************
