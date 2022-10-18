C*********************************************************************
C subroutine Convection Vector 2 NPLot format files ******************
C            -          -      - ---                ******************
C Steve Gibbons  5. 8.99                                             C
C____________________________________________________________________C
C                                                                    C
C  Takes a vector CV of dimension NDIM with NDIM = NR * NH           C
C  and the value of harmonic IH at grid node IRN is stored in        C
C  the [ ( IRN - 1 )*NH + IH ]th element of CV.                      C
C                                                                    C
C  Will output files readable by the nplot family of programs        C
C  these will be a .inters file (a converted version of the MHT,     C
C  MHL, MHM and MHC variables.                                       C
C                                                                    C
C  These will be called either .scal.inters , .mag.inters  or        C
C  .vel.inters depending upon specifications.                        C
C                                                                    C
C  The poloidal and toroidal radial functions read by the nplot      C
C  programs differ from the convection codes by a factor of R.       C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes                        C
C     NH        : Number of spherical harmonics                      C
C     NDIM      : Dimension of vector                                C
C     LU        : Logical file unit number                           C
C     IWR       : Write flag; = 2 for caution; = 3 for overwrite ..  C
C     LH        : Maximum degree spherical harmonic                  C
C                                                                    C
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C     MHM       : Dimension ( NH )                                   C
C     MHC       : Dimension ( NH )                                   C
C                                                                    C
C  mht, mhl, mhm and mhc define the spherical harmonics present      C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C   MHT( I ) = 4 for a poloidal magnetic field vector                C
C   MHT( I ) = 5 for a poloidal magnetic field vector                C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C   IFLARR      : Flag array. Dimension ( * )                        C
C                 A set of coded instructions to limit and guide     C
C                 output. Elements are as follows ...                C
C                                                                    C
C     IFLARR( 1 ) = 0 if scalars are not required.                   C
C     IFLARR( 1 ) = 1 if scalar radial functions are to be output.   C
C     IFLARR( 1 ) = 2 if .scal.inters file required in addition.     C
C                                                                    C
C     IFLARR( 2 ) = 0 if velocity vectors are not required.          C
C     IFLARR( 2 ) = 1 if velocity vectors are required.              C
C     IFLARR( 2 ) = 2 if vel. vec.s are required with .vel.inters    C
C                                                                    C
C     IFLARR( 3 ) = 0 if magnetic field vectors are not required.    C
C     IFLARR( 3 ) = 1 if magnetic field vectors are required.        C
C     IFLARR( 3 ) = 2 if mag. vec.s are required with .mag.inters    C
C                                                                    C
C     IFLARR( 4 ) = 0 if .scal.nplot file is not required.           C
C     IFLARR( 4 ) = 1 if .scal.nplot file is required.               C
C                                                                    C
C     IFLARR( 5 ) = 0 if .vel.nplot file is not required.            C
C     IFLARR( 5 ) = 1 if .vel.nplot file is required.                C
C                                                                    C
C     IFLARR( 6 ) = 0 if .mag.nplot file is not required.            C
C     IFLARR( 6 ) = 1 if .mag.nplot file is required.                C
C                                                                    C
C     IFLARR(  7 ) = IRAD   for .scal.nplot                          C
C     IFLARR(  8 ) = LAT    for .scal.nplot                          C
C     IFLARR(  9 ) = LON    for .scal.nplot                          C
C     IFLARR( 10 ) = IPLOT  for .scal.nplot                          C
C     IFLARR( 11 ) = IFN    for .scal.nplot                          C
C                                                                    C
C     IFLARR( 12 ) = IRAD   for .vel.nplot                           C
C     IFLARR( 13 ) = LAT    for .vel.nplot                           C
C     IFLARR( 14 ) = LON    for .vel.nplot                           C
C     IFLARR( 15 ) = IPLOT  for .vel.nplot                           C
C     IFLARR( 16 ) = IFN    for .vel.nplot                           C
C                                                                    C
C     IFLARR( 17 ) = IRAD   for .mag.nplot                           C
C     IFLARR( 18 ) = LAT    for .mag.nplot                           C
C     IFLARR( 19 ) = LON    for .mag.nplot                           C
C     IFLARR( 20 ) = IPLOT  for .mag.nplot                           C
C     IFLARR( 21 ) = IFN    for .mag.nplot                           C
C                                                                    C
C     IFLARR( 22 ) = IAS    Axisymmetric flag                        C
C                          = 0 for 3D and = 1 for axisym part.       C
C  Guide to IPLOT                                                    C
C  --------------                                                    C
C    iplot = 1       -->  Constant radius plot at r = IRAD           C
C    iplot = 2       -->  Constant theta plot with LAT = theta       C
C                                               IN DEGREES !!        C
C    iplot = 3       -->  Constant phi plot with LON = phi           C
C                                               IN DEGREES !!        C
C  Guide to IFN                                                      C
C  ------------                                                      C
C    streamfn{r,the,phi} = 1-3                                       C
C    vector compt{r,the,phi} = 4-6                                   C
C    vector curl{r,the,phi} = 7-9;                                   C
C    helicity = 10; |v|^2 = 11; |w|^2 = 12.                          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     CV        : Solution Vector as described above. Dim ( NDIM )   C
C     CVWORK    : Workspace. Dim ( NDIM )                            C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     ROOT      : Dimension (*). Forms the stem of all output fnames C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CV2NPL ( CV, CVWORK, NR, NH, MHT, MHL, MHM, MHC,
     1                    NDIM, RI, RO, LU, IWR, ROOT, LH, IFLARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH, NDIM, LU, IWR, LH, IFLARR( * ),
     1        MHT( NH ), MHM( NH ), MHC( NH ), MHL( NH )
      DOUBLE PRECISION RI, RO, CV( NDIM ), CVWORK( NDIM )
      CHARACTER *(*) ROOT
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L, M, ICS, IH, IRAD, NSCALH, NVELH, NMAGH, MAXL, I,
     1        INDCH, LFNAME, IICS, NT, NOH, IND1, IND2, IRN
      PARAMETER ( MAXL = 80 )
      DOUBLE PRECISION H, RAD
      CHARACTER *(MAXL) FNAME, LINE
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of NDIM against NH and NR etc.
C
      IF ( NDIM.NE.NH*NR ) THEN
        PRINT *,' Subroutine CV2NPL, bad array size'
        PRINT *,' NDIM = ',NDIM
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C Set working parameters
C
      H = (RO-RI)/DBLE( NR - 1 )
C
C Find the location in ROOT at which the ' ' occurs
C
      INDCH = 0
      DO I = 1, MAXL
        IF ( ROOT(I:I).EQ.' ' ) THEN
           INDCH = I
           GOTO 500
        ENDIF
      ENDDO
 500  CONTINUE
      IF ( INDCH.EQ.0 ) THEN
         PRINT *,' Subroutine CV2NPL, there is no space '
         PRINT *,' character in the first ',MAXL,' chars '
         PRINT *,' of ROOT. Either supply a shorter '
         PRINT *,' file name stem, or recompile with bigger '
         PRINT *,' value of MAXL.'
         STOP
      ENDIF
C
      FNAME(1:INDCH) = ROOT(1:INDCH)
C
C Add up the number of harmonics for each category
C in the solution vector. These are stored in
C NSCALH, NVELH, NMAGH for their obvious types ...
C
      NSCALH = 0
      NVELH = 0
      NMAGH = 0
C
      DO IH = 1, NH
        IF ( MHT( IH ).EQ.1 .OR. MHT( IH ).EQ.2 ) 
     1                  NVELH = NVELH + 1
        IF ( MHT( IH ).EQ.3 ) NSCALH = NSCALH + 1
        IF ( MHT( IH ).EQ.4 .OR. MHT( IH ).EQ.5 ) 
     1                  NMAGH = NMAGH + 1
      ENDDO
C
C Check value of IWR
C
      IF ( IWR.NE.2 .AND. IWR.NE.3 ) THEN
        PRINT *,' Subroutine CV2NPL: IWR = ', IWR
        STOP
      ENDIF
C
C First write out the .nplot files
C First of these is the scalar file ...
C
      IF ( IFLARR(4).EQ.1 ) THEN
         LFNAME = INDCH + 10
         FNAME(INDCH:LFNAME) = '.scal.nplot'
         CALL FOPEN  ( LU, FNAME(1:LFNAME), IWR )
C
C
C
         WRITE ( LU, 781 )
         WRITE ( LU , * ) '           TITLE'
         WRITE ( LU , * )
         WRITE ( LU, 782 )
         WRITE ( LU , * ) '           HARMONICS FILE '
         LINE(1:INDCH) = ROOT(1:INDCH)
         I = INDCH + 11
         LINE(INDCH:I) = '.scal.inters'
         WRITE ( LU, 783 ) LINE(1:I)
         WRITE ( LU, 782 )
         WRITE ( LU, * ) '           DATA FILE '
         WRITE ( LU, * )
         LINE(1:INDCH) = ROOT(1:INDCH)
         I = INDCH + 5
         LINE(INDCH:I) = '.scalf'
         WRITE ( LU, 783 ) LINE(1:I)
         WRITE ( LU , * ) 'irad    lat     lon    iplot    ifn    ias'
         WRITE ( LU , 684 ) IFLARR(  7 ), IFLARR(  8 ), IFLARR(  9 ),
     1                      IFLARR( 10 ), IFLARR( 11 ), IFLARR( 22 )
         WRITE ( LU, 782 )
         WRITE ( LU , * ) 'lh     nh     nr    inorm      ri        ro'
         WRITE ( LU , 784 ) LH, NSCALH, NR, RI, RO, RI/RO
         WRITE ( LU, 782 )
         WRITE (LU, *) 'levels      colours     contours     labels',
     1                '     colourbar     edges'
         WRITE (LU, *) '  12           0           1            0  ',
     1                '         0           0  '
         WRITE ( LU, 782 )
         WRITE ( LU, * ) ' info man_cont  cont_min  cont_max '
         WRITE ( LU, * ) '   0      0       10.0      10.0   '
         WRITE ( LU, 782 )
         WRITE ( LU, * ) 'iplot  =  projection : radial = 1; '
         WRITE ( LU, * ) ' equatorial = 2; meridional = 3. '
         WRITE ( LU, * ) 'ifn    =  fn plotted : '
         WRITE ( LU, * ) 'streamfn{r,the,phi} = 1-3;'
         WRITE ( LU, * ) 'vector compt{r,the,phi} = 4-6; '
         WRITE ( LU, * ) 'vector curl{r,the,phi} = 7-9; '
         WRITE ( LU, * ) '     helicity = 10; |v|^2 = 11; |w|^2 = 12. '
         WRITE ( LU, 781 )
C
C
C
         CALL FCLOSE ( LU, FNAME(1:LFNAME), 'error' )
      ENDIF
C
C Now write out .vel.nplot
C
      IF ( IFLARR(5).EQ.1 ) THEN
         LFNAME = INDCH + 9
         FNAME(INDCH:LFNAME) = '.vel.nplot'
         CALL FOPEN  ( LU, FNAME(1:LFNAME), IWR )
C
C
C
         WRITE ( LU, 781 )
         WRITE ( LU , * ) '           TITLE'
         WRITE ( LU , * )
         WRITE ( LU, 782 )
         WRITE ( LU , * ) '           HARMONICS FILE '
         LINE(1:INDCH) = ROOT(1:INDCH)
         I = INDCH + 10
         LINE(INDCH:I) = '.vel.inters'
         WRITE ( LU, 783 ) LINE(1:I)
         WRITE ( LU, 782 )
         WRITE ( LU, * ) '           DATA FILE '
         LINE(1:INDCH) = ROOT(1:INDCH)
         I = INDCH + 9
         LINE(INDCH:I) = '.vel.evecs'
         WRITE ( LU, 783 ) LINE(1:I)
         LINE(1:INDCH) = ROOT(1:INDCH)
         I = INDCH + 10
         LINE(INDCH:I) = '.scal.evecs'
         WRITE ( LU, 783 ) LINE(1:I)
         WRITE ( LU , * ) 'irad    lat     lon    iplot    ifn    ias'
         WRITE ( LU , 684 ) IFLARR( 12 ), IFLARR( 13 ), IFLARR( 14 ),
     1                      IFLARR( 15 ), IFLARR( 16 ), IFLARR( 22 )
         WRITE ( LU, 782 )
         WRITE ( LU , * ) 'lh     nh     nr    inorm      ri        ro'
         WRITE ( LU , 784 ) LH, NVELH, NR, RI, RO, RI/RO
         WRITE ( LU, 782 )
         WRITE (LU, *) 'levels      colours     contours     labels',
     1                '     colourbar     edges'
         WRITE (LU, *) '  12           0           1            0  ',
     1                '         0           0  '
         WRITE ( LU, 782 )
         WRITE ( LU, * ) ' info man_cont  cont_min  cont_max '
         WRITE ( LU, * ) '   0      0       10.0      10.0   '
         WRITE ( LU, 782 )
         WRITE ( LU, * ) 'iplot  =  projection : radial = 1; '
         WRITE ( LU, * ) ' equatorial = 2; meridional = 3. '
         WRITE ( LU, * ) 'ifn    =  fn plotted : '
         WRITE ( LU, * ) 'streamfn{r,the,phi} = 1-3;'
         WRITE ( LU, * ) 'vector compt{r,the,phi} = 4-6; '
         WRITE ( LU, * ) 'vector curl{r,the,phi} = 7-9; '
         WRITE ( LU, * ) '     helicity = 10; |v|^2 = 11; |w|^2 = 12. '
         WRITE ( LU, 781 )
C
C
C
         CALL FCLOSE ( LU, FNAME(1:LFNAME), 'error' )
      ENDIF         
C
C Now output .mag.nplot
C
      IF ( IFLARR(6).EQ.1 ) THEN
         LFNAME = INDCH + 9
         FNAME(INDCH:LFNAME) = '.mag.nplot'
         CALL FOPEN  ( LU, FNAME(1:LFNAME), IWR )
C
C
C
         WRITE ( LU, 781 )
         WRITE ( LU , * ) '           TITLE'
         WRITE ( LU , * )
         WRITE ( LU, 782 )
         WRITE ( LU , * ) '           HARMONICS FILE '
         LINE(1:INDCH) = ROOT(1:INDCH)
         I = INDCH + 10
         LINE(INDCH:I) = '.mag.inters'
         WRITE ( LU, 783 ) LINE(1:I)
         WRITE ( LU, 782 )
         WRITE ( LU, * ) '           DATA FILE '
         LINE(1:INDCH) = ROOT(1:INDCH)
         I = INDCH + 9
         LINE(INDCH:I) = '.mag.evecs'
         WRITE ( LU, 783 ) LINE(1:I)
         LINE(1:INDCH) = ROOT(1:INDCH)
         I = INDCH + 10
         LINE(INDCH:I) = '.scal.evecs'
         WRITE ( LU, 783 ) LINE(1:I)
         WRITE ( LU , * ) 'irad    lat     lon    iplot    ifn    ias'
         WRITE ( LU , 684 ) IFLARR( 17 ), IFLARR( 18 ), IFLARR( 19 ),
     1                      IFLARR( 21 ), IFLARR( 22 ), IFLARR( 22 )
         WRITE ( LU, 782 )
         WRITE ( LU , * ) 'lh     nh     nr    inorm      ri        ro'
         WRITE ( LU , 784 ) LH, NMAGH, NR, RI, RO, RI/RO
         WRITE ( LU, 782 )
         WRITE (LU, *) 'levels      colours     contours     labels',
     1                '     colourbar     edges'
         WRITE (LU, *) '  12           0           1            0  ',
     1                '         0           0  '
         WRITE ( LU, 782 )
         WRITE ( LU, * ) ' info man_cont  cont_min  cont_max '
         WRITE ( LU, * ) '   0      0       10.0      10.0   '
         WRITE ( LU, 782 )
         WRITE ( LU, * ) 'iplot  =  projection : radial = 1; '
         WRITE ( LU, * ) ' equatorial = 2; meridional = 3. '
         WRITE ( LU, * ) 'ifn    =  fn plotted : '
         WRITE ( LU, * ) 'streamfn{r,the,phi} = 1-3;'
         WRITE ( LU, * ) 'vector compt{r,the,phi} = 4-6; '
         WRITE ( LU, * ) 'vector curl{r,the,phi} = 7-9; '
         WRITE ( LU, * ) '     helicity = 10; |v|^2 = 11; |w|^2 = 12. '
         WRITE ( LU, 781 )
C
C
C
         CALL FCLOSE ( LU, FNAME(1:LFNAME), 'error' )
      ENDIF
C
C Now output .scal.inters
C
      IF ( IFLARR( 1 ).EQ.2 ) THEN
         IF ( NSCALH.EQ.0 ) GOTO 401
C
         LFNAME = INDCH + 11
         FNAME(INDCH:LFNAME) = '.scal.inters'
         CALL FOPEN  ( LU, FNAME(1:LFNAME), IWR )
C
         NT = 1
         DO IH = 1, NH
           IF ( MHT( IH ).EQ.3 .AND. MHL( IH ).NE.0 ) THEN
             L    = MHL( IH )
             M    = MHM( IH )
             IICS = MHC( IH )
             ICS  = 3 - IICS
C note nplot uses ics =1 for my ics=2 and vice versa
             WRITE ( LU, 791 )  NT, L, M, ICS
           ENDIF
         ENDDO
C
C
C
         CALL FCLOSE ( LU, FNAME(1:LFNAME), 'error' )
C
 401  CONTINUE
      ENDIF
C
C Now output .vel.inters
C
      IF ( IFLARR( 2 ).EQ.2 ) THEN
         IF ( NVELH.EQ.0 ) GOTO 402
C
         LFNAME = INDCH + 10
         FNAME(INDCH:LFNAME) = '.vel.inters'
         CALL FOPEN  ( LU, FNAME(1:LFNAME), IWR )
C
         NT = 2
         DO IH = 1, NH
           IF ( MHT( IH ).EQ.2 ) THEN
             L    = MHL( IH )
             M    = MHM( IH )
             IICS = MHC( IH )
             ICS  = 3 - IICS
C note nplot uses ics =1 for my ics=2 and vice versa
             WRITE ( LU, 791 )  NT, L, M, ICS
           ENDIF
         ENDDO
C
         NT = 1
         DO IH = 1, NH
           IF ( MHT( IH ).EQ.1 ) THEN
             L    = MHL( IH )
             M    = MHM( IH )
             IICS = MHC( IH )
             ICS  = 3 - IICS
C note nplot uses ics =1 for my ics=2 and vice versa
             WRITE ( LU, 791 )  NT, L, M, ICS
           ENDIF
         ENDDO
 402  CONTINUE
      ENDIF
C
C Now output .mag.inters
C
      IF ( IFLARR( 3 ).EQ.2 ) THEN
         IF ( NMAGH.EQ.0 ) GOTO 403
C
         LFNAME = INDCH + 10
         FNAME(INDCH:LFNAME) = '.mag.inters'
         CALL FOPEN  ( LU, FNAME(1:LFNAME), IWR )
C
         NT = 2
         DO IH = 1, NH
           IF ( MHT( IH ).EQ.5 ) THEN
             L    = MHL( IH )
             M    = MHM( IH )
             IICS = MHC( IH )
             ICS  = 3 - IICS
C note nplot uses ics =1 for my ics=2 and vice versa
             WRITE ( LU, 791 )  NT, L, M, ICS
           ENDIF
         ENDDO
C
         NT = 1
         DO IH = 1, NH
           IF ( MHT( IH ).EQ.4 ) THEN
             L    = MHL( IH )
             M    = MHM( IH )
             IICS = MHC( IH )
             ICS  = 3 - IICS
C note nplot uses ics =1 for my ics=2 and vice versa
             WRITE ( LU, 791 )  NT, L, M, ICS
           ENDIF
         ENDDO
C
C
         CALL FCLOSE ( LU, FNAME(1:LFNAME), 'error' )
C
 403  CONTINUE
      ENDIF
C
C Now write out .scalf file
C
      IF ( IFLARR( 1 ).EQ.1 .OR. IFLARR( 1 ).EQ.2 .AND.
     1                  NSCALH.GT.0          ) THEN
C
C Let's rearrange the elements into CVWORK array
C         
        NOH = 0
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.3 .AND. MHL( IH ).EQ.0 ) THEN
            NOH = NOH + 1
            DO IRN = 1, NR
              IND1 = ( IRN - 1 )*NH + IH
              IND2 = ( NOH - 1 )*NR + IRN
              CVWORK( IND2 ) = CV( IND1 )
            ENDDO
          ENDIF
        ENDDO
C
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.3 .AND. MHL( IH ).NE.0 ) THEN
            NOH = NOH + 1
            DO IRN = 1, NR
              IND1 = ( IRN - 1 )*NH + IH
              IND2 = ( NOH - 1 )*NR + IRN
              CVWORK( IND2 ) = CV( IND1 )
            ENDDO
          ENDIF
        ENDDO
C         
        IF ( NOH.NE.NSCALH ) THEN
           PRINT *,' Subroutine CV2NPL, NOH = ', NOH
           PRINT *,' NSCALH = ', NSCALH,'. Error.'
           STOP
        ENDIF
C         
        NOH = NOH*NR
C         
        LFNAME = INDCH + 6
        FNAME(INDCH:LFNAME) = '.scalf'
        CALL FOPEN  ( LU, FNAME(1:LFNAME), IWR )
        WRITE ( LU, 100 )
        WRITE ( LU, * ) '  0.0  0.0  0.0  0.0  0.0 '
        WRITE ( LU, 373 ) LH, NSCALH-1, NR, RI/RO
        WRITE (LU,792) (CVWORK(I),I=1,NOH)
        CALL FCLOSE ( LU, FNAME(1:LFNAME), 'error' )
C         
      ENDIF
C
C Now write out .vel.evecs file
C
      IF ( IFLARR( 2 ).EQ.1 .OR. IFLARR( 2 ).EQ.2 .AND.
     1                  NVELH.GT.0          ) THEN
C
C Let's rearrange the elements into CVWORK array
C
        NOH = 0
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.2 ) THEN
            NOH = NOH + 1
            DO IRN = 1, NR
              RAD = RI + DBLE( IRN - 1 )*H
              IND1 = ( IRN - 1 )*NH + IH
              IND2 = ( NOH - 1 )*NR + IRN
              CVWORK( IND2 ) = CV( IND1 )*RAD
            ENDDO
          ENDIF
        ENDDO
C
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.1 ) THEN
            NOH = NOH + 1
            DO IRN = 1, NR
              RAD = RI + DBLE( IRN - 1 )*H
              IND1 = ( IRN - 1 )*NH + IH
              IND2 = ( NOH - 1 )*NR + IRN
              CVWORK( IND2 ) = CV( IND1 )*RAD
            ENDDO
          ENDIF
        ENDDO
C
        IF ( NOH.NE.NVELH ) THEN
           PRINT *,' Subroutine CV2NPL, NOH = ', NOH
           PRINT *,' NVELH = ', NVELH,'. Error.'
           STOP
        ENDIF
C
        NOH = NOH*NR
C
        LFNAME = INDCH + 9
        FNAME(INDCH:LFNAME) = '.vel.evecs'
        CALL FOPEN  ( LU, FNAME(1:LFNAME), IWR )
        WRITE ( LU, 100 )
        WRITE ( LU, * ) '  0.0  0.0  0.0  0.0  0.0 '
        WRITE ( LU, 373 ) LH, NVELH, NR, RI/RO
        WRITE (LU,792) (CVWORK(I),I=1,NOH)
        CALL FCLOSE ( LU, FNAME(1:LFNAME), 'error' )
C
      ENDIF
C
C Now write out .mag.evecs file
C
      IF ( IFLARR( 3 ).EQ.1 .OR. IFLARR( 3 ).EQ.2 .AND.
     1                  NMAGH.GT.0          ) THEN
C
C Let's rearrange the elements into CVWORK array
C
        NOH = 0
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.5 ) THEN
            NOH = NOH + 1
            DO IRN = 1, NR
              RAD = RI + DBLE( IRN - 1 )*H
              IND1 = ( IRN - 1 )*NH + IH
              IND2 = ( NOH - 1 )*NR + IRN
              CVWORK( IND2 ) = CV( IND1 )*RAD
            ENDDO  
          ENDIF
        ENDDO
C
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.4 ) THEN
            NOH = NOH + 1
            DO IRN = 1, NR
              RAD = RI + DBLE( IRN - 1 )*H
              IND1 = ( IRN - 1 )*NH + IH
              IND2 = ( NOH - 1 )*NR + IRN
              CVWORK( IND2 ) = CV( IND1 )*RAD
            ENDDO
          ENDIF
        ENDDO
C
        IF ( NOH.NE.NMAGH ) THEN
           PRINT *,' Subroutine CV2NPL, NOH = ', NOH
           PRINT *,' NMAGH = ', NMAGH,'. Error.'
           STOP
        ENDIF   
C
        NOH = NOH*NR
C
        LFNAME = INDCH + 9
        FNAME(INDCH:LFNAME) = '.mag.evecs'
        CALL FOPEN  ( LU, FNAME(1:LFNAME), IWR )
        WRITE ( LU, 100 )
        WRITE ( LU, * ) '  0.0  0.0  0.0  0.0  0.0 '
        WRITE ( LU, 373 ) LH, NMAGH, NR, RI/RO
        WRITE (LU,792) (CVWORK(I),I=1,NOH)
        CALL FCLOSE ( LU, FNAME(1:LFNAME), 'error' )
C
      ENDIF
C
C
C
      RETURN
 373  FORMAT (i6,i6,i6,' 0.0  0.0  ',f10.4,' 0.0  0.0 ')
 100  FORMAT (' ')
 781  FORMAT ('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx',
     1        'xxxxxxxxxxxxxxxxxxxx')
 782  FORMAT ('----------------------------------------------------',
     1        '--------------------')
 783  FORMAT (a)
 784  FORMAT (i3,'    ',i3,'    ',i4,'     1    ',f6.2,'     ',f6.2,
     1        '    ',f6.2)
 684  FORMAT (i7,i7,i7,i7,i7,i7)
 791  FORMAT (4i2)
 792  FORMAT ( 5e15.7 )
      END
C*********************************************************************

