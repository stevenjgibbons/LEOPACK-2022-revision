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
C subroutine HarMonic File ReaD **************************************
C            -  -     -    -  - **************************************
C Steve Gibbons Sat Nov 13 12:50:33 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Reads in the integer indices of spherical harmonic sets incl.      C
C the appropriate boundary conditions from a file.                   C
C It carefully assesses the boundary conditions and creates or       C
C appends the MHIBC, MHOBC and LARR arrays which must be sent to     C
C SVFDCF to calculate the finite difference coefficients.            C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Number of vector spherical harmonics. (Output)     C
C     NHMAX     : Maximum number of vector spherical harmonics.      C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NHMAX          C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NHMAX          C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NHMAX          C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NHMAX          C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     NCUDS     : Number of currently used finite diff. schemes.     C
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
C     LARR      : Spherical harmonic degree, L. Dim. ( NDCS ).       C
C                 Array to be passed to SVFDCF.                      C
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
      SUBROUTINE HMFRD( NH, NHMAX, MHT, MHL, MHM, MHP, NCUDS, NDCS,
     1                  MHIBC, MHOBC, LARR, LU, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, NHMAX, MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     1        MHP( NHMAX ), NCUDS, NDCS,
     1        MHIBC( NDCS ), MHOBC( NDCS ), LARR( NDCS ), LU
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH, IWR, IS, IIBF, IOBF, ND1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
C Read in number of spherical harmonics
C
      READ ( LU, * ) NH
C
C Check value of NH
C
      IF ( NH.GT.NHMAX ) THEN
        PRINT *, ' Subroutine HMFRD.'
        PRINT *, ' From file, NH = ', NH
        PRINT *, ' NHMAX = ', NHMAX
        PRINT *, ' Program aborted.'
        STOP
      ENDIF
C
C     Now loop around the harmonics and read in their
C     properties
C
      DO IH = 1, NH
        READ ( LU, * ) MHT( IH ), MHL( IH ), MHM( IH ), IIBF, IOBF
C       .
C       . We now need to see whether we already have a
C       . suitable finite difference scheme for this outcome.
C       .
        ND1 = MIN( NCUDS, NDCS )
        DO IS = 1, ND1
         IF ( MHIBC( IS ).EQ.IIBF .AND. MHOBC( IS ).EQ.IOBF ) THEN
C          .
C          . Simple case where LARR needn't be applied
C          .
           IF ( IIBF.NE.7 .AND. IOBF.NE.7 ) THEN
             MHP( IH ) = IS
             GOTO 60
           ENDIF
C          .
C          . OK - check to see if LARR value is correct
C          .
           IF ( LARR( IS ).EQ.MHL( IH ) ) THEN
             MHP( IH ) = IS
             GOTO 60
           ENDIF
C          .
         ENDIF
        ENDDO
C       .
C       . OK - we don't already have this diff. scheme
C       .
        NCUDS = NCUDS + 1
        IF ( NCUDS.GT.NDCS ) THEN
          PRINT *, ' Subroutine HMFRD.'
          PRINT *, NDCS,' schemes are available.'
          PRINT *, ' and these have now been exhausted.'
          PRINT *, ' Program aborted.'
          STOP
        ENDIF
C       .
        MHIBC( NCUDS ) = IIBF
        MHOBC( NCUDS ) = IOBF
        IF ( IIBF.NE.7 .AND. IOBF.NE.7 ) THEN
          LARR( NCUDS ) = 0
        ELSE
          LARR( NCUDS ) = MHL( IH )
        ENDIF
        MHP( IH ) = NCUDS
C       .
 60     CONTINUE
C       .
      ENDDO
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
C Abort program if more finite difference schemes were
C requested by the harmonic sets than were allowed
C
      IF ( NCUDS.LT.NDCS ) THEN
        DO IS = NCUDS + 1, NDCS
          LARR( IS ) = -1
        ENDDO
      ENDIF
C
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine Solution Vector File ReaD *******************************
C            -        -      -    -  - *******************************
C Steve Gibbons Sat Nov 13 14:30:24 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Reads in a solution vector from a file.                            C
C IMPORTANT. The set of spherical harmonics must ALREADY BE KNOWN    C
C before calling SVFRD. If SVFRD reads in a file with NH different   C
C to the NH input in INARR( 3 ), an error will be reported.          C
C                                                                    C
C The only check on NR is that it is not greater than NRMAX.         C
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
C     NRMAX     : Maximum permitted radial grid nodes.               C
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
      SUBROUTINE SVFRD( INARR, LU, NRMAX, SV, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), LU, NRMAX
      CHARACTER *(*) FNAME
      DOUBLE PRECISION SV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IWR, ILEN, IFORMF, NR, NH, NH1, IFORM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NH     = INARR( 3 )
C
C Open file for reading
C
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
C Read iformf, nr, nh, iform
C  
       READ ( LU, * ) IFORMF, NR, NH1, IFORM
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine SVFRD.'
        PRINT *,' IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check that value of NH is legal
C
      IF ( NH.NE.NH1 ) THEN
        PRINT *,' Subroutine SVFRD.'
        PRINT *,' INARR( 3 ) = ', NH
        PRINT *,' File contains ',NH1,' harmonics.'
        PRINT *,' Load in correct indices before calling SVFRD.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check that value of NR is legal
C
      IF ( NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine SVFRD.'
        PRINT *,' In file, NR = ', NR
        PRINT *,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INARR( 1 ) = IFORMF
      INARR( 2 ) = NR
C
      ILEN   = NR*NH
C
C OK, so read SV values ...
C
      IF ( IFORM.EQ.1 ) READ ( LU, 41 ) ( SV( I ), I = 1, ILEN )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 41   FORMAT(5(1PD16.7))
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine X value ARRay ReaD **************************************
C            -       ---   -  - **************************************
C Steve Gibbons Fri Nov 12 07:58:51 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Reads in the XARR array of abscissae from a file.                  C
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
C     NRMAX     : Maximum number of radial grid nodes.               C
C                                                                    C
C     LU        : Logical file unit number.                          C
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
      SUBROUTINE XARRRD( NR, NRMAX, XARR, LU, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NRMAX, LU
      CHARACTER *(*) FNAME
      DOUBLE PRECISION XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IWR, IFORM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
C Read number of radial grid nodes
C
      READ ( LU, * ) NR, IFORM
C
C Check NR
C
      IF ( NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine XARRRD.'
        PRINT *,' From file, NR = ', NR
        PRINT *,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine XARRRD.'
        PRINT *,' From file, IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C OK, so read X values ...
C
      IF ( IFORM.EQ.1 ) READ ( LU, 41 ) ( XARR( I ), I = 1, NR )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 41   FORMAT(5(1PD16.7))
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
C subroutine Zeroes of Chebyshev Polynomial Abscissa Allocation 2 ****
C            -         -         -          -        -          - ****
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
      SUBROUTINE ZCPAA2( NR, XARR, RI, RO)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR
      DOUBLE PRECISION RI, RO, XARR( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER N, I
      DOUBLE PRECISION PH, DN, C, SI, DI, CSX, DM, DC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NR.LT.4 ) THEN
         PRINT *,' Subroutine ZCPAA2.'
         PRINT *,' NR = ', NR
         STOP
      ENDIF
C
      IF ( RI.GE.RO ) THEN
         PRINT *,' Subroutine ZCPAA2.'
         PRINT *,' RI = ', RI
         PRINT *,' RO = ', RO
         STOP
      ENDIF
C
      DM = RO - RI
      DC = 2.0d0*RI - RO
C
      XARR(  1 ) = RI
      XARR( NR ) = RO
      N = NR - 2
      XARR(2) = 0.D0
      PH = 1.57079632679489661923D0
      DN = DFLOAT(2*N)
      C  = PH/DN
      SI = -1.D0
      DO 10 I = 1, N
         DI = DFLOAT(I)
         CSX = DCOS(C*(2.D0*DI-1.D0))
         XARR(N-I+2) = (1.0d0+CSX)*DM + DC
         SI = -SI
 10   CONTINUE
C
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
