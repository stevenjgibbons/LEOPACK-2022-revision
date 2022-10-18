C*********************************************************************
C subroutine Imposed Temperature Vector Auxiliary Arrays Form ********
C            -       -           -      -         -      -    ********
C Steve Gibbons Fri Dec 15 10:22:21 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C There are NH3 temperature harmonics in the time-stepping solution  C
C vector. They are indexed by the integer arrays ML3 and MM3.        C
C Each are represented at NR radial grid points and the element for  C
C node ir and harmonic ih is stored in SV3( ( IH - 1 )*NR + IR ).    C
C The radial spacings are given by the array XARR.                   C
C                                                                    C
C Now the "foreign" temperature is stored in SVC with indexing       C
C by INARRC, MHTC, MHLC, MHMC and with grid nodes stored in the      C
C array XARRC. Now XARRC and XARR may be different but the first and C
C last points must really be the same! So, we have a variable called C
C DTOL which checks how close the inner and outer values are         C
C between the two arrays. If ( DABS( RI - RIC ).GT.DTOL   ) or       C
C ( DABS( RO - ROC ).GT.DTOL   ) then the program will abort.        C
C                                                                    C
C The routine uses the subroutine SVRINT to interpolate the 0th (F), C
C 1st (DF) and 2nd (DDF) derivatives of SVC and put into the         C
C appropriate positions of the arrays DV3C, SV3C and DSV3C           C
C respectively                                                       C
C                                                                    C
C SCAL*CD*DL( L , R , F , DF , DDF ),                                C
C SCAL*F                                                             C
C SCAL*DF.                                                           C
C                                                                    C
C DV3C, SV3C and DSV3C are blanked upon entry.                       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH3       : Number of temperature harmonics in main soln.      C
C     NR        : Number of radial grid nodes in main soln.          C
C     ML3       : Dim ( NH3 ). Sph. harm. degrees of temp. harms.    C
C     MM3       : Dim ( NH3 ). Sph. harm. orders of temp. harms.     C
C                                                                    C
C     INARRC    : Indexing array (dim 3) for foreign temperature.    C
C                   INARRC( 1 ) = iformfc                            C
C                   INARRC( 2 ) = NRC                                C
C                   INARRC( 3 ) = NHC                                C
C     MHTC      : Type array (dim NHC) for foreign temperature.      C
C     MHLC      : l array (dim NHC) for foreign temperature.         C
C     MHMC      : m array (dim NHC) for foreign temperature.         C
C                                                                    C
C     NNDS      : Number of nodes for interpolation. (atleast 3 )    C
C     IWORK     : Dimension ( NNDS ). Work array.                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim (NR). Radial grid nodes of main solution.      C
C     XARRC     : Dim (NRC). Radial grid nodes of foreign temp.      C
C     SVC       : Dim (NRC*NHC). foreign temp. solution vector.      C
C                                                                    C
C     CD        : Coefficient of diffusion for temperature.          C
C     SCAL      : Multiplier for inhomogeneous part of the temp.     C
C                                                                    C
C     DTOL      : Tolerance allowed for difference of inner and      C
C                 outer boundaries between the two spacings arrays.  C
C                                                                    C
C     DV3C      : Diffusion forcing term from inhomog. temp (NR*NH3) C
C     SV3C      : Zero^th derivative interpolation      Dim (NR*NH3) C
C     DSV3C     : First derivative interpolation      Dim (NR*NH3)   C
C                                                                    C
C     WORK1     : Dimension ( NNDS ). Work array.                    C
C     WORK2     : Dimension ( NNDS ). Work array.                    C
C     WORKM     : Dimension ( NNDS, NNDS ). Work array.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITVAAF( NH3, NR, ML3, MM3, INARRC, MHTC, MHLC, MHMC,
     1                  NNDS, IWORK, XARR, XARRC, SVC, CD, SCAL,
     2                  DTOL, DV3C, SV3C, DSV3C, WORK1, WORK2, WORKM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NH3, NR, ML3( NH3 ), MM3( NH3 ), INARRC( 3 ),
     1                 MHTC( * ), MHLC( * ), MHMC( * ), NNDS,
     2                 IWORK( NNDS )
      DOUBLE PRECISION XARR( NR ), XARRC( * ), SVC( * ), CD, SCAL,
     1                 DTOL, DV3C( NR*NH3 ), SV3C( NR*NH3 ),
     2                 DSV3C( NR*NH3 ), WORK1( NNDS ), WORK2( NNDS ),
     3                 WORKM( NNDS, NNDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IR, IH, IHC, IND, NRC, NHC, IOP, IHARM, L, N
      DOUBLE PRECISION RAD, RI, RIC, RO, ROC, ZERO, F, DF, DDF, DL
      EXTERNAL         DL
      PARAMETER      ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Zero arrays.
C
      IOP     = 0
      N       = NR*NH3
      CALL VECOP(  DV3C, ZERO, N, IOP )
      CALL VECOP(  SV3C, ZERO, N, IOP )
      CALL VECOP( DSV3C, ZERO, N, IOP )
C
      NRC     = INARRC( 2 )
      NHC     = INARRC( 3 )
C
      RI      = XARR(   1 )
      RO      = XARR(  NR )
C
      RIC     = XARRC(   1 )
      ROC     = XARRC( NRC )
C
      IF ( DABS( RI-RIC ).GT.DTOL  .OR.
     1     DABS( RO-ROC ).GT.DTOL      ) THEN
        PRINT *,' Subroutine ITVAAF.'
        PRINT *,' RI = ', RI,' RIC = ', RIC
        PRINT *,' RO = ', RO,' ROC = ', ROC
        PRINT *,' Incompatible radial functions.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO IHC = 1, NHC
        IF ( MHTC( IHC ).NE.3 ) GOTO 50
C       .
C       . OK so IHC is a temperature harmonic -
C       . so let's try to find the corresponding harmonic.
C       .
        L     = MHLC( IHC )
        IHARM = 0
        DO IH = 1, NH3
          IF (     ML3( IH ).EQ.MHLC( IHC )   .AND.
     1             MM3( IH ).EQ.MHMC( IHC )        ) IHARM = IH
        ENDDO
        IF ( IHARM.EQ.0 ) THEN
          PRINT *,' Subroutine ITVAAF.'
          PRINT *,' Solution vector contains no harmonic, ih,'
          PRINT *,' with ML3( ih ) = ', MHLC( IHC )
          PRINT *,' and MM3( ih ) = ', MHMC( IHC )
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
        DO IR = 2, NR - 1
          IND    = ( IHARM - 1 )*NR + IR
          RAD    = XARR( IR )
          CALL SVRINT( RAD, SVC, XARRC, INARRC, IHC, NNDS, WORK1,
     1                   IWORK, WORK2, WORKM )
C         .
          F      = WORK1( 1 )
          DF     = WORK1( 2 )
          DDF    = WORK1( 3 )
C         .
          DV3C ( IND ) = SCAL*CD*DL( L , RAD , F , DF , DDF )
          SV3C ( IND ) = SCAL*F
          DSV3C( IND ) = SCAL*DF
C         .
        ENDDO
C       .
 50   CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************
