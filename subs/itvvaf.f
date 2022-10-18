C*********************************************************************
C subroutine Imposed Toroidal Velocity Vector Array Form *************
C            -       -        -        -      -     -    *************
C Steve Gibbons Wed Dec 20 08:51:58 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C There are NH2 toroidal velocity harm.s in the time-stepping soln.  C
C vector. They are indexed by the integer arrays ML2 and MM2.        C
C Each are represented at NR radial grid points and the element for  C
C node ir and harmonic ih is stored in SV2( ( IH - 1 )*NR + IR ).    C
C The radial spacings are given by the array XARR.                   C
C                                                                    C
C Now the "foreign" toroidal vel. is stored in SVC with indexing     C
C by INARRC, MHTC, MHLC, MHMC and with grid nodes stored in the      C
C array XARRC. Now XARRC and XARR may be different but the first and C
C last points must really be the same! So, we have a variable called C
C DTOL which checks how close the inner and outer values are         C
C between the two arrays. If ( DABS( RI - RIC ).GT.DTOL   ) or       C
C ( DABS( RO - ROC ).GT.DTOL   ) then the program will abort.        C
C                                                                    C
C The routine uses the subroutine SVRINT to interpolate the 0th (F)  C
C derivative of SVC and puts into the appropriate                    C
C positions of the array DIVT                                        C
C                                                                    C
C SCAL*DSQRT( DBLE( L*L + L ) )*(-1.0d0)*F                           C
C                                                                    C
C DIVT is blanked upon entry.                                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH2       : Number of tor. vel. harmonics in main soln.        C
C     NR        : Number of radial grid nodes in main soln.          C
C     ML2       : Dim ( NH2 ). Sph. harm. degrees of tor. v. harms.  C
C     MM2       : Dim ( NH2 ). Sph. harm. orders of tor. v. harms.   C
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
C     SCAL      : Multiplier for inhomogeneous part of the temp.     C
C                                                                    C
C     DTOL      : Tolerance allowed for difference of inner and      C
C                 outer boundaries between the two spacings arrays.  C
C                                                                    C
C     DIVT      : Toroidal array of imposed v_0. Dim (NR*NH2)        C
C                                                                    C
C     WORK1     : Dimension ( NNDS ). Work array.                    C
C     WORK2     : Dimension ( NNDS ). Work array.                    C
C     WORKM     : Dimension ( NNDS, NNDS ). Work array.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITVVAF( NH2, NR, ML2, MM2, INARRC, MHTC, MHLC, MHMC,
     1                   NNDS, IWORK, XARR, XARRC, SVC, SCAL, DTOL,
     2                   DIVT, WORK1, WORK2, WORKM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NH2, NR, ML2( NH2 ), MM2( NH2 ), INARRC( 3 ),
     1                 MHTC( * ), MHLC( * ), MHMC( * ), NNDS,
     2                 IWORK( NNDS )
      DOUBLE PRECISION XARR( NR ), XARRC( * ), SVC( * ), SCAL,
     1                 DTOL, DIVT( NR*NH2 ), WORK1( NNDS ),
     2                 WORK2( NNDS ), WORKM( NNDS, NNDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IR, IH, IHC, IND, NRC, NHC, IOP, IHARM, L, N
      DOUBLE PRECISION RAD, RI, RIC, RO, ROC, ZERO, F, DF
      PARAMETER      ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Zero arrays.
C
      IOP     = 0
      N       = NR*NH2
      CALL VECOP( DIVT, ZERO, N, IOP )
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
        PRINT *,' Subroutine ITVVAF.'
        PRINT *,' RI = ', RI,' RIC = ', RIC
        PRINT *,' RO = ', RO,' ROC = ', ROC
        PRINT *,' Incompatible radial functions.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO IHC = 1, NHC
        IF ( MHTC( IHC ).NE.2 ) GOTO 50
C       .
C       . OK so IHC is a toroidal velocity harmonic -
C       . so let's try to find the corresponding harmonic.
C       .
        L     = MHLC( IHC )
        IHARM = 0
        DO IH = 1, NH2
          IF (     ML2( IH ).EQ.MHLC( IHC )   .AND.
     1             MM2( IH ).EQ.MHMC( IHC )        ) IHARM = IH
        ENDDO
        IF ( IHARM.EQ.0 ) THEN
          PRINT *,' Subroutine ITVVAF.'
          PRINT *,' Solution vector contains no harmonic, ih,'
          PRINT *,' with ML2( ih ) = ', MHLC( IHC )
          PRINT *,' and MM2( ih ) = ', MHMC( IHC )
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
C         .
          DIVT( IND ) = SCAL*DSQRT( DBLE( L*L + L ) )*(-1.0d0)*F
C         .
        ENDDO
C       .
 50   CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************
