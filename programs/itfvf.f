C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C             Inhomogeneous Temperature Function Vector Form         C
C             -             -           -        -      -            C
C                                                                    C
C Mon Mar 27 12:19:14 GMT 2000                                       C
C                                                                    C
C Inputs a set of boundary conditions for the temperature            C
C (i.e. fixed temperature/fixed heat flux at inner/outer boundary)   C
C and reads in sets of spherical harmonic coefficients for           C
C the inhomogeneity of the inner and outer boundary functions.       C
C                                                                    C
C A solution vector is then written to file with a temperature       C
C which satisfies exactly this requested b.c.                        C
C                                                                    C
C These coefficients must be supplied in a separate file.            C
C Each line of this file must begin either with an asterisk (which   C
C comments out the remainder of the line) with the characters IB     C
C (followed by L, M, ICS, COEF - see INDSHC) for a harmonic for the  C
C inner boundary and similarly for the outer boundary.               C
C                                                                    C
C Finally, a .xarr, .ints and .vecs file is written out with         C
C the filename stem ROOT.                                            C
C                                                                    C
C*********************************************************************
      PROGRAM itfvf
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, NHMAX, ISVMAX, LHMAX, LHLH2M, NDCS, NITHMX
      PARAMETER ( NRMAX = 100, LHMAX = 64, LHLH2M = LHMAX*(LHMAX+2),
     1            NHMAX = LHLH2M+1, ISVMAX = NHMAX*NRMAX, NDCS = 1,
     2            NITHMX = NHMAX )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ), MHP( NHMAX ),
     1        MHI( NHMAX ), MHIBC( NDCS ), MHOBC( NDCS )
      DOUBLE PRECISION XARR( NRMAX ), VEC( ISVMAX ),
     1                 CAFIT( 3, NITHMX ), SHCI( LHLH2M ), SHCIM,
     2                 SHCO( LHLH2M ), SHCOM, SHCI2( LHLH2M ),
     3                 SHCO2( LHLH2M )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NH, IREAD, L, M, ICS, MM, INDSHC, ILEN, I, NR, ISP,
     1        IFORMF, INARR( 3 ), IND, ITHEBC, NITH, KIB, KOB,
     2        IITH, IH, LU, INDFUN, IOP, IFORM, IHD, LH
      DOUBLE PRECISION ZERO, COEF, RI, RO, DERV( 1 ), CAK, CBK, CCK,
     1                 RAD, EPSI, EPSO
      CHARACTER *(2)   BCH
      CHARACTER *(80)  FNAME, ROOT
      CHARACTER *(200) LINE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( ZERO = 0.0d0, IREAD = 1, IOP = 0 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      LU = 99
C
      DO I = 1, 200
        LINE(I:I)  = ' '
      ENDDO
C
      DO I = 1, 80
        FNAME(I:I) = ' '
        ROOT(I:I)  = ' '
      ENDDO
C
 80   FORMAT(A)
C
      PRINT *,' Please enter root for output files.'
 21   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 21
      DO I = 1, 80
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 49
        ENDIF
      ENDDO
 49   CONTINUE
      ROOT = LINE(1:ILEN)
C
      SHCIM = ZERO
      SHCOM = ZERO
      CALL VECOP( SHCI, ZERO, LHLH2M, IOP )
      CALL VECOP( SHCO, ZERO, LHLH2M, IOP )
C
      MHIBC( 1 ) = 1
      MHOBC( 1 ) = 1
C
C Enter name of file containing boundary coeff.s
C
      PRINT *,' Please enter filename of coefficients file.'
 22   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 22
      DO I = 1, 80
        IF ( LINE(I:I).EQ.' ' ) THEN
          NH = I-1
          GOTO 48
        ENDIF
      ENDDO
 48   CONTINUE
      FNAME = LINE(1:NH)
      NH = 0
C
C Now let's read in boundary coefficients
C
      CALL FOPEN( LU, FNAME, IREAD )
 92   CONTINUE
      READ ( LU, 80, END = 93 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 92
C
      BCH = LINE(1:2)
      READ ( LINE(3:80), * ) L, M, ICS, COEF
      IF ( BCH.NE.'IB' .AND. BCH.NE.'OB' ) THEN
        PRINT *,' Start of boundary coefficient line'
        PRINT *,' must either be IB (inner) or '
        PRINT *,' OB (outer) boundary.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ICS.EQ.2 ) THEN
        MM = -M
      ELSE
        MM = M
      ENDIF
C
C See if we have read in this harmonic before
C
      IF ( NH.EQ.0 ) GOTO 76
      DO I = 1, NH
        IF (  MHL( I ).EQ.L .AND. MHM( I ).EQ.MM ) GOTO 77
      ENDDO
 76   CONTINUE
      NH = NH + 1
      IF ( NH.GT.NHMAX ) THEN
        PRINT *,' Maximum NH = ', NHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      MHT( NH ) = 3
      MHL( NH ) = L
      MHM( NH ) = MM
      MHP( NH ) = 1
 77   CONTINUE
C
      IF ( L.GT.LHMAX ) THEN
        PRINT *,' Boundary coefficient with L = ', L
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IND = INDSHC( L, M, ICS )
      IF ( BCH.EQ.'IB' ) THEN
        IF ( IND.EQ.0 ) THEN
          SHCIM = COEF
        ELSE
          SHCI( IND ) = COEF
        ENDIF
      ELSE
        IF ( IND.EQ.0 ) THEN
          SHCOM = COEF
        ELSE
          SHCO( IND ) = COEF
        ENDIF
      ENDIF
C
      GOTO 92
 93   CONTINUE
      CALL FCLOSE( LU, FNAME, 'Error' )
C
C Enter Format data
C
      PRINT *,' Enter NR, ISP, IFORMF, RI, RO, ITHEBC '
      PRINT *,' ----------------------------- '
      PRINT *,' nr - number of radial grid nodes.'
      PRINT *,' isp = 1 --> uniform nodes.'
      PRINT *,' isp = 2 --> Chebshev nodes.'
      PRINT *,' iformf = 3 --> ind = (ir-1)*nh + ih '
      PRINT *,' iformf = 4 --> ind = (ih-1)*nr + ir '
      PRINT *,' ri = inner radius value '
      PRINT *,' ro = outer radius value '
      PRINT *,' ithebc = 1 --> fixed temp inner and outer '
      PRINT *,' ithebc = 2 --> fixed temp inner, fixed flux outer '
      PRINT *,' ithebc = 3 --> fixed flux inner, fixed temp outer '
      PRINT *,' ----------------------------- '
      READ ( 5, * ) NR, ISP, IFORMF, RI, RO, ITHEBC
C
      IF ( NR.LT.1 .OR. NR.GT.NRMAX ) THEN
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
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
        PRINT *,' IFORMF = ', IFORMF
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
      PRINT *,' Enter EPSI, EPSO : Scaling of g^2 '
      READ (5, * ) EPSI, EPSO
C
C Fill radial value arrays
C
      IF ( ISP.EQ.1 ) THEN
        CALL ESNAAS( NR, XARR, RI, RO )
      ELSE
        CALL ZCPAAS( NR, XARR, RI, RO )
      ENDIF
C
      INARR( 1 ) = IFORMF
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C
      LH = 0
      DO I = 1, NH
        IF ( MHL( I ).GT.LH ) LH = MHL( I )
      ENDDO
C
      IF ( ITHEBC.EQ.1 ) THEN
        KIB        = 1
        KOB        = 1
      ENDIF
C
      IF ( ITHEBC.EQ.2 ) THEN
        KIB        = 1
        KOB        = 2
      ENDIF
C
      IF ( ITHEBC.EQ.3 ) THEN
        KIB        = 2
        KOB        = 1
      ENDIF
C
      CALL SHCANC( LH, SHCI, SHCI2, EPSI )
      CALL SHCANC( LH, SHCO, SHCO2, EPSO )
      NITH = 0
      CALL ITHCAR( KIB, KOB, NITH, NITHMX, NH, MHT, MHL, MHM,
     1             MHI, LH, RI, RO, SHCI2, SHCIM, SHCO2, SHCOM,
     2             CAFIT )
C     .
      IHD = 0
      DO IH = 1, NH
        IITH   = MHI( IH )
        CAK    = CAFIT( 1, IITH )
        CBK    = CAFIT( 2, IITH )
        CCK    = CAFIT( 3, IITH )
        DO I = 1, NR
          RAD = XARR( I )
C         .
          IND = INDFUN( I, IH, INARR )
          DERV( 1 ) = 0.0d0
C         .
          CALL ITFA( RAD, RI, RO, CAK, CBK, CCK, DERV, IHD )
C         .
          VEC( IND ) = DERV( 1 )
C         .
        ENDDO
      ENDDO
C
C Write out the harmonic integers file
C
      FNAME(1:ILEN) = ROOT(1:ILEN)
      FNAME(ILEN+1:ILEN+5) = '.ints'
      FNAME = FNAME(1:ILEN+5)
      CALL HMFWT( NH, MHT, MHL, MHM, MHP, NDCS, MHIBC, MHOBC,
     1            LU, FNAME )
C
C Write out the eigenvectors
C
      FNAME(ILEN+1:ILEN+5) = '.vecs'
      FNAME = FNAME(1:ILEN+5)
      IFORM = 1
      CALL SVFWT( INARR, LU, IFORM, VEC, FNAME )
C
C Write out radial node data
C
      FNAME(ILEN+1:ILEN+5) = '.xarr'
      FNAME = FNAME(1:ILEN+5)
      CALL XARRWT( NR, XARR, LU, FNAME, IFORM )
C     .
      STOP
      END
C*********************************************************************
