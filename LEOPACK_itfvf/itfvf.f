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
C subroutine Spherical Harmonic Coef. Array Normalised Copy **********
C            -         -        -     -     -          -    **********
C Steve Gibbons Thu Mar 16 17:51:14 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Let the function g( \theta, \phi ) be expressed as the sum of      C
C Schmidt normalised spherical harmonics:                            C
C                                                                    C
C  g = \sum_{l,m} [  c_{l,mc} P_l^m( cos theta ) cos (m phi)         C
C               + c_{l,ms} P_l^m( cos theta ) sin (m phi)            C
C                                                                    C
C with the coefficients (ordered by the function INDSHC) given in    C
C the array SHC.                                                     C
C                                                                    C
C If all coefficients are zero, then SHCANC returns SHCN as a zero   C
C array. Otherwise, SHCN returns the coefficients scaled such that   C
C the spherical surface integral of g^2 is equal to VALN.            C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SHC       : Dim ( LH*(LH+2) ). Input array coefficients.       C
C     SHCN      : Dim ( LH*(LH+2) ). Output array coefficients.      C
C     VALN      : Value for normalisation.                           C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHCANC( LH, SHC, SHCN, VALN )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH
      DOUBLE PRECISION SHC( * ), SHCN( * ), VALN
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L, M, ICS, I, NH
      DOUBLE PRECISION DLOW, ZERO, PI, RNORM, COEF, FAC
      PARAMETER ( DLOW = 1.0d-8, ZERO = 0.0d0,
     1            PI = 3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NH = LH*(LH+2)
C
      RNORM = ZERO
      DO I = 1, NH
        CALL LMFIND( I, L, M, ICS )
        COEF = SHC( I )
        FAC  = 2.0d0*DBLE(L) + 1.0d0
        RNORM = RNORM + COEF*COEF*4.0d0*PI/FAC
      ENDDO
C
      IF ( RNORM.LT.DLOW ) THEN
        FAC = 0.0d0
      ELSE
        FAC = VALN/SQRT( RNORM )
      ENDIF
C
      DO I = 1, NH
        SHCN( I ) = SHC( I )*FAC
      ENDDO
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Inhomogeneous Temperature Harmonic Coefficient Allocate *
C            -             -           -        -           -        *
C Steve Gibbons Wed Jan 19 08:46:34 GMT 2000              Routine    C
C____________________________________________________________________C
C                                                                    C
C ITHCAR receives arrays defining sets of harmonics; these are       C
C                                                                    C
C  MHT( ih ) = 1 for a poloidal velocity harmonic                    C
C  MHT( ih ) = 2 for a toroidal velocity harmonic                    C
C  MHT( ih ) = 3 for a temperature harmonic                          C
C  MHT( ih ) = 4 for a poloidal magnetic field harmonic              C
C  MHT( ih ) = 5 for a toroidal magnetic field harmonic              C
C                                                                    C
C  MHL( ih ) = l, degree of spherical harmonic                       C
C  MHM( ih ) = m (order) for cos ( m phi ) dependence or             C
C              -m for sin ( m phi ) dependence.                      C
C                                                                    C
C  A fourth array MHI is assigned values by ITHCAR: If MHT( ih ) = 3 C
C then ITHCAR will look at the arrays HMIB and HMOB (harmonic map    C
C for inner/outer boundary) and the flags KIB and KOB and decide     C
C what value of temperature/temp. gradient that harmonic radial      C
C function should achieve at that boundary.                          C
C                                                                    C
C MHIB and MHOB are indexed by INDSHC.                               C
C                                                                    C
C These values are supplied to the routine ITFCF as VALIB and VALOB  C
C and along with RI and RO will be used to provide the coefficients  C
C CA, CB and CC which define the inhomogeneous temperature function  C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C                                                                    C
C The coefficients CA, CB and CC are stored in the array             C
C CAFIT (coefficient array for inhomogeneous temperature) which has  C
C dimensions ( 3, NITHMX ) where NITHMX limits the possible number   C
C of temperature harmonics.                                          C
C                                                                    C
C NITH is the number of inhomogeneous temperature harmonics which    C
C have already been assigned coefficients by ITHCAR.                 C
C                                                                    C
C For example if IH is a harmonic radial function with MHT( IH ) = 3 C
C and ITFCF calculates that f(r) should take the coefficients        C
C CA, CB and CC; then                                                C
C                                                                    C
C  CA is stored in CAFIT( 1, IITH )                                  C
C  CB is stored in CAFIT( 2, IITH )                                  C
C  CC is stored in CAFIT( 3, IITH )                                  C
C                                                                    C
C and the index IITH is stored in MHI( IH )                          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     KIB      : 1 if T( r ) is to be fixed at the inner boundary.   C
C                2 if dT/dr(r) is to be fixed at the inner boundary. C
C                                                                    C
C     KOB      : 1 if T( r ) is to be fixed at the outer boundary.   C
C                2 if dT/dr(r) is to be fixed at the outer boundary. C
C                                                                    C
C     NITH     : Number of inhomogeneous temperature harmonics       C
C                 with coeff.s stored in CAFIT.                      C
C                                                                    C
C     NITHMX   : Limit on NITH.                                      C
C                                                                    C
C     NH       : Number of harmonics                                 C
C     MHT      : Dim( * ). Harmonic type - see above.                C
C     MHL      : Dim( * ). Harmonic degree - see above.              C
C     MHM      : Dim( * ). Harmonic order  - see above.              C
C     MHI      : Dim( * ). 2nd index of CAFIT where coeffs are storedC
C                                                                    C
C     LH       : Maximum spherical harmonic degree.                  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RI       : Radius of the inner boundary.                       C
C     RO       : Radius of the outer boundary.                       C
C                                                                    C
C     HMIB     : Harmonic map for inner boundary. Dim( LH*(LH+2) )   C
C     HMIBMT   : Harmonic map for inner boundary monopole term       C
C     HMOB     : Harmonic map for outer boundary. Dim( LH*(LH+2) )   C
C     HMOBMT   : Harmonic map for outer boundary monopole term       C
C                                                                    C
C                 MHIB and MHOB can refer to either temperature      C
C                 or temperature gradient depending upon the         C
C                 values of KIB and KOB.                             C
C                                                                    C
C                 The arrays are indexed by INDSHC.                  C
C                                                                    C
C     CAFIT    : Coeff arr. for inhomog. temp. Dim( 3, NITHMX ).     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITHCAR( KIB, KOB, NITH, NITHMX, NH, MHT, MHL, MHM,
     1                   MHI, LH, RI, RO, HMIB, HMIBMT, HMOB, HMOBMT,
     2                   CAFIT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KIB, KOB, NITH, NITHMX, NH, MHT( * ), MHL( * ),
     1        MHM( * ), MHI( * ), LH
      DOUBLE PRECISION RI, RO, HMIB( LH*(LH+2) ), HMIBMT,
     1                  HMOB( LH*(LH+2) ), HMOBMT, CAFIT( 3, NITHMX )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IHARM, L, M, ICS, IH, INDSHC, NIT
      DOUBLE PRECISION VALIB, VALOB, CA, CB, CC, DLOW
      LOGICAL OK
      PARAMETER ( DLOW = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      IF ( NITH.LT.0 ) THEN
        PRINT *,' Subroutine ITHCAR.'
        PRINT *,' NITH = ', NITH,' on entry.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      OK = .TRUE.
C     .
      DO IH = 1, NH
        IF ( MHT( IH ).NE.3 ) THEN
          MHI( IH ) = -1
          GOTO 50
        ENDIF
C       .
        L = MHL( IH )
        IF ( MHM( IH ).LT.0 ) THEN
          M   = -MHM( IH )
          ICS = 2
        ELSE
          M   = MHM( IH )
          ICS = 1
        ENDIF
        IHARM = INDSHC( L, M, ICS )
C       .
        IF ( L.GT.LH ) THEN
          PRINT *,' Subroutine ITHCAR.'
          PRINT *,' Harmonic ',IH,' has L = ', L
          PRINT *,' Maximum degree = ', LH
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
        IF ( IHARM.EQ.0 ) THEN
          VALIB = HMIBMT
          VALOB = HMOBMT
        ELSE
          VALIB = HMIB( IHARM )
          VALOB = HMOB( IHARM )
        ENDIF
C       .
        CALL ITFCF( KIB, KOB, RI, RO, VALIB, VALOB, CA, CB, CC )
C       .
C       . OK - we have calculated values
C       . for CA, CB and CC, so let's loop around
C       . CAFIT to see if we have already stored
C       . such values.
C       .
        DO NIT = 1, NITH
          IF ( DABS( CA - CAFIT( 1, NIT ) ).LT.DLOW .AND.
     1         DABS( CB - CAFIT( 2, NIT ) ).LT.DLOW .AND.
     2         DABS( CC - CAFIT( 3, NIT ) ).LT.DLOW   ) THEN
            MHI( IH ) = NIT
            GOTO 50
          ENDIF
        ENDDO
C       .
        CALL CNTRIC( NITH, NITHMX, OK )
C       .
        IF ( OK ) THEN
          CAFIT( 1, NITH ) = CA
          CAFIT( 2, NITH ) = CB
          CAFIT( 3, NITH ) = CC
          MHI( IH ) = NITH
        ENDIF
C       .
 50     CONTINUE
      ENDDO
C     .
      IF ( OK ) RETURN
      PRINT *,' Subroutine ITHCAR.'
      PRINT *,' Number of requested coefficients = ', NITH
      PRINT *,' NITHMX = ', NITHMX
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
C*********************************************************************
C subroutine Inhomogeneous Temperature Function Add ******************
C            -             -           -        -   ******************
C Steve Gibbons Mon Jan 17 11:25:35 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C The function f( r ) for ri .le. r .le. ro is defined by            C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C                                                                    C
C Then clearly,    f( ri )  = CB 2 ( ri-ro )/pi + CC                 C
C                  f'( ri ) = CA 0.5 pi/(ro-ri)                      C
C                                                                    C
C                  f( ro )  = CA                + CC                 C
C                  f'( ro ) = CB                                     C
C                                                                    C
C For the coefficients CA, CB and CC (see ITFCF) ITFA will return    C
C derivatives 0 to IHD of f( r ) with f[nd]( r ) in DERV( nd + 1 )   C
C                                                                    C
C Currently IHD may be no larger than 4 although the routine could   C
C easily be modified for higher derivatives if the need arose.       C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IHD      : The number of the highest derivative requested.     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD      : Current radius.                                     C
C     RI       : Radius of the inner boundary.                       C
C     RO       : Radius of the outer boundary.                       C
C     CA       : Coefficient of term in f(r). See above.             C
C     CB       : Coefficient of term in f(r). See above.             C
C     CC       : Coefficient of term in f(r). See above.             C
C     DERV     : Array of length atleast ( IHD + 1 ).                C
C                 DERV( nd + 1 ) returned with nd^{th} deriv. of f.  C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITFA( RAD, RI, RO, CA, CB, CC, DERV, IHD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IHD
      DOUBLE PRECISION RAD, RI, RO, CA, CB, CC, DERV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION PI, ROMRI, RMRI, LOW, FAC, TERMA, TERMB,
     1                 DOPRND
      PARAMETER ( PI=3.14159265358979312D0, LOW = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      IF (      DABS( CA ).LT.LOW      .AND.
     1          DABS( CB ).LT.LOW      .AND.
     2          DABS( CC ).LT.LOW     )      RETURN
C     .
      RMRI  = RAD - RI
      ROMRI = RO - RI
      IF ( DABS( ROMRI ).LT.LOW ) THEN
        PRINT *,' Subroutine ITFA.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check on bounds for IHD
C     .
      IF ( IHD.LT.0 .OR. IHD.GT.4 ) THEN
        PRINT *,' Subroutine ITFA.'
        PRINT *,' IHD = ', IHD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      FAC    = 0.5d0*PI/ROMRI
      DOPRND = FAC*RMRI
C
      TERMA  = CA*DSIN( DOPRND )
      TERMB  = CB*DCOS( DOPRND )*(-1.0d0)/FAC
      DERV( 1 ) = DERV( 1 ) + TERMA + TERMB + CC
      IF ( IHD.EQ.0 ) RETURN
C     .
      TERMA  = CA*FAC*DCOS( DOPRND )
      TERMB  = CB*DSIN( DOPRND )
      DERV( 2 ) = DERV( 2 ) + TERMA + TERMB
      IF ( IHD.EQ.1 ) RETURN
C     .
      TERMA  = (-1.0d0)*CA*FAC*FAC*DSIN( DOPRND )
      TERMB  = CB*FAC*DCOS( DOPRND )
      DERV( 3 ) = DERV( 3 ) + TERMA + TERMB
      IF ( IHD.EQ.2 ) RETURN
C     .
      TERMA  = (-1.0d0)*CA*FAC*FAC*FAC*DCOS( DOPRND )
      TERMB  = (-1.0d0)*CB*FAC*FAC*DSIN( DOPRND )
      DERV( 4 ) = DERV( 4 ) + TERMA + TERMB
      IF ( IHD.EQ.3 ) RETURN
C     .
      TERMA  = CA*FAC*FAC*FAC*FAC*DSIN( DOPRND )
      TERMB  = (-1.0d0)*CB*FAC*FAC*FAC*DCOS( DOPRND )
      DERV( 5 ) = DERV( 5 ) + TERMA + TERMB
      IF ( IHD.EQ.4 ) RETURN
C     .
      STOP
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
C********************************************************************
C subroutine L and M FIND *******************************************
C            -     - ---- *******************************************
C Steve Gibbons 12.4.97                                             C
C Modified 12.6.97 to have IT = 1 for cos and 2 for sine            C
C___________________________________________________________________C
C Given a number of harmonic, N; LMFIND will return the level, L,   C
C the order, M, and IT - which is equal to 2 for sin                C
C spherical harmonics and 1 for cosine ones.                        C
C All the above are integers - no point in a variable list .....    C
C     N = L*L for M = 0, IT = 1                                     C
C     N = L*L + 2*M - 1 for non-zero M and IT = 1                   C
C     N = L*L + 2*M for non-zero M and IT = 2                       C
C If N = 0, we have a monopole: L = 0, M = 0 and IT = 1.            C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE LMFIND ( N, L, M, IT)
      IMPLICIT NONE
C___________________________________________________________________C
C Variable Declarations - Parameters ...............................C
      INTEGER N,L,M,IT
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      INTEGER LL,IDIFF,N2,ITWO
      PARAMETER (ITWO=2)
C___________________________________________________________________C
C First put N into N2 so that N is not altered
      N2=N
      IF ( N2.EQ.0 ) THEN
         L=0
         M=0
         IT=1
         RETURN
      ENDIF
      IF ( N2.LT.1 ) THEN 
         WRITE (6,989)
         STOP
 989  FORMAT (' Subroutine LMFIND. N<1 - Program stopped.')
      ENDIF
      IF ( N2.EQ.1 ) THEN
         L=1
         M=0
         IT=1
         RETURN
      ENDIF
      L=1
 500  CONTINUE
      LL=L*L
      IDIFF = N2 - LL
      IF ( IDIFF.GT.0 ) THEN
         L=L+1
         GOTO 500
      ENDIF
      IF ( IDIFF.EQ.0 ) THEN
         M=0
         IT=1
         RETURN
      ENDIF
      L=L-1
      LL=L*L
      N2=N2-LL
C so we know now that N2 is equal to either 2*M or 2*M-1
C corresponding to IT=1 and IT=2 respectively
      IDIFF = MOD ( N2, ITWO )
      IF ( IDIFF.EQ.1) THEN
         IT = 1
         M = (N2+1)/ITWO
      ELSE
         IT = 2
         M = N2/ITWO
      ENDIF
      RETURN
      END
C********************************************************************
C*********************************************************************
C subroutine Inhomogeneous Temperature Function Coefficient Find *****
C            -             -           -        -           -    *****
C Steve Gibbons Mon Jan 17 11:25:35 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C The function f( r ) for ri .le. r .le. ro is defined by            C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C                                                                    C
C Then clearly,    f( ri )  = CB 2 ( ri-ro )/pi + CC                 C
C                  f'( ri ) = CA 0.5 pi/(ro-ri)                      C
C                                                                    C
C                  f( ro )  = CA                + CC                 C
C                  f'( ro ) = CB                                     C
C                                                                    C
C ITFCF chooses the coefficients CA, CB and CC such that f will      C
C have the desired properties at ri and ro.                          C
C                                                                    C
C Set KIB (KOB) to 1 to fix the temperature at the inner (outer) bnd C
C Set KIB (KOB) to 2 to fix the heat flux at the inner (outer) bnd.  C
C                                                                    C
C VALIB (VALOB) is the actual value which f( r ) or f'( r ) is       C
C to attain at the inner (outer) boundary, depending upon the        C
C setting of KIB (KOB).                                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     KIB      : 1 if T( ri ) is to be set to the value VALIB.       C
C                2 if dT/dr( ri ) is to be set to the value VALIB.   C
C                                                                    C
C     KOB      : 1 if T( ro ) is to be set to the value VALOB.       C
C                2 if dT/dr( ro ) is to be set to the value VALOB.   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RI       : Radius of the inner boundary.                       C
C     RO       : Radius of the outer boundary.                       C
C                                                                    C
C     VALIB    : Value for temp./temp. gradient at inner boundary.   C
C     VALOB    : Value for temp./temp. gradient at outer boundary.   C
C                                                                    C
C     CA       : Coefficient of term in f(r). See above.             C
C     CB       : Coefficient of term in f(r). See above.             C
C     CC       : Coefficient of term in f(r). See above.             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITFCF( KIB, KOB, RI, RO, VALIB, VALOB, CA, CB, CC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KIB, KOB
      DOUBLE PRECISION RI, RO, VALIB, VALOB, CA, CB, CC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION PI, ROMRI, LOW
      PARAMETER ( PI=3.14159265358979312D0, LOW = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      ROMRI = RO - RI
      IF ( DABS( ROMRI ).LT.LOW ) THEN
        PRINT *,' Subroutine ITFCF.'
        PRINT *,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . First deal with the case of fixed temperature
C     . at both inner and outer boundaries
C     .
      IF ( KIB.EQ.1 .AND. KOB.EQ.1 ) THEN
        CA = VALOB
        CB = VALIB*PI*(-0.5d0)*ROMRI
        CC = 0.0d0
        RETURN
      ENDIF
C     .
C     . Fixed temperature at inner boundary
C     . Fixed flux at outer boundary
C     .
      IF ( KIB.EQ.1 .AND. KOB.EQ.2 ) THEN
        CA = 0.0d0
        CB = VALOB
        CC = VALIB - CB*2.0d0*(RI-RO)/PI
        RETURN
      ENDIF
C     .
C     . Fixed flux at inner boundary
C     . Fixed temperature at outer boundary
C     .
      IF ( KIB.EQ.2 .AND. KOB.EQ.1 ) THEN
        CA = VALIB*ROMRI*2.0d0/PI
        CB = 0.0d0
        CC = VALOB - CA
        RETURN
      ENDIF
C     .
C     . Fixed flux at inner and outer boundaries
C     .
      IF ( KIB.EQ.2 .AND. KOB.EQ.2 ) THEN
        CA = VALIB*ROMRI*2.0d0/PI
        CB = VALOB
        CC = 0.0d0
        RETURN
      ENDIF
C     .
      PRINT *,' Subroutine ITFCF.'
      PRINT *,' KIB = ', KIB,' KOB = ',KOB
      PRINT *,' Illegal values for KIB and KOB.'
      PRINT *,' Program aborted.'
C     .
      STOP
      END
C*********************************************************************

C*********************************************************************
C subroutine CouNTeR Increment and Check *****************************
C            -  -- - -             -     *****************************
C Steve Gibbons Mon Jan 10 10:16:42 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Adds 1 to the number IC.                                           C
C If the resulting number is greater than ICMAX then the logical     C
C flag OK is set to .FALSE.                                          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IC       : Value of counter variable.                          C
C     ICMAX    : Maximum permitted value of IC.                      C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     OK       : Unaltered if IC.le.ICMAX.                           C
C                Set to .FALSE. if IC.gt.ICMAX.                      C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CNTRIC( IC, ICMAX, OK )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IC, ICMAX
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      IC = IC + 1
      IF ( IC.GT.ICMAX ) OK = .FALSE.
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
C*********************************************************************
C function INDex for Spherical Harmonic Coefficient ******************
C          ---       -         -        -           ******************
C Steve Gibbons 25.4.97                                              C
C____________________________________________________________________C
C Inputs are integers, L and M are obvious ICS = 1 for a cosine harm C
C and ICS = 2 for a sine harmonic.                                   C
C____________________________________________________________________C
      FUNCTION INDSHC ( L, M, ICS )
      IMPLICIT NONE
      INTEGER INDSHC, L, M, ICS
C 
      IF ( L.EQ.0 .AND. M.EQ.0 .AND. ICS.EQ.1 ) THEN
         INDSHC = 0
         RETURN
      ENDIF
C 
      IF ( L.LT.1 ) THEN
        PRINT *,' Function INDSHC. L = ', L
        STOP
      ENDIF
C 
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' Function INDSHC. M invalid. Program Aborted.'
         PRINT *,' L = ', L
         PRINT *,' M = ', M
         STOP
      ENDIF
      IF ( ICS.NE.1 .AND. ICS.NE.2 ) THEN
         PRINT *,' Function INDSHC. ICS = ', ICS
         PRINT *,' ICS must be 1 or 2.'
         PRINT *,' Program Aborted.'
         STOP
      ENDIF
      IF ( M.EQ.0 .AND. ICS.EQ.1 ) THEN
         INDSHC = L*L
         RETURN
      ENDIF
      IF ( M.EQ.0 .AND. ICS.NE.1 ) THEN
         PRINT *,' Function INDSHC. M = ', M
         PRINT *,' ICS = ', ICS,' Program aborted.'
         STOP
      ENDIF
      INDSHC = L*L + 2*M - 2 + ICS
C
      RETURN
      END
C*********************************************************************

